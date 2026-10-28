"""Multi-copy nuclear markers - rDNA (18S-ITS1-5.8S-ITS2-28S) and histones (H1, H2A, H2B, H3, H4): GetOrganelle
assembly of the tandem unit, MitoM annotation (GenBank .gbf), GetOrganelle vs SPAdes-MitoM comparison.

Why both: GetOrganelle walks the multi-copy unit in its graph and returns it whole (gene order, spacers); SPAdes
scatters the copies over many contigs, which MitoM merges into one consensus per gene. The same genome must give
the same genes, so differences point at assembly or annotation problems (MJ20260820041: the SPAdes histone cluster
is split into 84-134 nt contigs; GetOrganelle from H1-H4 seeds gave the cluster H3 whole).

Histones (MJ20260820041, CHY006): a label of H3 only kept just the H3 piece of the cluster (2.2 kb); H1-H4 cluster
copies as seed / label (resources/histone) gave all five genes.

Reads (user, 2026-10-28): multi-copy units need little data - round 1 uses at most 5 Gb; when the unit does not
close (no circular path), round 2 uses 10 Gb with round 1's result as the seed (the iteration of mt_iterate.py).
An assembly the results directory already has (rDNA of wgs_skim_pipeline) counts as round 1.

Subcommands (also used by the workflow rules, workflow/rules/nuclear.smk):
  assemble  GetOrganelle -F anonym of a multi-copy unit (--kind rdna | histone) from clean reads, in rounds
  annotate  `mitom nuc annotate` of a GetOrganelle result or SPAdes scaffolds -> <sample>.<kind>.gbf + MitoM tables
  compare   GetOrganelle vs SPAdes-MitoM per gene -> TSV
  batch     every sample of an existing results directory (LAYOUTS), resumable: reuses finished assemblies, redoes an
            annotation when MitoM changed (mitom.version)
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

RES = Path(__file__).resolve().parents[2] / "resources"
GENES = {"rdna": ["18S", "ITS1", "5.8S", "ITS2", "28S"], "histone": ["H1", "H2A", "H2B", "H3", "H4"]}
# bowtie2 options for GetOrganelle's seed reads: histone DNA of another family is too far from the bundled seed for
# end-to-end mapping (HNP1sub, a Pisione: 391 bytes of seed reads, "Too few seed reads", no result; local
# --very-sensitive-local: 153 kb of seed reads, H2A H2B H3 H4 complete on 9.7 kb); rDNA: GetOrganelle's default
BOWTIE2 = {"rdna": "--very-fast -t", "histone": "--very-sensitive-local -t"}
SEEDS = {"rdna": (RES / "rdna" / "seed.fasta", RES / "rdna" / "label.fasta"),
         "histone": (RES / "histone" / "seed.fasta", RES / "histone" / "label.fasta")}
COLS = ["sample", "marker", "gene", "getorganelle_nt", "getorganelle_partial", "getorganelle_class", "spades_nt",
        "spades_partial", "spades_class", "overlap_nt", "mismatches", "indel_nt", "identity_pct", "verdict"]


# ── assembly ───────────────────────────────────────────────────────────────
def result_fasta(d: Path):
    """GetOrganelle's final sequence: complete graph first, scaffolds only as fallback (as collect_fasta.py)."""
    fas = sorted(Path(d).glob("*.path_sequence.fasta")) if Path(d).is_dir() else []
    return next((f for f in fas if "scaffolds" not in f.name), fas[0] if fas else None)


def reads_per_file(r1, gb: float) -> int:
    """GetOrganelle --max-reads (reads per file) for gb Gb of paired reads; mean length of the first 2000 reads"""
    import gzip
    n = tot = 0
    with (gzip.open if str(r1).endswith(".gz") else open)(r1, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                n, tot = n + 1, tot + len(line.strip())
                if n == 2000:
                    break
    return int(gb * 1e9 / 2 / max(tot / max(n, 1), 1))


def is_complete(fasta) -> bool:
    """the whole tandem unit: GetOrganelle closed it into a circle"""
    return bool(fasta) and any(x.startswith(">") and "circular" in x for x in open(fasta))


def label_bp(fasta, label, blastn="blastn") -> int:
    """bases of the result covered by the label genes (blastn, e <= 1e-5): how much of the target unit it holds"""
    if not fasta:
        return 0
    out = subprocess.run([blastn, "-query", str(label), "-subject", str(fasta), "-evalue", "1e-5",
                          "-outfmt", "6 sseqid sstart send"], capture_output=True, text=True).stdout
    pos = set()
    for line in out.splitlines():
        sid, a, b = line.split("\t")
        a, b = sorted((int(a), int(b)))
        pos.update((sid, x) for x in range(a, b + 1))
    return len(pos)


def keep(src: Path, dst: Path):
    """result, graph and log of a GetOrganelle directory"""
    dst.mkdir(parents=True, exist_ok=True)
    for f in [result_fasta(src), *src.glob("*.selected_graph.gfa"), src / "get_org.log.txt"]:
        if f and Path(f).exists():
            shutil.copy(f, dst)


def getorganelle(kind, r1, r2, out, threads=16, max_gb=(5, 10), seed=None, label=None, given=None,
                 exe="get_organelle_from_reads.py", kmers="21,45,65,85,105", rounds=10):
    """GetOrganelle -F anonym of one multi-copy unit in rounds (module doc): round k uses max_gb[k] Gb of reads,
    round 1 the seed (or `given`, an existing GetOrganelle directory, as round 1), later rounds the previous round's
    result as seed; stops at the first complete unit that holds as much of the label genes as the best round so far
    (HNP1sub histones: round 2 "completed" a 155-bp tandem repeat while round 1 held H2A-H4 on 9.7 kb). Chosen: the
    round holding the most label bases; a complete one when it holds >= 90 % of that. The chosen round's result /
    graph / log are copied to out/, out/rounds.tsv lists the rounds; returns the result FASTA or None. A finished
    run (out/done) is reused. Bulky
    intermediates (extended reads, the extension's SPAdes run, the seed index) are removed as the user's
    remove_getorganelle_templefile*.py did."""
    out = Path(out)
    if (out / "done").exists():
        return result_fasta(out)
    seed, label = seed or SEEDS[kind][0], label or SEEDS[kind][1]
    env = dict(os.environ)                       # GetOrganelle finds bowtie2 / SPAdes / BLAST beside itself
    if Path(exe).parent != Path("."):
        env["PATH"] = f"{Path(exe).parent}{os.pathsep}{env['PATH']}"
    hist, prev = [], None
    blastn = str(Path(exe).parent / "blastn") if (Path(exe).parent / "blastn").exists() else "blastn"
    for k, gb in enumerate(max_gb, 1):
        d = out / f"round{k}"
        if k == 1 and given and result_fasta(given):
            keep(Path(given), d)
            hist.append((k, f"given ({given})", "", result_fasta(d), label_bp(result_fasta(d), label, blastn)))
        elif r1 and Path(r1).exists():
            cap = reads_per_file(r1, gb)
            cmd = [exe, "-1", str(r1), "-2", str(r2), "-o", str(d), "-F", "anonym", "-s", str(prev or seed),
                   "--genes", str(label), "-P", "0", "-t", str(threads), "-R", str(rounds), "-k", kmers,
                   "--max-reads", str(cap), "--reduce-reads-for-coverage", "inf", "--overwrite",
                   "--bowtie2-options", BOWTIE2[kind]]
            d.mkdir(parents=True, exist_ok=True)
            with open(d / "getorganelle.stdout", "w") as log:
                rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env).returncode
            # + the unzipped reads and their index, which GetOrganelle deletes itself only when it finishes
            temp = list(d.glob("[0-9]-*.fastq")) + list(d.glob("*temp.indices*"))
            for x in list(d.glob("extended_*.fq")) + temp + [d / "extended_spades", d / "seed"]:
                shutil.rmtree(x, ignore_errors=True) if x.is_dir() else x.unlink(missing_ok=True)
            f = result_fasta(d)
            hist.append((k, f"{gb} Gb ({cap} reads per file), seed {'round ' + str(k - 1) if prev else 'bundled'}",
                         rc, f, label_bp(f, label, blastn)))
        else:
            break
        top = max(h[4] for h in hist)
        if hist[-1][3] and hist[-1][4] >= 0.9 * top:
            prev = hist[-1][3]                      # the next round grows the best result, not a stray one
        if is_complete(hist[-1][3]) and hist[-1][4] >= 0.9 * top:
            break
    done = [h for h in hist if h[3] and h[4]]
    top = max((h[4] for h in done), default=0)
    best = next((h for h in done if is_complete(h[3]) and h[4] >= 0.9 * top),
                max(done, key=lambda h: h[4]) if done else None)
    if best:
        keep(best[3].parent, out)
    with open(out / "rounds.tsv", "w") as fh:
        fh.write("round\treads\texit\tresult\tlabel_bp\tcomplete\tchosen\n")
        for k, how, rc, f, bp in hist:
            chosen = best is not None and k == best[0]
            fh.write(f"{k}\t{how}\t{rc}\t{f.name if f else ''}\t{bp}\t{is_complete(f)}\t{chosen}\n")
    (out / "done").write_text(f"round {best[0]}: {best[3].name}, complete {is_complete(best[3])}\n" if best
                              else "no result\n")
    return result_fasta(out) if best else None


# ── annotation ─────────────────────────────────────────────────────────────
def mitom_version(mitom: Path) -> str:
    r = subprocess.run(["git", "-C", str(mitom), "rev-parse", "--short", "HEAD"], capture_output=True, text=True)
    return r.stdout.strip() or "unknown"


def annotate(fasta, out, kind, sample, mitom, python=sys.executable, threads=4, gbf=True):
    """`mitom nuc annotate --genes <kind>` into out/; out/<sample>.<kind>.gbf (GetOrganelle results). Skipped when
    out/mitom.version is this MitoM's commit; returns the consensus line of MitoM's log."""
    out, mitom = Path(out), Path(mitom)
    ver = mitom_version(mitom)
    stamp = out / "mitom.version"
    if not (stamp.exists() and stamp.read_text().strip() == ver and (out / "consensus.tsv").exists()):
        shutil.rmtree(out, ignore_errors=True)
        out.mkdir(parents=True)
        env = dict(os.environ, PATH=f"{Path(python).parent}{os.pathsep}{os.environ['PATH']}")
        with open(out / "mitom_nuc.log", "w") as log:
            subprocess.run([str(python), str(mitom / "mitom"), "nuc", "annotate", "-i", str(fasta), "-o", str(out),
                            "--genes", kind, "--threads", str(threads)], stdout=log, stderr=subprocess.STDOUT,
                           env=env, check=True)
        stamp.write_text(ver + "\n")
    if gbf and (out / "nuc_annotation.gb").exists():
        shutil.copy(out / "nuc_annotation.gb", out / f"{sample}.{kind}.gbf")
    line = [x for x in (out / "mitom_nuc.log").read_text().splitlines() if "consensus:" in x]
    return line[-1].split("consensus:")[1].split("(")[0].strip() if line else "no genes"


# ── comparison ─────────────────────────────────────────────────────────────
def consensus(d):
    """{gene: (sequence, partial, class)} of a MitoM output directory"""
    from Bio import SeqIO
    d, out = Path(d), {}
    if not (d / "consensus.tsv").exists():
        return out
    for r in csv.DictReader(open(d / "consensus.tsv"), delimiter="\t"):
        f = d / "consensus" / f"{r['gene']}.fasta"
        if f.exists():
            seq = str(next(SeqIO.parse(str(f), "fasta")).seq).upper()
            part = "/".join(x for x, k in (("5'", "partial_5"), ("3'", "partial_3")) if r.get(k) == "1")
            out[r["gene"]] = (seq, part, r.get("class", ""))
    return out


def overlap(a, b):
    """(aligned columns, mismatches, indel nt) of an end-gap-free alignment: one sequence may be partial"""
    from Bio import Align
    al = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-3, open_gap_score=-5,
                               extend_gap_score=-2, end_gap_score=0)
    ca, cb = al.align(a, b)[0].coordinates
    blocks = list(zip(zip(ca[:-1], ca[1:]), zip(cb[:-1], cb[1:])))
    inner = [i for i, ((a0, a1), (b0, b1)) in enumerate(blocks) if a1 > a0 and b1 > b0]
    cols = mis = indel = 0
    for (a0, a1), (b0, b1) in blocks[inner[0]:inner[-1] + 1] if inner else []:
        if a1 > a0 and b1 > b0:
            cols += a1 - a0
            mis += sum(x != y for x, y in zip(a[a0:a1], b[b0:b1]))
        else:
            indel += max(a1 - a0, b1 - b0)
    return cols, mis, indel


def compare(sample, kind, go_dir, spades_dir):
    """One row per gene: GetOrganelle vs SPAdes-MitoM consensus. verdict: identical | contained (the shorter,
    partial one is identical inside the longer) | differ | only GetOrganelle | only SPAdes-MitoM"""
    go, sp = consensus(go_dir) if go_dir else {}, consensus(spades_dir) if spades_dir else {}
    rows = []
    for g in GENES[kind]:
        if g not in go and g not in sp:
            continue
        r = {"sample": sample, "marker": kind, "gene": g}
        for k, d in (("getorganelle", go), ("spades", sp)):
            if g in d:
                r.update({f"{k}_nt": len(d[g][0]), f"{k}_partial": d[g][1], f"{k}_class": d[g][2]})
        if g in go and g in sp:
            cols, mis, indel = overlap(go[g][0], sp[g][0])
            r.update(overlap_nt=cols, mismatches=mis, indel_nt=indel,
                     identity_pct=f"{100 * (cols - mis) / max(cols + indel, 1):.2f}",
                     verdict="differ" if mis or indel else "identical" if len(go[g][0]) == len(sp[g][0])
                     else "contained")
        else:
            r["verdict"] = "only GetOrganelle" if g in go else "only SPAdes-MitoM"
        rows.append(r)
    return rows


def write_rows(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)


# ── batch over an existing results directory ───────────────────────────────
# paths relative to the results directory, {s} = sample. wgs_skim: this workflow's own layout (the SPAdes MitoM run
# of rule nuc_markers is reused); wgs_skim_pipeline: Temp_scripts/wgs_skim_pipeline (01.fastp, 02.organelle/03_nr,
# 03.spades).
LAYOUTS = {
    "wgs_skim": {"samples": "qc/*", "r1": "qc/{s}/{s}.R1.fq.gz", "r2": "qc/{s}/{s}.R2.fq.gz",
                 "rdna": "rdna/{s}/getorganelle", "scaffolds": "assembly/{s}/{s}.scaffolds.fasta",
                 "spades_mitom": "nuclear/{s}/mitom_nuc"},
    "wgs_skim_pipeline": {"samples": "01.fastp/*", "r1": "01.fastp/{s}/{s}.R1.clean.fastq.gz",
                          "r2": "01.fastp/{s}/{s}.R2.clean.fastq.gz", "rdna": "02.organelle/03_nr/{s}",
                          "scaffolds": "03.spades/{s}/scaffolds.fasta"},
}


def layout(results):
    res = Path(results)
    for name, lay in LAYOUTS.items():
        dirs = sorted(d.name for d in res.glob(lay["samples"]) if d.is_dir() and not d.name.startswith("_"))
        if dirs:
            return name, {s: {k: res / v.format(s=s) for k, v in lay.items() if k != "samples"} for s in dirs}
    sys.exit(f"[nuclear] {res}: no known results layout ({', '.join(LAYOUTS)})")


def one_sample(job):
    """Every step for one sample; returns (printed lines, comparison rows). Steps whose output exists are reused."""
    s, p, a = job["sample"], job["paths"], job["args"]
    out, lines, rows = Path(a.outdir), [], []
    try:
        for kind in a.markers:
            go_dir = out / kind / "getorganelle" / s
            # an assembly the results directory already has (rDNA) is round 1
            fasta = getorganelle(kind, p["r1"] if p["r1"].exists() else None, p["r2"], go_dir, a.threads, a.max_gb,
                                 given=p.get(kind), exe=a.getorganelle)
            ann = out / kind / "annotation" / s
            got = annotate(fasta, ann, kind, s, a.mitom, a.mitom_python) if fasta else "no GetOrganelle result"
            sp = p.get("spades_mitom") if p.get("spades_mitom") and Path(p["spades_mitom"]).exists() else None
            if not sp and p["scaffolds"].exists():
                sp = out / kind / "spades_mitom" / s
                annotate(p["scaffolds"], sp, kind, s, a.mitom, a.mitom_python, threads=12, gbf=False)
            rows += compare(s, kind, ann if fasta else None, sp)
            lines.append(f"{s} {kind}: GetOrganelle {got}" + ("" if sp else " (no SPAdes assembly)"))
    except Exception as e:      # one sample failing must not stop the batch; reported
        lines.append(f"{s}: error: {type(e).__name__}: {e}")
    return lines, rows


README = """多拷贝核基因：rDNA（18S-ITS1-5.8S-ITS2-28S）与组蛋白（H1 H2A H2B H3 H4）  wgs_skim nuclear

summary.tsv   每个样本 x 每个基因：GetOrganelle 与 SPAdes-MitoM 的长度、是否不完整（5'/3'）、组蛋白类别、
              重叠区错配 / 插入缺失、一致度、结论（identical 相同 | contained 不完整的一条完全包含于另一条 |
              differ 有差异 | only GetOrganelle | only SPAdes-MitoM）
rdna/ histone/
  getorganelle/<样本>/  GetOrganelle 组装（*.path_sequence.fasta、*.selected_graph.gfa、get_org.log.txt）
  annotation/<样本>/    MitoM 注释：<样本>.rdna.gbf / <样本>.histone.gbf（GenBank），nuc_features.tsv，
                        consensus/（每个基因一条序列），genes/（每个拷贝），submission/（NCBI table2asn）
  spades_mitom/<样本>/  MitoM 对 SPAdes 全基因组组装的注释（有 SPAdes 的样本）
组蛋白类别：replication-dependent 串联基因簇拷贝 | replacement 替换变体（H3.3、含内含子）| single-copy 单拷贝
"""


def batch(a):
    name, samples = layout(a.results)
    if a.samples:
        keep = set(a.samples.split(","))
        samples = {s: p for s, p in samples.items() if s in keep}
    print(f"[nuclear] {a.results}: layout {name}, {len(samples)} sample(s), markers {','.join(a.markers)}, "
          f"MitoM {mitom_version(a.mitom)}", flush=True)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "README.txt").write_text(README)
    rows = []
    with ProcessPoolExecutor(max_workers=max(1, a.jobs)) as ex:
        for lines, r in ex.map(one_sample, [{"sample": s, "paths": p, "args": a} for s, p in samples.items()]):
            print("\n".join(lines), flush=True)
            rows += r
    old = []                                            # keep rows of samples not run this time
    if (out / "summary.tsv").exists():
        old = [r for r in csv.DictReader(open(out / "summary.tsv"), delimiter="\t") if r["sample"] not in samples]
    write_rows(out / "summary.tsv", sorted(old + rows, key=lambda r: (r["sample"], r["marker"])))
    v = {}
    for r in rows:
        v[r["verdict"]] = v.get(r["verdict"], 0) + 1
    print(f"[nuclear] {len(rows)} gene rows: " + ", ".join(f"{k} {n}" for k, n in sorted(v.items()))
          + f" -> {out / 'summary.tsv'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    tools = argparse.ArgumentParser(add_help=False)
    tools.add_argument("--mitom", type=Path, help="MitoM checkout (with `mitom nuc`)")
    tools.add_argument("--mitom-python", default=sys.executable, help="python of MitoM's environment")
    tools.add_argument("--getorganelle", default="get_organelle_from_reads.py")
    h = sub.add_parser("assemble", parents=[tools])
    h.add_argument("--kind", choices=GENES, required=True)
    h.add_argument("-1", dest="r1", required=True)
    h.add_argument("-2", dest="r2", required=True)
    h.add_argument("-o", "--outdir", required=True)
    h.add_argument("-t", "--threads", type=int, default=16)
    h.add_argument("--max-gb", default="5,10", help="Gb of reads per round (default 5,10: a second round only when the "
                   "first gives no complete unit, seeded with it)")
    h.add_argument("--kmers", default="21,45,65,85,105")
    h.add_argument("--rounds", type=int, default=10)
    h.add_argument("--seed", help="seed FASTA (default resources/<kind>)")
    h.add_argument("--label", help="label FASTA (default resources/<kind>)")
    an = sub.add_parser("annotate", parents=[tools])
    an.add_argument("-i", "--fasta", required=True)
    an.add_argument("-o", "--outdir", required=True)
    an.add_argument("--kind", choices=GENES, required=True)
    an.add_argument("--sample", required=True)
    an.add_argument("-t", "--threads", type=int, default=4)
    c = sub.add_parser("compare")
    c.add_argument("--sample", required=True)
    c.add_argument("--getorganelle-rdna")
    c.add_argument("--getorganelle-histone")
    c.add_argument("--spades", help="MitoM output on the SPAdes scaffolds (both markers or one)")
    c.add_argument("-o", "--out", required=True)
    b = sub.add_parser("batch", parents=[tools])
    b.add_argument("--results", required=True, help="results directory of a known layout")
    b.add_argument("-o", "--outdir", required=True)
    b.add_argument("--samples", help="only these samples (comma-separated)")
    b.add_argument("--markers", default="rdna,histone")
    b.add_argument("-j", "--jobs", type=int, default=3, help="samples at a time (GetOrganelle 16 threads each)")
    b.add_argument("-t", "--threads", type=int, default=16)
    b.add_argument("--max-gb", default="5,10", help="Gb of reads per GetOrganelle round")
    a = ap.parse_args()
    if hasattr(a, "max_gb"):
        a.max_gb = [float(x) for x in a.max_gb.split(",") if x.strip()]
    if a.cmd == "assemble":
        f = getorganelle(a.kind, a.r1, a.r2, a.outdir, a.threads, a.max_gb, a.seed, a.label, exe=a.getorganelle,
                         kmers=a.kmers, rounds=a.rounds)
        print(f or "no result")
    elif a.cmd == "annotate":
        print(annotate(a.fasta, a.outdir, a.kind, a.sample, a.mitom, a.mitom_python, a.threads))
    elif a.cmd == "compare":
        rows = []
        for kind, d in (("rdna", a.getorganelle_rdna), ("histone", a.getorganelle_histone)):
            rows += compare(a.sample, kind, d, a.spades)
        write_rows(a.out, rows)
    else:
        a.markers = [m for m in a.markers.split(",") if m in GENES]
        batch(a)


if __name__ == "__main__":
    main()

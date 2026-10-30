"""Multi-copy nuclear markers - rDNA (18S-ITS1-5.8S-ITS2-28S) and histones (H1, H2A, H2B, H3, H4): GetOrganelle
assembly of the tandem unit, MitoM annotation (GenBank .gbf), GetOrganelle vs SPAdes-MitoM comparison.

Why both: GetOrganelle walks the multi-copy unit in its graph and returns it whole (gene order, spacers); SPAdes
scatters the copies over many contigs, which MitoM merges into one consensus per gene. The same genome must give
the same genes, so differences point at assembly or annotation problems (MJ20260820041: the SPAdes histone cluster
is split into 84-134 nt contigs; GetOrganelle from H1-H4 seeds gave the cluster H3 whole).

Histones (MJ20260820041, CHY006): a label of H3 only kept just the H3 piece of the cluster (2.2 kb); H1-H4 cluster
copies as seed / label (resources/histone) gave all five genes.

Reads (user, 2026-10-28): multi-copy units need little data - round 1 uses at most 5 Gb; when the unit is not
complete, round 2 uses 10 Gb with round 1's result as the seed (the iteration of mt_iterate.py). Complete = the
path is circular, or every core gene of the unit (rDNA: 18S, 5.8S, 28S; histones: H2A, H2B, H3, H4 - H1 is often
outside the cluster) is covered end to end by a label sequence (blastn, >= 95 % of the label): a linear rDNA unit
with all genes is complete (MJ20260820041: 18 of 22 rDNA units are linear with every gene whole; a 10 Gb round
on them would cost 15-40 min each for nothing). An assembly the results directory already has (rDNA of
wgs_skim_pipeline) counts as round 1.

Subcommands (also used by the workflow rules, workflow/rules/nuclear.smk):
  assemble  GetOrganelle -F anonym of a multi-copy unit (--kind rdna | histone) from clean reads, in rounds
  annotate  `mitom nuc annotate` of a GetOrganelle result or SPAdes scaffolds -> <sample>.<kind>.gbf + MitoM tables
  compare   GetOrganelle vs SPAdes-MitoM per gene -> TSV
  merge     per-sample comparison tables -> summary.tsv
  status    what a batch has done so far
Batches: workflow/nuclear_batch.smk (every sample of an existing results directory, LAYOUTS; `wgs_skim nuclear`)
and workflow/rules/nuclear.smk (within the workflow) - Snakemake schedules, logs and resumes.
"""
import argparse
import csv
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from collect_fasta import clean, getorganelle_cmd, result_fasta  # noqa: E402,F401  (result_fasta: also from here)
from runtools import blastn, output, run, tool_env  # noqa: E402
from seqtools import is_circular as circular_header, iter_fasta  # noqa: E402
from tables import write_tsv  # noqa: E402

RES = Path(__file__).resolve().parents[2] / "resources"


def log(msg):
    import time
    print(f"[nuclear {time.strftime('%H:%M:%S')}] {msg}", flush=True)



GENES = {"rdna": ["18S", "ITS1", "5.8S", "ITS2", "28S"], "histone": ["H1", "H2A", "H2B", "H3", "H4"]}
# bowtie2 options for GetOrganelle's seed reads: histone DNA of another family is too far from the bundled seed for
# end-to-end mapping (HNP1sub, a Pisione: 391 bytes of seed reads, "Too few seed reads", no result; local
# --very-sensitive-local: 153 kb of seed reads, H2A H2B H3 H4 complete on 9.7 kb); rDNA: GetOrganelle's default
BOWTIE2 = {"rdna": "--very-fast -t", "histone": "--very-sensitive-local -t"}
SEEDS = {"rdna": (RES / "rdna" / "seed.fasta", RES / "rdna" / "label.fasta"),
         "histone": (RES / "histone" / "seed.fasta", RES / "histone" / "label.fasta")}
CORE = {"rdna": ["18S", "5.8S", "28S"], "histone": ["H2A", "H2B", "H3", "H4"]}
COLS = ["sample", "marker", "gene", "getorganelle_nt", "getorganelle_partial", "getorganelle_class", "spades_nt",
        "spades_partial", "spades_class", "overlap_nt", "mismatches", "indel_nt", "identity_pct", "verdict"]


# ── assembly ───────────────────────────────────────────────────────────────
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


def is_circular(fasta) -> bool:
    """any record of the FASTA is circular"""
    return bool(fasta) and any(circular_header(n) for n, _ in iter_fasta(fasta, upper=False))


def label_gene(name: str, kind: str):
    """gene of a label / seed sequence from its header: rDNA '18SribosomalRNA ...', histones 'H3|sample|...' or an
    NCBI title 'histone H3 (H3) gene'; None for whole-cluster records and unnamed genes"""
    first = name.split()[0]
    if kind == "rdna":
        return next((g for g in ("5.8S", "18S", "28S") if first.startswith(g)), None)
    m = re.match(r"(H1|H2A|H2B|H3|H4)\|", first) or re.search(r"histone (H1|H2A|H2B|H3|H4)\b", name)
    return m.group(1) if m and "cluster" not in name and name.count("histone") <= 1 else None


def genes_covered(fasta, label, kind, blastn_exe="blastn"):
    """{gene: best label coverage %} of the result: how much of each gene it holds (blastn qcovs of the labels)"""
    if not fasta or not Path(fasta).stat().st_size:
        return {}
    names = {n.split()[0]: n for n, _ in iter_fasta(label, upper=False)}
    cov = {}
    for h in blastn(label, fasta, "qseqid qcovs", 1e-5, exe=blastn_exe):
        g = label_gene(names.get(h["qseqid"], h["qseqid"]), kind)
        if g:
            cov[g] = max(cov.get(g, 0), int(h["qcovs"]))
    return cov


def is_complete(fasta, kind=None, label=None, blastn_exe="blastn", cov=None) -> bool:
    """the whole tandem unit: circular, or every core gene covered end to end by a label (module doc)"""
    if is_circular(fasta):
        return True
    if not fasta or not kind:
        return False
    cov = cov if cov is not None else genes_covered(fasta, label or SEEDS[kind][1], kind, blastn_exe)
    return all(cov.get(g, 0) >= 95 for g in CORE[kind])


def label_bp(fasta, label, blastn_exe="blastn") -> int:
    """bases of the result covered by the label genes (blastn, e <= 1e-5): how much of the target unit it holds"""
    if not fasta or not Path(fasta).stat().st_size:
        return 0
    pos = set()
    for h in blastn(label, fasta, "sseqid sstart send", 1e-5, exe=blastn_exe):
        a, b = sorted((int(h["sstart"]), int(h["send"])))
        pos.update((h["sseqid"], x) for x in range(a, b + 1))
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
    env = tool_env(exe)                          # GetOrganelle finds bowtie2 / SPAdes / BLAST beside itself
    hist, prev = [], None
    bn = str(Path(exe).parent / "blastn") if (Path(exe).parent / "blastn").exists() else "blastn"
    for k, gb in enumerate(max_gb, 1):
        d = out / f"round{k}"
        if k == 1 and given and result_fasta(given):
            keep(Path(given), d)
            f = result_fasta(d)
            hist.append((k, f"given ({given})", "", f, label_bp(f, label, bn), genes_covered(f, label, kind, bn)))
        elif r1 and Path(r1).exists():
            cap = reads_per_file(r1, gb)
            cmd = getorganelle_cmd(r1, r2, d, "anonym", threads, prev or seed, rounds, kmers, cap, genes=label,
                                   extra=["-P", "0", "--bowtie2-options", BOWTIE2[kind]], exe=exe)
            d.mkdir(parents=True, exist_ok=True)
            seed_how = f"round {k - 1}" if prev else "bundled"
            log(f"{out.parent.name}/{out.name}: {kind} round {k}, {gb} Gb, seed {seed_how}")
            rc = run(cmd, d / "getorganelle.stdout", env=env)
            clean(d)
            f = result_fasta(d)
            hist.append((k, f"{gb} Gb ({cap} reads per file), seed {'round ' + str(k - 1) if prev else 'bundled'}",
                         rc, f, label_bp(f, label, bn), genes_covered(f, label, kind, bn)))
        else:
            break
        top = max(h[4] for h in hist)
        if hist[-1][3] and hist[-1][4] > 0 and hist[-1][4] >= 0.9 * top:
            prev = hist[-1][3]                      # the next round grows the best result, not a stray one
        if is_complete(hist[-1][3], kind, cov=hist[-1][5]) and hist[-1][4] > 0 and hist[-1][4] >= 0.9 * top:
            break
    done = [h for h in hist if h[3] and h[4]]
    top = max((h[4] for h in done), default=0)
    best = next((h for h in done if is_complete(h[3], kind, cov=h[5]) and h[4] >= 0.9 * top),
                max(done, key=lambda h: h[4]) if done else None)
    if not hist or (not best and all(h[2] not in ("", 0) for h in hist)):
        raise RuntimeError(f"{kind}: no reads and no given assembly" if not hist else
                           f"{kind}: GetOrganelle failed in every round (exit {[h[2] for h in hist]}); see "
                           f"{out}/round*/getorganelle.stdout")          # no `done` stamp: the next run tries again
    out.mkdir(parents=True, exist_ok=True)
    if best:
        keep(best[3].parent, out)
    with open(out / "rounds.tsv", "w") as fh:
        fh.write("round\treads\texit\tresult\tlabel_bp\tgenes_covered_pct\tcircular\tcomplete\tchosen\n")
        for k, how, rc, f, bp, cov in hist:
            chosen = best is not None and k == best[0]
            genes = ",".join(f"{g}:{cov.get(g, 0)}" for g in CORE[kind] + (["H1"] if kind == "histone" else []))
            fh.write(f"{k}\t{how}\t{rc}\t{f.name if f else ''}\t{bp}\t{genes}\t{is_circular(f)}\t"
                     f"{is_complete(f, kind, cov=cov)}\t{chosen}\n")
    (out / "done").write_text(f"round {best[0]}: {best[3].name}, complete {is_complete(best[3], kind, cov=best[5])}, "
                              f"circular {is_circular(best[3])}\n" if best else "no result\n")
    return result_fasta(out) if best else None


# ── annotation ─────────────────────────────────────────────────────────────
def mitom_version(mitom: Path) -> str:
    """git commit of the MitoM checkout; for a plain copy the size and time of mito_nuc.py"""
    commit = output(["git", "-C", mitom, "rev-parse", "--short", "HEAD"], check=False).strip()
    if commit:
        return commit
    f = Path(mitom) / "mito_nuc.py"
    return f"mito_nuc.py:{f.stat().st_size}:{int(f.stat().st_mtime)}" if f.exists() else "unknown"


def annotate(fasta, out, kind, sample, mitom, python=sys.executable, threads=4, gbf=True):
    """`mitom nuc annotate --genes <kind>` into out/; out/<sample>.<kind>.gbf (GetOrganelle results). Skipped when
    out/mitom.version is this MitoM's commit; returns the consensus line of MitoM's log."""
    out, mitom = Path(out), Path(mitom)
    st = Path(fasta).stat()
    ver = f"{mitom_version(mitom)} input:{st.st_size}:{int(st.st_mtime)}"   # redone when MitoM or the input changed
    stamp = out / "mitom.version"
    old = stamp.read_text().strip() if stamp.exists() else ""
    if not (old in (ver, ver.split(" ")[0]) and (out / "consensus.tsv").exists()):
        shutil.rmtree(out, ignore_errors=True)
        out.mkdir(parents=True)
        log(f"{sample}: MitoM {kind} on {Path(fasta).name} -> {out.parent.name}/{out.name}")
        run([python, mitom / "mitom", "nuc", "annotate", "-i", fasta, "-o", out, "--genes", kind, "--threads",
             threads], out / "mitom_nuc.log", env=tool_env(python), check=True)
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
            out.setdefault(r["gene"], (seq, part, r.get("class", "")))   # first organism = first record
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
    # a gene MitoM lists with an empty consensus was not assembled (MJ20260820041 CHY002 / CHY008 / CHY010: SPAdes
    # ITS1 / ITS2 of 0 nt - the spacer between contigs; aligning it raised "sequence has zero length" and lost the
    # sample's whole rDNA comparison)
    empty = [f"{k} {g}" for k, d in (("GetOrganelle", go), ("SPAdes-MitoM", sp)) for g in list(d) if not d[g][0]]
    go, sp = {g: v for g, v in go.items() if v[0]}, {g: v for g, v in sp.items() if v[0]}
    rows = []
    for g in GENES[kind]:
        if g not in go and g not in sp:
            continue
        r = {"sample": sample, "marker": kind, "gene": g}
        for k, d in (("getorganelle", go), ("spades", sp)):
            if g in d:
                r.update({f"{k}_nt": len(d[g][0]), f"{k}_partial": d[g][1], f"{k}_class": d[g][2]})
        if g in go and g in sp:
            a, b = go[g][0], sp[g][0]
            cols, mis, indel = overlap(a, b)
            short = min(len(a), len(b))
            # free end gaps: unrelated sequences align over a few bases with 0 mismatches - the shorter sequence
            # must be covered over >= 90 % of its length before "identical" / "contained" means anything
            r.update(overlap_nt=cols, mismatches=mis, indel_nt=indel,
                     identity_pct=f"{100 * (cols - mis) / max(short, 1):.2f}",
                     verdict="no overlap" if cols < 0.5 * short else "differ" if mis or indel or cols < 0.9 * short
                     else "identical" if len(a) == len(b) else "contained")
        else:
            r["verdict"] = "only GetOrganelle" if g in go else "only SPAdes-MitoM"
            if any(e.endswith(f" {g}") for e in empty):
                r["verdict"] += " (" + ", ".join(e for e in empty if e.endswith(f" {g}")) + ": empty consensus)"
        rows.append(r)
    return rows


def write_rows(path, rows):
    write_tsv(path, rows, COLS)


# ── batch over an existing results directory ───────────────────────────────
# paths relative to the results directory, {s} = sample. wgs_skim: this workflow's own layout (the SPAdes MitoM run
# of rule nuc_markers is reused); wgs_skim_pipeline: Temp_scripts/wgs_skim_pipeline (01.fastp, 02.organelle/03_nr,
# 03.spades).
LAYOUTS = {
    "wgs_skim": {"samples": "rdna/*", "r1": "fastp/{s}/{s}.R1.fq.gz", "r2": "fastp/{s}/{s}.R2.fq.gz",
                 "rdna": "rdna/{s}/getorganelle", "histone": "histone/{s}/getorganelle",
                 "scaffolds": "spades/{s}/{s}.scaffolds.fasta",
                 "spades_mitom": "nuclear_markers/{s}/mitom_nuc"},
    "wgs_skim_0.1": {"samples": "qc/*", "r1": "qc/{s}/{s}.R1.fq.gz", "r2": "qc/{s}/{s}.R2.fq.gz",   # wgs_skim <= 0.1
                     "rdna": "rdna/{s}/getorganelle", "histone": "histone/{s}/getorganelle",
                     "scaffolds": "assembly/{s}/{s}.scaffolds.fasta",
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


README = """多拷贝核基因：rDNA（18S-ITS1-5.8S-ITS2-28S）与组蛋白（H1 H2A H2B H3 H4）  wgs_skim nuclear

summary.tsv   每个样本 x 每个基因：GetOrganelle 与 SPAdes-MitoM 的长度、是否不完整（5'/3'）、组蛋白类别、
              重叠区错配 / 插入缺失、一致度、结论（identical 相同 | contained 不完整的一条完全包含于另一条 |
              differ 有差异 | only GetOrganelle | only SPAdes-MitoM）
rdna/ histone/
  getorganelle/<样本>/  GetOrganelle 组装（*.path_sequence.fasta、*.selected_graph.gfa、get_org.log.txt）
  annotation/<样本>/    MitoM 注释：<样本>.rdna.gbf / <样本>.histone.gbf（GenBank），nuc_features.tsv，
                        consensus/（每个基因一条序列），genes/（每个拷贝），submission/（NCBI table2asn）
spades_mitom/<样本>/    MitoM 对 SPAdes 全基因组组装的注释，rDNA 与组蛋白一次完成（有 SPAdes 的样本）
组蛋白类别：replication-dependent 串联基因簇拷贝 | replacement 替换变体（H3.3、含内含子）| single-copy 单拷贝
"""


def status(out, markers=("rdna", "histone"), samples=None):
    """What is done per sample, from the files (works while a batch runs): one line per sample"""
    out = Path(out)
    names = samples or sorted({d.name for k in markers for sub in ("getorganelle", "annotation")
                               for d in (out / k / sub).glob("*") if d.is_dir()})
    done = {}
    for t in [out / "summary.tsv", *sorted((out / "compare").glob("*.tsv"))]:
        if t.exists() and t.stat().st_size:
            for r in csv.DictReader(open(t), delimiter="\t"):
                done.setdefault((r["sample"], r["marker"]), []).append(r["verdict"])
    lines = []
    for s in names:
        parts = []
        for k in markers:
            g, ann = out / k / "getorganelle" / s / "done", out / k / "annotation" / s / "mitom.version"
            sp = next((d for d in (out / k / "spades_mitom" / s, out / "spades_mitom" / s)
                       if (d / "mitom.version").exists()), out / "none")
            sp = sp / "mitom.version"
            v = done.get((s, k), [])
            step = ("error" if any(x.startswith("error") for x in v) else "compared" if v else
                    "annotated" if ann.exists() else "assembled" if g.exists()
                    else "started" if (out / k / "getorganelle" / s).exists() else "waiting")
            parts.append(f"{k} {step}" + (" +spades" if sp.exists() else "")
                         + (f" [{', '.join(sorted({x[:40] for x in v}))}]" if v else ""))
        lines.append(f"{s:16s} " + " | ".join(parts))
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    tools = argparse.ArgumentParser(add_help=False)
    tools.add_argument("--mitom", type=Path, help="MitoM checkout (with `mitom nuc`)")
    tools.add_argument("--mitom-python", default=sys.executable, help="python of MitoM's environment")
    tools.add_argument("--getorganelle", default="get_organelle_from_reads.py")
    h = sub.add_parser("assemble", parents=[tools])
    h.add_argument("--kind", choices=GENES, required=True)
    h.add_argument("-1", dest="r1", help="clean reads (optional with --given)")
    h.add_argument("-2", dest="r2")
    h.add_argument("--given", help="an existing GetOrganelle directory of this unit: round 1")
    h.add_argument("-o", "--outdir", required=True)
    h.add_argument("-t", "--threads", type=int, default=16)
    h.add_argument("--max-gb", default="5,10", help="Gb of reads per round (default 5,10: a second round only when the "
                   "first gives no complete unit, seeded with it)")
    h.add_argument("--kmers", default="21,45,65,85,105")
    h.add_argument("--rounds", type=int, default=10)
    h.add_argument("--seed", help="seed FASTA (default resources/<kind>)")
    h.add_argument("--label", help="label FASTA (default resources/<kind>)")
    an = sub.add_parser("annotate", parents=[tools])
    an.add_argument("-i", "--fasta", required=True,
                    help="FASTA, or a GetOrganelle directory (its result; none: a 'no result' stamp)")
    an.add_argument("-o", "--outdir", required=True)
    an.add_argument("--kind", required=True, help="rdna, histone or rdna,histone (one MitoM run for both)")
    an.add_argument("--sample", required=True)
    an.add_argument("-t", "--threads", type=int, default=4)
    an.add_argument("--no-gbf", action="store_true", help="no <sample>.<kind>.gbf (SPAdes scaffolds)")
    c = sub.add_parser("compare")
    c.add_argument("--sample", required=True)
    c.add_argument("--getorganelle-rdna")
    c.add_argument("--getorganelle-histone")
    c.add_argument("--spades", help="MitoM output on the SPAdes scaffolds (both markers or one)")
    c.add_argument("--markers", default="rdna,histone")
    c.add_argument("-o", "--out", required=True)
    m = sub.add_parser("merge", help="per-sample comparison tables -> one summary.tsv (+ README.txt)")
    m.add_argument("tables", nargs="*")
    m.add_argument("-o", "--out", required=True)
    st = sub.add_parser("status", help="what a batch has done so far (from its output directory)")
    st.add_argument("-o", "--outdir", required=True)
    st.add_argument("--markers", default="rdna,histone")
    a = ap.parse_args()
    if a.cmd == "status":
        print("\n".join(status(a.outdir, a.markers.split(","))))
        return
    if hasattr(a, "max_gb"):
        a.max_gb = [float(x) for x in a.max_gb.split(",") if x.strip()]
    if a.cmd == "assemble":
        f = getorganelle(a.kind, a.r1 if a.r1 and Path(a.r1).exists() else None, a.r2, a.outdir, a.threads,
                         a.max_gb, a.seed, a.label, given=a.given, exe=a.getorganelle, kmers=a.kmers, rounds=a.rounds)
        print(f or "no result")
    elif a.cmd == "annotate":
        fasta = result_fasta(a.fasta) if Path(a.fasta).is_dir() else Path(a.fasta)
        if not fasta or not fasta.stat().st_size:                # nothing assembled: a stamp, no MitoM run
            Path(a.outdir).mkdir(parents=True, exist_ok=True)
            (Path(a.outdir) / "mitom.version").write_text("no result\n")
            print("no result")
            return
        print(annotate(fasta, a.outdir, a.kind, a.sample, a.mitom, a.mitom_python, a.threads, gbf=not a.no_gbf))
    elif a.cmd == "compare":
        rows = []
        for kind, d in (("rdna", a.getorganelle_rdna), ("histone", a.getorganelle_histone)):
            if kind in a.markers.split(","):
                rows += compare(a.sample, kind, d, a.spades)
        if not rows:                                             # an empty table still says what was compared
            rows = [{"sample": a.sample, "marker": m, "gene": "", "verdict": "no genes"} for m in a.markers.split(",")]
        write_rows(a.out, rows)
    elif a.cmd == "merge":
        rows = [r for t in sorted(a.tables) if Path(t).stat().st_size
                for r in csv.DictReader(open(t), delimiter="\t")]
        write_rows(a.out, rows)
        (Path(a.out).parent / "README.txt").write_text(README)
        print(f"{len(rows)} gene rows from {len(a.tables)} sample(s) -> {a.out}")


if __name__ == "__main__":
    main()

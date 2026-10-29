"""Compare the results of tests/all_features.sh with the truth of tests/sim/simulate.py.

    check_truth.py SIM_DIR OUT_DIR        -> OUT_DIR/truth_check.tsv and truth_check.md

For every function directory (OUT_DIR/NN_name/) and sample, whatever it produced:
  mitogenome   mito/<s>/<s>.mt.fasta: length, topology, % of the true mitogenome covered and identity (blastn)
  rDNA/histone rdna|histone/<s>/*.fasta: % of the true unit covered and identity
  assembly     assembly/<s>/<s>.scaffolds.fasta: total Mb, N50 (truth: nuclear_mb)
  genome size  genome_size/<s>/<s>.genomescope.txt: GenomeScope haploid length (truth: nuclear_mb + organelles)
  MitoM        annotation/<s>/<s>.mt.gb: CDS / rRNA / tRNA genes (truth: 13 / 2 / 22)
  MitoFinder   mitofinder/<src>/mitofinder_summary.tsv: genes found (truth: 15 = 13 CDS + 2 rRNA)
  contamination mitofinder/<src>/blast/samples_check.tsv: flagged genes (truth: the rna_contaminant column)
verdict: ok | check (differs from the truth; the note says how).
"""
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow" / "scripts"))
from seqtools import fasta_lengths, is_circular, iter_fasta, n50  # noqa: E402
from tables import read_tsv, write_tsv  # noqa: E402

COLS = ["function", "sample", "item", "result", "truth", "verdict", "note"]


def covered(truth_fa, fasta):
    """(% of the truth covered by the fasta, identity %) - blastn, HSPs merged on the truth"""
    out = subprocess.run(["blastn", "-query", str(truth_fa), "-subject", str(fasta), "-evalue", "1e-20",
                          "-outfmt", "6 qstart qend pident length qlen"], capture_output=True, text=True).stdout
    pos, ident, aln, qlen = set(), 0.0, 0, 1
    for line in out.splitlines():
        a, b, pid, ln, qlen = line.split("\t")
        a, b = sorted((int(a), int(b)))
        pos.update(range(a, b + 1))
        ident += float(pid) * int(ln)
        aln += int(ln)
        qlen = int(qlen)
    return 100 * len(pos) / int(qlen), (ident / aln if aln else 0.0)


def row(fn, s, item, result, truth, ok, note=""):
    return dict(function=fn, sample=s, item=item, result=result, truth=truth, verdict="ok" if ok else "check",
                note=note)


def main():
    sim, out = Path(sys.argv[1]), Path(sys.argv[2])
    truth = {r["sample"]: r for r in read_tsv(sim / "truth.tsv")}
    rows = []
    for d in sorted(p for p in out.iterdir() if p.is_dir() and re.match(r"\d\d_", p.name)):
        fn = d.name
        for s, t in truth.items():
            tdir = sim / "truth" / s
            f = d / "mito" / s / f"{s}.mt.fasta"
            if f.exists() and f.stat().st_size:
                recs = list(iter_fasta(f))
                L = sum(len(x) for _, x in recs)
                circ = all(is_circular(n) for n, _ in recs) and len(recs) == 1
                cov, pid = covered(tdir / "mt.fasta", f)
                topo = "circular" if circ else "linear"
                ok = cov >= 99 and pid >= 99.5 and topo == t["mt_topology"]
                rows.append(row(fn, s, "mitogenome", f"{L} bp {topo}, {len(recs)} seq, {cov:.1f}% covered, "
                                f"{pid:.2f}% identity", f"{t['mt_length']} bp {t['mt_topology']}", ok))
            for kind in ("rdna", "histone"):
                fas = [x for x in sorted((d / kind / s).glob("*.fasta")) if x.stat().st_size] \
                    if (d / kind / s).is_dir() else []
                if fas:
                    cov, pid = covered(tdir / f"{kind}.fasta", fas[0])
                    rows.append(row(fn, s, kind, f"{fas[0].name}: {cov:.1f}% of the unit, {pid:.2f}% identity",
                                    f"{t[f'{kind}_unit']} bp unit", cov >= 95 and pid >= 99))
            f = d / "assembly" / s / f"{s}.scaffolds.fasta"
            if f.exists() and f.stat().st_size:
                lens = fasta_lengths(f)
                tot = sum(lens) / 1e6
                rows.append(row(fn, s, "assembly", f"{tot:.1f} Mb, N50 {n50(lens)}, {len(lens)} scaffolds",
                                f"{t['nuclear_mb']} Mb nuclear", 0.8 * float(t["nuclear_mb"]) <= tot <=
                                1.3 * float(t["nuclear_mb"]) + 0.1))
            f = d / "genome_size" / s / f"{s}.genomescope.txt"
            if f.exists():
                m = re.search(r"Genome Haploid Length\s+([\d,]+) bp\s+([\d,]+) bp", f.read_text())
                mb = int(m.group(2).replace(",", "")) / 1e6 if m else 0
                rows.append(row(fn, s, "genome size", f"{mb:.1f} Mb" if m else "no fit", f"{t['nuclear_mb']} Mb",
                                m and 0.7 * float(t["nuclear_mb"]) <= mb <= 1.5 * float(t["nuclear_mb"])))
            f = d / "annotation" / s / f"{s}.mt.gb"
            if f.exists() and f.stat().st_size:
                txt = f.read_text()
                n = {k: len(re.findall(rf"^ {{5}}{k} ", txt, re.M)) for k in ("CDS", "rRNA", "tRNA")}
                rows.append(row(fn, s, "MitoM genes", f"{n['CDS']} CDS / {n['rRNA']} rRNA / {n['tRNA']} tRNA",
                                "13 / 2 / 22", n["CDS"] == 13 and n["rRNA"] == 2 and n["tRNA"] >= 20))
        for summ in sorted((d / "mitofinder").glob("*/mitofinder_summary.tsv")) if (d / "mitofinder").is_dir() else []:
            src = summ.parent.name
            for r in read_tsv(summ):
                if r["sample"] not in truth:
                    continue
                genes = [k for k, v in r.items() if k not in ("sample", "contigs", "length", "circular",
                                                              "protein_genes", "rRNAs", "tRNAs") and v != "N"]
                found = int(r["protein_genes"] or 0) + int(r["rRNAs"] or 0)
                rows.append(row(fn, r["sample"], f"MitoFinder {src}", f"{found} genes ({r['contigs']} contigs, "
                                f"{r['tRNAs']} tRNAs)", "15 genes", found >= 15, ",".join(genes)))
            check = summ.parent / "blast" / "samples_check.tsv"
            for r in read_tsv(check):
                t = truth.get(r["sample"])
                if not t:
                    continue
                # the contaminant is in the RNA reads only (simulate.py --contaminant)
                want = {g.replace("contaminant_", "").upper() for g in t["rna_contaminant"].split(",") if g} \
                    if src in ("rna", "trinity") else set()
                v = r["verdict"]
                flagged = {g.strip().upper() for g in v.split(":", 1)[1].split(",")} \
                    if v.startswith("possible contamination") else set()
                ok = flagged == want
                rows.append(row(fn, r["sample"], f"contamination {src}", r["verdict"],
                                f"contaminant: {','.join(sorted(want)) or 'none'}", ok, r.get("flagged", "")))
    write_tsv(out / "truth_check.tsv", rows, COLS)
    bad = [r for r in rows if r["verdict"] != "ok"]
    with open(out / "truth_check.md", "w") as fh:
        fh.write(f"# 与真值比较：{len(rows)} 项，{len(rows) - len(bad)} ok，{len(bad)} 需要检查\n\n")
        fh.write("| 功能 | 样品 | 项目 | 结果 | 真值 | 判断 | 备注 |\n|---|---|---|---|---|---|---|\n")
        for r in rows:
            fh.write("| " + " | ".join(str(r[c]).replace("|", "/") for c in COLS) + " |\n")
    print(f"{len(rows)} checks, {len(bad)} to check -> {out / 'truth_check.md'}")


if __name__ == "__main__":
    main()

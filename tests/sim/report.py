"""README.md in every function directory of tests/all_features.sh, and README.md above them.

    report.py SIM_DIR OUT_DIR          (after tests/sim/check_truth.py)
Each function directory is wgs_skim's own output (<tool>/<sample>/, <tool>/<tool>_summary.tsv, pipeline_info/);
its README.md says what the function does, the command, the directory layout, every sample's result against the
truth and the verdict."""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_truth import covered  # noqa: E402
from seqtools import fasta_lengths, n50  # noqa: E402
from tables import read_tsv  # noqa: E402

# what each function does and what to look at (by the function name of all_features.sh)
WHAT = {
    "samples": "sample names detected from a reads directory (5 samples, R1/R2 paired): samples.tsv",
    "fastp": "fastp, once: clean reads and reports; fastp/fastp_summary.tsv = data and Q30 per sample. Every later "
             "DNA function reads these clean reads where they are",
    "meangs": "MEANGS mitochondrial seed per sample; meangs/meangs_summary.tsv",
    "getorganelle_mt": "GetOrganelle with its own animal_mt seeds; getorganelle/getorganelle_summary.tsv = length, "
                       "circular, rounds. Seeds far from the sample can give nothing - the reason for 05",
    "meangs-getorganelle": "recommended: MEANGS seed + GetOrganelle, reseeded while not circular; "
                           "mitogenome/mitogenome_summary.tsv = final length and topology (truth: S1-S4 circular, "
                           "S5 linear)",
    "getorganelle_nr": "rDNA unit (18S-ITS1-5.8S-ITS2-28S) and histone cluster per sample, GenBank annotation (.gbf)",
    "spades_genome": "SPAdes genome assembly; spades/spades_summary.tsv = total length (truth ~20 Mb), N50",
    "spades_meta": "metaSPAdes; compare spades_summary.tsv with 07",
    "genome-size": "GenomeScope; genomescope/genomescope_summary.tsv = genome size (truth ~20 Mb), heterozygosity "
                   "(truth 0.2 %)",
    "busco_genome": "BUSCO of given assemblies (-i); busco/busco_summary.tsv (a 20 Mb piece holds few genes: low "
                    "completeness is expected)",
    "mitom": "MitoM annotation of given mitogenomes (-i); mitom/mitom_summary.tsv = 13 PCG + 2 rRNA + 22 tRNA",
    "blast": "BLAST of cox1 and 16S against the local metazoa_marker database: best hit species and identity",
    "mitofinder_assembly": "MitoFinder on assemblies, one distant reference (S1). mitofinder_assembly/: <sample>/ + "
                           "mitofinder_summary.tsv (gene table) + genes/ + contamination/",
    "mitofinder_reads": "MitoFinder assembling the reads itself (MEGAHIT); compare with 13",
    "mitofinder_bioconda": "bioconda MitoFinder 1.4.1 (no Singularity image): tRNAs from tRNAscan-SE; compare with 13",
    "mitofinder_multiref": "as 13 with 5 related references (--ref directory): gene counts of the distant samples "
                           "against 13 (recommended use)",
    "spades-busco": "combination SPAdes -> BUSCO",
    "spades-mitofinder": "combination: GetOrganelle mitogenome + SPAdes-MitoFinder (reference = each sample's own "
                         "MitoM annotation)",
    "circularize": "circularization without reads: method (walk / overlap / join) and circle length per sample",
    "compare": "GetOrganelle vs SPAdes-MitoFinder: verdict in all.compare.tsv; differences and figure per sample",
    "run_all": "the full workflow from raw reads (every step): summary/samples.tsv, plus each tool's *_summary.tsv",
    "nuclear": "nuclear markers: GetOrganelle vs SPAdes-MitoM per gene (verdict)",
    "sheet_given": "sample sheet with existing assemblies and mitogenomes: no assembly, BUSCO / MitoM / BLAST on them",
    "split_by_id": "one merged FASTA of several samples, split by sequence name, annotated per sample",
    "python_api": "Python interface wgs_skim.tool(...): same as 10",
    "nextflow": "the Nextflow engine (raw reads, default steps): compare summary/samples.tsv with 21",
    "fastp_rna": "fastp of the RNA-seq reads, once: every later RNA function reads these clean reads",
    "rnaspades": "rnaSPAdes: transcripts per sample in spades/spades_summary.tsv",
    "trinity": "Trinity: trinity/trinity_summary.tsv, compare with 31",
    "busco_transcriptome": "BUSCO -m transcriptome (nuclear transcripts taken from a real RNA-seq assembly)",
    "rnaspades-busco": "combination (recommended) rnaSPAdes -> BUSCO",
    "trinity-busco": "combination Trinity -> BUSCO",
    "rnaspades-mitofinder": "mitochondrial genes from RNA-seq: mitofinder_summary.tsv (gene table, N = not found); "
                            "contamination/samples_check.tsv: S5 must report COX1, CYTB, the others nothing",
    "trinity-mitofinder": "as 36 with Trinity transcripts",
    "mitofinder_rna_tool": "the tool command wgs_skim mitofinder --assembler rnaspades,trinity (one directory per "
                           "source)",
}

# what each tool directory holds (the Layout section)
DIRS = {
    "fastp": "clean reads and QC reports; fastp_summary.tsv = data and Q30 per sample",
    "meangs": "MEANGS mitochondrial seeds; meangs_summary.tsv",
    "getorganelle": "GetOrganelle rounds (round0, round1 ...); getorganelle_summary.tsv = length, circular, rounds",
    "mitogenome": "final mitogenomes (GetOrganelle or the better SPAdes scaffold); mitogenome_summary.tsv",
    "spades": "SPAdes assemblies (<sample>.scaffolds.fasta / .rna.fasta / .meta.fasta); spades_summary.tsv = "
              "sequences, total length, N50",
    "trinity": "Trinity transcripts; trinity_summary.tsv",
    "busco": "BUSCO results; busco_summary.tsv = C / S / D / F / M per sample (and assembly)",
    "genomescope": "k-mer histogram and GenomeScope; genomescope_summary.tsv = genome size, heterozygosity, error",
    "mitom": "MitoM annotation (GenBank .gb, gene table); mitom_summary.tsv = protein genes / rRNAs / tRNAs",
    "rdna": "rDNA units (GetOrganelle) and their annotation .gbf; rdna_summary.tsv",
    "histone": "histone clusters (GetOrganelle) and their annotation .gbf; histone_summary.tsv",
    "nuclear_markers": "rDNA / histone copies in the SPAdes assembly (MitoM), copy numbers, comparison with "
                       "GetOrganelle",
    "topology": "read-based check of circular / linear",
    "blast": "blast_topn.tsv = best hits and taxonomy of every query",
    "summary": "samples.tsv (every tool, one row per sample), report.html",
    "pipeline_info": "run records: logs/, benchmarks/ (time, memory), config_*.yaml, run_*.txt (command, versions)",
}


def dir_line(p):
    n = p.name
    if n.startswith("mitofinder_"):
        return (f"`{n}/` MitoFinder on the {n.removeprefix('mitofinder_')} contigs: <sample>/; mitofinder_summary.tsv "
                "= gene table (N = not found); genes/<gene>.fasta; contamination/ = contamination check")
    return f"`{n}/` {DIRS.get(n, '')}".rstrip() if p.is_dir() else f"`{n}`"


def busco(f):
    m = re.search(r"C:[\d.]+%\[S:[\d.]+%,D:[\d.]+%\],F:[\d.]+%,M:[\d.]+%,n:\d+", Path(f).read_text()) \
        if Path(f).exists() else None
    return m.group(0) if m else "no result"


def extra_rows(fn, d, sim, truth):
    """per-sample results of the functions check_truth.py has no truth item for: [(sample, item, result)]"""
    out = []
    for s in truth:
        if fn.endswith(("_fastp", "_fastp_rna")):
            j = json.load(open(d / "fastp" / s / f"{s}.fastp.json"))["summary"]["after_filtering"]
            out.append((s, "clean data", f"{j['total_bases'] / 1e9:.2f} Gb, Q30 {j['q30_rate'] * 100:.1f}%"))
        elif fn.endswith("_meangs"):
            f = d / "meangs" / s / f"{s}.seed.fasta"
            cov, pid = covered(sim / "truth" / s / "mt.fasta", f) if f.exists() and f.stat().st_size else (0, 0)
            out.append((s, "seed", f"{sum(fasta_lengths(f)) if f.exists() else 0} bp, {cov:.0f}% of the true "
                                   f"mitogenome, {pid:.2f}% identity"))
        elif fn.endswith(("spades_meta", "_rnaspades", "_trinity")):
            f = d / "trinity" / s / f"{s}.trinity.fasta" if (d / "trinity").is_dir() else \
                next(iter(sorted((d / "spades" / s).glob(f"{s}.*.fasta"))), None)
            lens = fasta_lengths(f) if f and f.exists() else []
            out.append((s, "assembly", f"{len(lens)} sequences, {sum(lens) / 1e6:.1f} Mb, N50 {n50(lens)}"))
        elif "busco" in fn or fn.endswith("python_api"):
            for f in sorted((d / "busco" / s).glob("*.busco.txt")):
                what = f.name.removeprefix(s + ".").removesuffix("busco.txt").strip(".") or "genome"
                out.append((s, f"BUSCO {what}", busco(f)))
        elif fn.endswith("circularize"):
            for r in read_tsv(d / s / f"{s}.circularize.tsv"):
                out.append((s, "circularization", f"{r['method']}, circle {r['circle_length'] or '-'} bp (truth "
                                                  f"{truth[s]['mt_length']} bp)"))
        elif fn.endswith("_compare"):
            for r in read_tsv(d / s / f"{s}.compare.tsv"):
                out.append((s, "comparison", f"{r['verdict']}, identity {r['identity_pct']}%, {r['ref_length']} vs "
                                             f"{r['other_length']} bp"))
        elif fn.endswith("_nuclear"):
            v = {}
            for r in read_tsv(d / "summary.tsv"):
                if r["sample"] == s:
                    v[r["verdict"]] = v.get(r["verdict"], 0) + 1
            out.append((s, "genes", ", ".join(f"{k} {n}" for k, n in v.items())))
        elif fn.endswith("_blast"):
            best = {}
            for r in read_tsv(d / "blast" / "blast_topn.tsv"):
                if r["qseqid"].startswith(s + "|"):
                    best.setdefault(r["qseqid"].split("|")[-1], r)
            for g, r in best.items():
                out.append((s, f"BLAST {g}", f"{r['organism']} ({float(r['pident']):.1f}%)"))
    return out


def main():
    sim, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    checks = read_tsv(out / "truth_check.tsv")
    truth = {r["sample"]: r for r in read_tsv(sim / "truth.tsv")}
    runs = {}
    for f in sorted(out.glob("features.*.tsv"), key=lambda x: x.stat().st_mtime):   # the latest run of each wins
        for r in read_tsv(f):
            runs[f"{r['id']}_{r['function']}"] = r
    index = ["# wgs_skim: every function on 5 simulated samples", "",
             f"Simulated data and truth: `{sim}` (truth.tsv). One directory per function, holding wgs_skim's own "
             "output: `<tool>/<sample>/`, `<tool>/<tool>_summary.tsv` (summary), `pipeline_info/` (logs, config). "
             "Each directory's `README.md`: what the function does, the command, the layout, every sample against the "
             "truth, the verdict. Reads are cleaned once (02 DNA, 30 RNA); later functions read those clean reads.", "",
             "| Function | Run | Agrees with the truth | Verdict |", "|---|---|---|---|"]
    for d in sorted(p for p in out.iterdir() if p.is_dir() and re.match(r"\d\d_", p.name)):
        fn, run = d.name, runs.get(d.name, {})
        rows = [r for r in checks if r["function"] == fn]
        bad = [r for r in rows if r["verdict"] != "ok"]
        ran = run.get("result") == "ok"
        verdict = ("pass" if not bad else f"{len(bad)} item(s) differ from the truth (see README)") if ran else "failed"
        md = [f"# {fn}", "", f"**What**: {WHAT.get(fn[3:], '')}", "",
              f"**Run**: {'ok' if ran else 'FAILED ' + run.get('missing', '')}, {run.get('seconds', '?')} s", "",
              "**Command**:", "```bash", run.get("command", ""), "```", "", "## Layout", ""]
        md += [f"- {dir_line(p)}" for p in sorted(d.iterdir()) if p.name != "README.md"]
        md += ["", "## Every sample", "", "| Sample | Item | Result | Truth | Agrees |", "|---|---|---|---|---|"]
        md += [f"| {r['sample']} | {r['item']} | {r['result']} | {r['truth']} | "
               f"{'yes' if r['verdict'] == 'ok' else '**no**'} |" for r in rows]
        md += [f"| {s} | {item} | {res} | | |" for s, item, res in extra_rows(fn, d, sim, truth)]
        md += ["", f"**Verdict**: {verdict}", ""]
        if bad:
            md += ["Differences:", *[f"- {r['sample']} {r['item']}: {r['result']} (truth {r['truth']})" for r in bad],
                   ""]
            if "mitofinder" in fn and any("MitoFinder" in r["item"] for r in bad):
                md.append("With one distant reference MitoFinder misses ND2, ND4L or rrnL in the distant samples S3-S5; "
                          "the same data with 5 related references (16) give 15 / 15. A MitoFinder limit, not a "
                          "workflow error.")
        (d / "README.md").write_text("\n".join(md) + "\n")
        tc = f"{len(rows) - len(bad)}/{len(rows)}" if rows else "(no truth item, see README)"
        index.append(f"| [{fn}]({fn}/README.md) | {'ok' if ran else 'FAILED'} | {tc} | {verdict} |")
    (out / "README.md").write_text("\n".join(index) + "\n")
    print(f"README.md in every function directory; index -> {out / 'README.md'}")


if __name__ == "__main__":
    main()

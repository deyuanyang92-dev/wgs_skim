"""Copy number of rDNA and histone genes from assembly coverage (rule nuc_copies).

copies = median k-mer coverage of the contigs carrying the gene / median coverage of the contigs carrying BUSCO
complete single-copy genes. Collapsed tandem copies (rDNA arrays, histone clusters) assemble into one contig whose
coverage is the sum of all copies, so the ratio estimates the copy number; contigs below 20% of the gene's best
coverage (fragments, contamination) are ignored.
"""
import csv
import statistics
import sys
from pathlib import Path

sm = snakemake  # noqa: F821
sys.stderr = open(sm.log[0], "w")


def cov(contig):
    try:
        return float(contig.rsplit("_cov_", 1)[1])
    except (IndexError, ValueError):
        return None


single = []
for f in Path(sm.params.busco_dir).glob("run_*/full_table.tsv"):
    for line in open(f):
        x = line.rstrip("\n").split("\t")
        if not line.startswith("#") and len(x) > 2 and x[1] == "Complete" and cov(x[2]) is not None:
            single.append(cov(x[2]))
base = statistics.median(single) if single else None
print(f"{len(single)} single-copy BUSCO contigs, median coverage {base}", file=sys.stderr)

genes = {}
if Path(sm.input.features).stat().st_size:
    for r in csv.DictReader(open(sm.input.features), delimiter="\t"):
        c = cov(r["contig"]) or (float(r["contig_cov"]) if r.get("contig_cov") else None)
        if c is not None:
            genes.setdefault(r["gene"], []).append(c)
cons = {}
if Path(sm.input.consensus).stat().st_size:
    cons = {r["gene"]: r for r in csv.DictReader(open(sm.input.consensus), delimiter="\t")}

with open(sm.output.tsv, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gene", "consensus_length", "partial", "contigs_with_gene", "coverage", "single_copy_coverage",
                "copy_number"])
    for g in ["18S", "ITS1", "5.8S", "ITS2", "28S", "H1", "H2A", "H2B", "H3", "H4"]:
        if g not in genes and g not in cons:
            continue
        covs = genes.get(g, [])
        good = [c for c in covs if c >= 0.2 * max(covs)] if covs else []
        c = statistics.median(good) if good else None
        r = cons.get(g, {})
        partial = "5'" * (r.get("partial_5") == "1") + "3'" * (r.get("partial_3") == "1")
        w.writerow([g, r.get("length", ""), partial, len(good), f"{c:.1f}" if c else "",
                    f"{base:.2f}" if base else "",
                    ("<1 (fragments)" if c / base < 0.5 else f"{c / base:.0f}") if c and base else ""])

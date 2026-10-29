# GetOrganelle vs other mitogenome assemblies (SPAdes-MitoFinder ...) for every sample of an existing results
# directory or a sample sheet - `wgs_skim compare --results / --sheet` runs it; Snakemake schedules by cores and
# memory, logs each sample and resumes. Per sample: scripts/compare_mt.py (circularize, align, structure, genes,
# figure); then all.*.tsv + all.compare.html.
#   snakemake -s workflow/compare_batch.smk --cores 16 --config outdir=/path/compare results=/path/results
#       [sheet=samples.tsv] [samples=A,B] [plot=true] tools=<dir with mafft + minimap2>
import sys
from pathlib import Path

sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
from compare_mt import jobs_from_results, jobs_from_sheet, missing_inputs  # noqa: E402

OUT = Path(config["outdir"]).resolve()
CMP = str(Path(workflow.basedir) / "scripts" / "compare_mt.py")
TOOLS = config.get("tools", "")
JOBS = jobs_from_results(config["results"]) if config.get("results") else jobs_from_sheet(config["sheet"])
if config.get("samples"):
    JOBS = [j for j in JOBS if j["sample"] in str(config["samples"]).split(",")]
SKIPPED, RUN = [], {}
for j in JOBS:
    miss = missing_inputs(j)                      # drops absent optional files, lists absent required ones
    if miss:
        SKIPPED.append(f"{j['sample']}:missing: {', '.join(miss)}")
    else:
        RUN[j["sample"]] = j


wildcard_constraints:
    sample="|".join(RUN) or "none",


rule all:
    input:
        OUT / "all.compare.tsv",


def args(w):
    j = RUN[w.sample]
    a = [f"-a {j['ref']}"] + [f"-b {o}" for o in j["other"]]
    a += [f"--{k} {j[k]}" for k in ("spades", "gb", "gff", "gfa") if j.get(k)]
    if j.get("graph_for"):
        a.append(f"--graph-for {j['graph_for']}")
    if str(config.get("plot", "true")).lower() in ("false", "no", "0"):
        a.append("--no-plot")
    return " ".join(a)


rule compare_sample:
    """one sample: every other assembly against its GetOrganelle result"""
    input:
        ref=lambda w: RUN[w.sample]["ref"],
    output:
        tsv=OUT / "{sample}.compare.tsv",
    params:
        args=args,
        path=f"export PATH={TOOLS}:$PATH; " if TOOLS else "",
    log:
        OUT / "logs" / "{sample}.log",
    benchmark:
        OUT / "benchmarks" / "{sample}.tsv"
    threads: 4
    resources:
        mem_mb=8000,                              # the SPAdes assembly graph of the circularization (~7 GB)
    shell:
        "{params.path}python {CMP} {params.args} -o {OUT} --prefix {wildcards.sample} > {log} 2>&1"


rule merge:
    input:
        expand(OUT / "{sample}.compare.tsv", sample=list(RUN)),
    output:
        OUT / "all.compare.tsv",
    params:
        samples=",".join(RUN),
        skipped=" ".join(f"'{x}'" for x in SKIPPED),
    log:
        OUT / "logs" / "merge.log",
    shell:
        "python {CMP} --merge '{params.samples}' --skipped {params.skipped} -o {OUT} > {log} 2>&1"

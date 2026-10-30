# Multi-copy nuclear markers (rDNA, histones) of every sample of an EXISTING results directory - one made by this
# workflow or by Temp_scripts/wgs_skim_pipeline (layouts: scripts/nuclear_markers.py LAYOUTS). `wgs_skim nuclear`
# runs it; Snakemake schedules the jobs by cores and memory, logs each one and resumes after an interruption.
#   snakemake -s workflow/nuclear_batch.smk --sdm conda --cores 64 --resources mem_mb=200000 \
#       --config results=/path/results outdir=/path/nuclear mitom=/path/Mitom [markers=rdna,histone max_gb=5,10]
# Per sample: GetOrganelle per marker (rounds of 5 / 10 Gb; an rDNA assembly the results directory has is round 1)
# -> MitoM per marker (.gbf) -> MitoM once on the SPAdes scaffolds (both markers) -> GetOrganelle vs SPAdes per
# gene -> summary.tsv. The steps are the functions of scripts/nuclear_markers.py, as in rules/nuclear.smk.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
from nuclear_markers import GENES, layout, result_fasta  # noqa: E402

OUT = Path(config["outdir"]).resolve()
MARKERS = [m for m in str(config.get("markers", "rdna,histone")).split(",") if m in GENES]
MAX_GB = str(config.get("max_gb", "5,10"))
MITOM = config["mitom"]
NM = str(Path(workflow.basedir) / "scripts" / "nuclear_markers.py")
LAYOUT, PATHS = layout(config["results"])
if config.get("samples"):
    PATHS = {s: p for s, p in PATHS.items() if s in str(config["samples"]).split(",")}
SAMPLES = sorted(PATHS)
print(f"[nuclear] {config['results']}: layout {LAYOUT}, {len(SAMPLES)} sample(s), markers {','.join(MARKERS)}",
      file=sys.stderr)


wildcard_constraints:
    sample="|".join(SAMPLES),
    kind="|".join(GENES),


def scaffolds(s):
    f = PATHS[s]["scaffolds"]
    return f if f.exists() and f.stat().st_size else None


def spades_mitom(s):
    """MitoM on the SPAdes scaffolds: the workflow's own run (rule nuc_markers), one of an earlier batch, or ours"""
    given = PATHS[s].get("spades_mitom")
    if given and Path(given, "consensus.tsv").exists():
        return Path(given)
    for d in [OUT / "spades_mitom" / s] + [OUT / k / "spades_mitom" / s for k in MARKERS]:
        if (d / "mitom.version").exists() and (d / "consensus.tsv").exists():
            return d
    return OUT / "spades_mitom" / s if scaffolds(s) else None


rule all:
    input:
        OUT / "summary.tsv",


rule assemble:
    """GetOrganelle -F anonym of one multi-copy unit (rounds of max_gb Gb of reads)"""
    input:
        reads=lambda w: [PATHS[w.sample]["r1"], PATHS[w.sample]["r2"]] if PATHS[w.sample]["r1"].exists() else [],
    output:
        done=OUT / "{kind}" / "getorganelle" / "{sample}" / "done",
    params:
        given=lambda w: f"--given {PATHS[w.sample][w.kind]}" if PATHS[w.sample].get(w.kind) and
        result_fasta(PATHS[w.sample][w.kind]) else "",
        reads=lambda w, input: f"-1 {input.reads[0]} -2 {input.reads[1]}" if input.reads else "",
        out=lambda w, output: str(Path(output.done).parent),
    log:
        OUT / "pipeline_info" / "logs" / "assemble" / "{kind}" / "{sample}.log",
    benchmark:
        OUT / "pipeline_info" / "benchmarks" / "assemble" / "{kind}" / "{sample}.tsv"
    threads: 16
    resources:
        mem_mb=int(max(8, 1.2 * max(float(x) for x in MAX_GB.split(",") if x)) * 1024),
    conda:
        "envs/getorganelle.yaml"
    shell:
        "python {NM} assemble --kind {wildcards.kind} {params.reads} {params.given} -o {params.out} "
        "-t {threads} --max-gb {MAX_GB} > {log} 2>&1"


rule annotate:
    """MitoM on the GetOrganelle unit -> <sample>.<kind>.gbf"""
    input:
        done=rules.assemble.output.done,
    output:
        stamp=OUT / "{kind}" / "annotation" / "{sample}" / "mitom.version",
    params:
        go=lambda w, input: str(Path(input.done).parent),
        out=lambda w, output: str(Path(output.stamp).parent),
    log:
        OUT / "pipeline_info" / "logs" / "annotate" / "{kind}" / "{sample}.log",
    benchmark:
        OUT / "pipeline_info" / "benchmarks" / "annotate" / "{kind}" / "{sample}.tsv"
    threads: 4
    resources:
        mem_mb=2000,
    conda:
        "envs/mitom.yaml"
    shell:
        "python {NM} annotate -i {params.go} -o {params.out} --kind {wildcards.kind} --sample {wildcards.sample} "
        "--mitom {MITOM} -t {threads} > {log} 2>&1"


rule spades_mitom:
    """MitoM once on the whole SPAdes assembly, both markers (hours on a 2 Gb assembly)"""
    input:
        scaffolds=lambda w: scaffolds(w.sample),
    output:
        stamp=OUT / "spades_mitom" / "{sample}" / "mitom.version",
    params:
        out=lambda w, output: str(Path(output.stamp).parent),
        kinds=",".join(MARKERS),
    log:
        OUT / "pipeline_info" / "logs" / "spades_mitom" / "{sample}.log",
    benchmark:
        OUT / "pipeline_info" / "benchmarks" / "spades_mitom" / "{sample}.tsv"
    threads: 12
    resources:
        mem_mb=8000,
    conda:
        "envs/mitom.yaml"
    shell:
        "python {NM} annotate -i {input.scaffolds} -o {params.out} --kind {params.kinds} --sample {wildcards.sample} "
        "--mitom {MITOM} -t {threads} --no-gbf > {log} 2>&1"


def compare_inputs(w):
    d = {k: str(OUT / k / "annotation" / w.sample / "mitom.version") for k in MARKERS}
    sp = spades_mitom(w.sample)
    if sp is not None and sp == OUT / "spades_mitom" / w.sample:
        d["spades"] = str(sp / "mitom.version")           # ours: a job; an existing run is just read
    return d


rule compare:
    """GetOrganelle vs SPAdes-MitoM, per gene"""
    input:
        unpack(compare_inputs),
    output:
        tsv=OUT / "compare" / "{sample}.tsv",
    params:
        go=lambda w: " ".join(f"--getorganelle-{k} {OUT / k / 'annotation' / w.sample}" for k in MARKERS),
        spades=lambda w: f"--spades {spades_mitom(w.sample)}" if spades_mitom(w.sample) else "",
        markers=",".join(MARKERS),
    log:
        OUT / "pipeline_info" / "logs" / "compare" / "{sample}.log",
    resources:
        mem_mb=1000,
    conda:
        "envs/mitom.yaml"
    shell:
        "python {NM} compare --sample {wildcards.sample} {params.go} {params.spades} --markers {params.markers} "
        "-o {output.tsv} > {log} 2>&1"


rule summary:
    """every sample x gene in one table (+ README.txt)"""
    input:
        expand(OUT / "compare" / "{sample}.tsv", sample=SAMPLES),
    output:
        OUT / "summary.tsv",
    log:
        OUT / "pipeline_info" / "logs" / "summary.log",
    resources:
        mem_mb=1000,
    shell:
        "python {NM} merge {input} -o {output} > {log} 2>&1"

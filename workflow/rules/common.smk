import json
import re
import sys
from pathlib import Path

from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
from samples import find_pairs, read_sheet  # noqa: E402

# ── Taxon presets (YZZ Team groups); explicit config keys win ────────────────
TAXA = {
    "annelid": dict(meangs_class="A-worms", genetic_code=5, mitoz_clade="Annelida-segmented-worms", busco_lineage="metazoa_odb12"),
    "mollusc": dict(meangs_class="Mollusca", genetic_code=5, mitoz_clade="Mollusca", busco_lineage="mollusca_odb12"),
    "fish": dict(meangs_class="Chordata", genetic_code=2, mitoz_clade="Chordata", busco_lineage="actinopterygii_odb12"),
    "crustacean": dict(meangs_class="Arthropoda", genetic_code=5, mitoz_clade="Arthropoda", busco_lineage="arthropoda_odb12"),
}
for key, value in TAXA[config["taxon"]].items():
    if config.get(key) in (None, ""):
        config[key] = value

# ── Samples: sample sheet or automatic detection in reads_dir ────────────────
if config.get("samples"):
    SAMPLES = read_sheet(config["samples"])
elif config.get("reads_dir"):
    SAMPLES = find_pairs(config["reads_dir"])
else:
    raise WorkflowError("give the reads directory: --config reads_dir=/path/to/rawdata (or samples=sheet.tsv)")

OUT = config["outdir"].rstrip("/")
RES = Path(workflow.basedir).parent / "resources"           # bundled files (read-only is fine)
DB = Path(config.get("db_dir") or RES).expanduser()            # downloaded databases (GetOrganelle, BUSCO)
NR_SEED = config.get("nr_seed") or str(RES / "rdna" / "seed.fasta")
NR_LABEL = config.get("nr_label") or str(RES / "rdna" / "label.fasta")

wildcard_constraints:
    sample="|".join(re.escape(s) for s in SAMPLES),


def reads(wildcards):
    r1, r2 = SAMPLES[wildcards.sample]
    return {"r1": r1, "r2": r2}


def clean_gb(json_file):
    """Clean data (Gb) from fastp.json; 0 if not available yet."""
    try:
        return json.load(open(json_file))["summary"]["after_filtering"]["total_bases"] / 1e9
    except (OSError, KeyError, ValueError):
        return 0.0


def spades_mem_mb(wildcards, input, attempt):
    """SPAdes -m is a hard limit: estimate from clean data and grow 1.5x on each retry."""
    est = clean_gb(input.json) * config["spades_gib_per_gb"]
    return int(max(config["spades_min_mem_gb"], est) * 1.5 ** (attempt - 1) * 1024)


def final_outputs():
    s = list(SAMPLES)
    out = [f"{OUT}/summary/samples.tsv"]
    if config["mito"]:
        out += expand(f"{OUT}/mito/{{sample}}/{{sample}}.mt.fasta", sample=s)
    if config["nr"]:
        out += expand(f"{OUT}/rdna/{{sample}}/{{sample}}.nr.fasta", sample=s)
    return out


def busco_lineage_dir(wildcards):
    lin = str(config["busco_lineage"])
    return lin if Path(lin).expanduser().is_dir() else str(DB / "busco_downloads" / "lineages" / lin)


def blast_inputs(wildcards):
    s = list(SAMPLES)
    fas = []
    if config["mito"]:
        fas += expand(rules.mt_iterate.output.fasta, sample=s)
    if config["nr"]:
        fas += expand(rules.getorganelle_nr.output.fasta, sample=s)
    return fas


def summary_inputs(wildcards):
    s = list(SAMPLES)
    d = {"fastp": expand(rules.fastp.output.json, sample=s)}
    if config["mito"]:
        d["mt"] = expand(rules.mt_iterate.output.status, sample=s)
        d["mt_fasta"] = expand(rules.mt_iterate.output.fasta, sample=s)
    if config["nr"]:
        d["nr"] = expand(rules.getorganelle_nr.output.fasta, sample=s)
    if config["assembly"]:
        d["scaffolds"] = expand(rules.spades.output.scaffolds, sample=s)
        d["busco"] = expand(rules.busco.output.summary, sample=s)
    if config["genome_size"]:
        d["genomescope"] = expand(rules.genomescope.output.summary, sample=s)
    if config["annotate"] and config["mito"]:
        d["annotation"] = expand(rules.mitoz_annotate.output.gb, sample=s)
    if config["topology"] and config["mito"]:
        d["topology"] = expand(rules.topology.output.tsv, sample=s)
    if config["blast"] and (config["mito"] or config["nr"]):
        d["blast"] = rules.blast.output.tsv
    return d

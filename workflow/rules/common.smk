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
    "annelid": dict(meangs_class="A-worms", genetic_code=5, clade="Annelida-segmented-worms", busco_lineage="metazoa_odb12"),
    "mollusc": dict(meangs_class="Mollusca", genetic_code=5, clade="Mollusca", busco_lineage="mollusca_odb12"),
    "fish": dict(meangs_class="Chordata", genetic_code=2, clade="Chordata", busco_lineage="actinopterygii_odb12"),
    "crustacean": dict(meangs_class="Arthropoda", genetic_code=5, clade="Arthropoda", busco_lineage="arthropoda_odb12"),
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

# ── Steps: `steps` runs only some parts (e.g. "fastp,spades,busco"); empty = everything switched on ──────
STEP_ALIASES = {"qc": "fastp", "meangs-getorganelle": "mito", "mitogenome": "mito", "rdna": "nr", "assembly": "spades",
                "jellyfish": "genome_size", "genomescope": "genome_size", "mitoz": "annotate", "mitom": "annotate",
                "nuc": "nuclear", "histone": "nuclear", "histones": "nuclear", "report": "summary"}
STEPS_ALL = ["fastp", "meangs", "mito", "getorganelle", "nr", "spades", "busco", "genome_size", "annotate", "topology",
             "blast", "nuclear", "summary"]
STEPS = [STEP_ALIASES.get(x, x) for x in re.split(r"[,\s]+", str(config.get("steps") or "").strip()) if x]
for x in STEPS:
    if x not in STEPS_ALL:
        raise WorkflowError(f"unknown step '{x}'; choose from: {', '.join(STEPS_ALL)}")
if STEPS:
    want = set(STEPS)
    genes_wanted = any(x.strip() not in ("", "full") for x in f"{config['blast_mt']},{config['blast_nr']}".split(","))
    config["annotate"] = bool(want & {"annotate", "topology"}) or ("blast" in want and genes_wanted)
    config["topology"] = "yes" if "topology" in want else "no"
    config["blast"] = "blast" in want
    config["mito"] = bool(want & {"mito", "getorganelle", "annotate", "topology", "blast"})
    config["nr"] = bool(want & {"nr", "blast"})
    config["assembly"] = bool(want & {"spades", "busco", "nuclear"})
    config["nuc_markers"] = "nuclear" in want
    config["genome_size"] = "genome_size" in want
    if "getorganelle" in want and not want & {"mito", "meangs"}:
        config["meangs"] = False                                     # GetOrganelle with its own animal_mt seeds

# BLAST queries: "full" (the whole assembled sequence) and/or genes cut from the annotations
BLAST_MT = [x.strip() for x in str(config["blast_mt"]).split(",") if x.strip()]
BLAST_NR = [x.strip() for x in str(config["blast_nr"]).split(",") if x.strip()]

# topology: auto (only linear mitogenomes) | yes | no; YAML yes/no/true/false arrive as booleans
TOPOLOGY = {True: "yes", False: "no"}.get(config["topology"], str(config["topology"]))
config["topology"] = TOPOLOGY != "no"

OUT = "" if config["outdir"].strip("./") == "" else config["outdir"].rstrip("/") + "/"   # path prefix; "./x" makes Snakemake warn
RES = Path(workflow.basedir).parent / "resources"           # bundled files (read-only is fine)
DB = Path(config.get("db_dir") or RES).expanduser()            # downloaded databases (GetOrganelle, BUSCO)
MITOM_DB = DB / "mitom"
MITOM_NUC_DB = DB / "mitom_nuc"
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
    if STEPS:
        targets = {
            "fastp": [rules.fastp.output.json],
            "meangs": [rules.meangs.output.seed],
            "mito": [rules.mt_select.output.fasta],
            "getorganelle": [rules.mt_iterate.output.fasta],
            "nr": [rules.getorganelle_nr.output.fasta],
            "spades": [rules.spades.output.scaffolds],
            "busco": [rules.busco.output.summary],
            "genome_size": [rules.genomescope.output.summary],
            "annotate": [rules.mitom_annotate.output.gb],
            "topology": [rules.topology.output.tsv],
            "nuclear": [rules.nuc_copies.output.tsv],
        }
        out = [f for x in STEPS for t in targets.get(x, []) for f in expand(t, sample=s)]
        if "blast" in STEPS:
            out.append(rules.blast.output.tsv)
        if "summary" in STEPS:
            out.append(rules.summary.output.tsv)
        return out
    out = [f"{OUT}summary/samples.tsv"]
    if config["mito"]:
        out += expand(f"{OUT}mito/{{sample}}/{{sample}}.mt.fasta", sample=s)
    if config["nr"]:
        out += expand(f"{OUT}rdna/{{sample}}/{{sample}}.nr.fasta", sample=s)
    return out


def mt_seed(wildcards):
    """MEANGS seed, or an empty file (GetOrganelle then uses its animal_mt seeds) when meangs is switched off."""
    if config["meangs"]:
        return rules.meangs.output.seed.format(sample=wildcards.sample)
    return str(RES / "no_seed.fasta")


def mt_select_inputs(wildcards):
    d = {"fasta": rules.mt_iterate.output.fasta, "status": rules.mt_iterate.output.status, "seed": mt_seed(wildcards)}
    if config["assembly"] and config["mt_from_spades"]:
        d["scaffolds"] = rules.spades.output.scaffolds
    return d


def busco_lineage_dir(wildcards):
    lin = str(config["busco_lineage"])
    return lin if Path(lin).expanduser().is_dir() else str(DB / "busco_downloads" / "lineages" / lin)


def blast_inputs(wildcards):
    """mt: COX1 from the MitoM annotation (the barcode; whole mitogenome if not annotated); nr: the rDNA."""
    s = list(SAMPLES)
    d = {"mt": [], "gb": [], "nr": []}
    if config["mito"]:
        d["mt"] = expand(rules.mt_select.output.fasta, sample=s)
        if config["annotate"] and any(x != "full" for x in BLAST_MT):
            d["gb"] = expand(rules.mitom_annotate.output.gb, sample=s)
    if config["nr"]:
        d["nr"] = expand(rules.getorganelle_nr.output.fasta, sample=s)
        if config["annotate"] and any(x != "full" for x in BLAST_NR):
            d["nr_gb"] = expand(rules.nr_annotate.output.gb, sample=s)
    return d


def summary_inputs(wildcards):
    s = list(SAMPLES)
    d = {"fastp": expand(rules.fastp.output.json, sample=s)}
    if config["mito"]:
        d["mt"] = expand(rules.mt_select.output.status, sample=s)
        d["mt_fasta"] = expand(rules.mt_select.output.fasta, sample=s)
    if config["nr"]:
        d["nr"] = expand(rules.getorganelle_nr.output.fasta, sample=s)
    if config["assembly"]:
        d["scaffolds"] = expand(rules.spades.output.scaffolds, sample=s)
        d["busco"] = expand(rules.busco.output.summary, sample=s)
    if config["assembly"] and config["nuc_markers"]:
        d["nuclear"] = expand(rules.nuc_copies.output.tsv, sample=s)
    if config["genome_size"]:
        d["genomescope"] = expand(rules.genomescope.output.summary, sample=s)
    if config["annotate"] and config["mito"]:
        d["annotation"] = expand(rules.mitom_annotate.output.gb, sample=s)
    if config["topology"] and config["mito"]:
        d["topology"] = expand(rules.topology.output.tsv, sample=s)
    if config["blast"] and (config["mito"] or config["nr"]):
        d["blast"] = rules.blast.output.tsv
    return d

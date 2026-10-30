import json
import os
import re
import sys
from pathlib import Path

from snakemake.utils import validate

# `--config key=` (an empty value, e.g. mitofinder_sif= for the bioconda MitoFinder) arrives as None: it means ""
for _k, _v in list(config.items()):
    if _v is None:
        config[_k] = ""
validate(config, schema="../schemas/config.schema.yaml")
PY = sys.executable           # Snakemake's Python 3, for standard-library scripts run inside older tool environments
if not config.get("mitofinder_sif") and config.get("mitofinder_trna") == "mitfi":
    # bioconda MitoFinder has no MiTFi (envs/mitofinder.post-deploy.sh); tRNAscan-SE is also the first choice of
    # the user's method (Yang et al. 2024: tRNAscan-SE > MiTFi > ARWEN)
    print("[wgs_skim] MitoFinder without the Singularity image (mitofinder_sif): tRNAs with tRNAscan-SE, not MiTFi",
          file=sys.stderr)
    config["mitofinder_trna"] = "trnascan"

sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
from samples import find_pairs, read_sheet, sheet_columns  # noqa: E402

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
    GIVEN = sheet_columns(config["samples"], ["scaffolds", "mitogenome"])   # optional: existing files per sample
elif config.get("reads_dir"):
    SAMPLES = find_pairs(config["reads_dir"])
    GIVEN = {}
else:
    raise WorkflowError("give the reads directory: --config reads_dir=/path/to/rawdata (or samples=sheet.tsv)")
# samples with an assembly / a mitogenome made elsewhere: linked instead of assembled (rules assembly_given,
# mt_given), so BUSCO, MitoFinder, MitoM and BLAST run on them without reads (wgs_skim busco / mitofinder / mitom -i)
GIVEN_ASM = sorted(s for s, g in GIVEN.items() if "scaffolds" in g)
GIVEN_MT = sorted(s for s, g in GIVEN.items() if "mitogenome" in g)
NO_READS = sorted(s for s, (r1, r2) in SAMPLES.items() if not r1)


def one_of(names):
    """wildcard constraint matching exactly these samples (nothing when empty)"""
    return "|".join(re.escape(x) for x in names) or "(?!x)x"

# ── Steps: `steps` runs only some parts (e.g. "fastp,spades,busco"); empty = everything switched on ──────
STEP_ALIASES = {"qc": "fastp", "mitogenome": "mito", "rdna": "nr", "assembly": "spades",
                "jellyfish": "genome_size", "genomescope": "genome_size", "mitoz": "annotate", "mitom": "annotate",
                "nuc": "nuclear", "histone": "nr", "histones": "nr", "report": "summary",
                "rna": "transcriptome"}
# combinations (the recommended ones): name -> steps
COMBOS = {"meangs-getorganelle": ["mito"],                     # mitogenome, reseeded rounds until circular
          "spades-mitofinder": ["spades", "mitofinder"],       # mitogenome from the whole-genome assembly
          "rnaspades-mitofinder": ["rnaspades", "mitofinder"], # mt genes from RNA-seq (rnaSPAdes transcripts)
          "trinity-mitofinder": ["trinity", "mitofinder"],     # mt genes from RNA-seq (Trinity transcripts)
          "spades-busco": ["spades", "busco"],                 # whole genome
          "fastp-spades-busco": ["fastp", "spades", "busco"],
          "rnaspades-busco": ["rnaspades", "busco"],           # transcriptome (recommended)
          "trinity-busco": ["trinity", "busco"]}               # transcriptome
STEPS_ALL = ["fastp", "meangs", "mito", "getorganelle", "nr", "spades", "rnaspades", "busco", "genome_size",
             "annotate", "topology", "mitofinder", "blast", "nuclear", "trinity", "transcriptome", "summary"]
STEPS = []
for x in re.split(r"[,\s]+", str(config.get("steps") or "").strip()):
    if x:
        STEPS += COMBOS.get(x, [STEP_ALIASES.get(x, x)])
for x in STEPS:
    if x not in STEPS_ALL:
        raise WorkflowError(f"unknown step '{x}'; choose from: {', '.join(STEPS_ALL)}; "
                            f"combinations: {', '.join(COMBOS)}")
if STEPS:
    want = set(STEPS)
    genes_wanted = any(x.strip().lower() not in ("", "full", "none")
                       for x in f"{config['blast_mt']},{config['blast_nr']}".split(","))
    config["annotate"] = bool(want & {"annotate", "topology"}) or ("blast" in want and genes_wanted)
    config["topology"] = "yes" if "topology" in want else "no"
    config["blast"] = "blast" in want
    blast_mt, blast_nr = (str(config[k]).strip().lower() not in ("", "none") and "blast" in want
                          for k in ("blast_mt", "blast_nr"))
    config["mito"] = bool(want & {"mito", "getorganelle", "annotate", "topology"}) or blast_mt
    config["nr"] = bool(want & {"nr", "nuclear"}) or blast_nr
    config["assembly"] = bool(want & {"spades", "rnaspades", "busco", "nuclear"})
    config["mitofinder"] = "mitofinder" in want
    if config["mitofinder"] and config["mito"] and not config["mitofinder_ref"]:
        config["annotate"] = True                       # the sample's MitoM annotation is MitoFinder's reference
    config["nuc_markers"] = "nuclear" in want
    config["genome_size"] = "genome_size" in want
    if "getorganelle" in want and not want & {"mito", "meangs"}:
        config["meangs"] = False                                     # GetOrganelle with its own animal_mt seeds

# BLAST queries: "full" (the whole assembled sequence) and/or genes cut from the annotations
BLAST_MT = [x.strip() for x in str(config["blast_mt"]).split(",") if x.strip().lower() not in ("", "none")]
BLAST_NR = [x.strip() for x in str(config["blast_nr"]).split(",") if x.strip().lower() not in ("", "none")]

# topology: auto (only linear mitogenomes) | yes | no; YAML yes/no/true/false arrive as booleans
TOPOLOGY = {True: "yes", False: "no"}.get(config["topology"], str(config["topology"]))
config["topology"] = TOPOLOGY != "no"

OUT = os.path.join(os.path.normpath(config["outdir"]), "").removeprefix(f".{os.sep}")   # "./x" makes Snakemake warn
OUT_PATH = os.path.normpath(config["outdir"])
RES = Path(workflow.basedir).parent / "resources"           # bundled files (read-only is fine)
DB = Path(config.get("db_dir") or RES).expanduser()            # downloaded databases (GetOrganelle, BUSCO)
MITOM_DB = DB / "mitom"
MITOM_NUC_DB = DB / "mitom_nuc"
NR_SEED = config.get("nr_seed") or str(RES / "rdna" / "seed.fasta")
NR_LABEL = config.get("nr_label") or str(RES / "rdna" / "label.fasta")
HIST_SEED = config.get("histone_seed") or str(RES / "histone" / "seed.fasta")
HIST_LABEL = config.get("histone_label") or str(RES / "histone" / "label.fasta")
# SPAdes modes to run (genome = rule spades; others = rule spades_mode, transcriptome.smk); `--steps rna` = rnaSPAdes
MODES = [m.strip() for m in str(config["spades_modes"]).split(",") if m.strip()]
if STEPS and "transcriptome" in STEPS:
    MODES = list(dict.fromkeys([m for m in MODES if m != "genome"] + ["rna"]))
if STEPS and "rnaspades" in STEPS:                    # rnaSPAdes alone (the rna mode); `spades` too keeps genome
    MODES = list(dict.fromkeys((MODES if "spades" in STEPS else [m for m in MODES if m != "genome"]) + ["rna"]))
OTHER_MODES = [m for m in MODES if m != "genome"]

if STEPS and set(STEPS) == {"summary"}:
    def have(pattern):
        return all(Path(OUT_PATH, pattern.format(sample=s)).exists() for s in SAMPLES)
    config["mito"] = have("mitogenome/{sample}/{sample}.mt.fasta")
    config["nr"] = have("rdna/{sample}/{sample}.nr.fasta")
    config["assembly"] = have("busco/{sample}/{sample}.busco.txt")
    config["nuc_markers"] = have("nuclear_markers/{sample}/{sample}.copy_number.tsv")
    config["genome_size"] = have("genomescope/{sample}/{sample}.genomescope.txt")
    config["annotate"] = have("mitom/{sample}/{sample}.mt.gb")
    config["topology"] = "yes" if have("topology/{sample}/{sample}.topology.tsv") else "no"
    config["blast"] = Path(OUT_PATH, "blast/blast_topn.tsv").exists()

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
    if GIVEN.get(wildcards.sample, {}).get("scaffolds"):
        return 1000
    est = clean_gb(input.json) * config["spades_gib_per_gb"]
    return int(max(config["spades_min_mem_gb"], est) * 1.5 ** (attempt - 1) * 1024)


# Memory per job (resources: mem_mb) so that `snakemake --resources mem_mb=N` (wgs_skim run --mem-gb) schedules by
# memory as well as cores. Measured peaks (benchmarks/*/*.tsv max_rss, 0.6 Gb test sample 2026-10-28): fastp 0.4 GB,
# MEANGS 7.6 (fixed 2 M reads), GetOrganelle mt 3.2, rDNA 0.8, histones 1.6, SPAdes 6.9, BUSCO 2.0, jellyfish 6.8
# (-s 2G hash, independent of data), MitoM 0.2, mitom nuc on scaffolds 2.0, topology 0.9, BLAST 48.6 (nt pages mapped
# in memory; counted but shared). GetOrganelle grows with the reads it uses: ~1 GB per Gb (36 GB for 48 Gb, CHY006).
def getorganelle_mem_mb(gb_cap):
    """mem_mb for GetOrganelle rounds that use up to gb_cap Gb of reads (0 = all clean reads)"""
    def f(wildcards, input, attempt):
        gb = clean_gb(input.json)
        used = min(gb, gb_cap) if gb_cap else gb
        return int(max(8, 1.2 * used) * 1024 * attempt)
    return f


MT_ITER_MEM = getorganelle_mem_mb(config["mt_iter_max_data_gb"])
NUC_MEM = getorganelle_mem_mb(max(float(x) for x in str(config["nuclear_max_gb"]).split(",") if x.strip()))


MT_FASTA = object()             # targets marker: each sample's final mitogenome (mt_of)
SCAFFOLDS = object()            # targets marker: each sample's genome assembly (scaffolds_of)


def scaffolds_of(sample):
    """the genome assembly of a sample: SPAdes, or the one given in the sample sheet (linked under pipeline_info)"""
    return (rules.assembly_given if sample in GIVEN_ASM else rules.spades).output.scaffolds.format(sample=sample)


def mt_of(sample, what="fasta"):
    """the final mitogenome (fasta) or its record (status) of a sample: chosen by mt_select, or given"""
    return getattr((rules.mt_given if sample in GIVEN_MT else rules.mt_select).output, what).format(sample=sample)


def final_outputs():
    s = list(SAMPLES)
    if STEPS:
        targets = {
            "fastp": [rules.fastp.output.html],
            "meangs": [rules.meangs.output.seed],
            "mito": [MT_FASTA],
            "getorganelle": [MT_FASTA],
            "nr": [rules.nr_annotate.output.gb, rules.histone_annotate.output.gb],
            "spades": ([SCAFFOLDS] if "genome" in MODES else [])
            + [rules.spades_mode.output.fasta.replace("{mode}", m) for m in OTHER_MODES],
            "busco": ([rules.busco.output.summary] if "genome" in MODES else [])
            + [rules.busco_assembly.output.summary.replace("{asm}", m) for m in OTHER_MODES],
            "rnaspades": [rules.spades_mode.output.fasta.replace("{mode}", "rna")],
            "trinity": [rules.trinity.output.fasta],

            "genome_size": [rules.genomescope.output.summary],
            "annotate": [rules.mitom_annotate.output.gb],
            "topology": [rules.topology.output.tsv],
            "nuclear": [rules.nuc_copies.output.tsv, rules.nuc_compare.output.tsv],
            "transcriptome": [rules.busco_assembly.output.summary.replace("{asm}", a)
                              for a in [m for m in OTHER_MODES if m in ("rna", "rnaviral")] + ["trinity"]],
        }
        if "busco" in STEPS and "trinity" in STEPS:
            targets["busco"].append(rules.busco_assembly.output.summary.replace("{asm}", "trinity"))
        out = [f for x in STEPS for t in targets.get(x, []) for f in
               ([mt_of(y) for y in s] if t is MT_FASTA else [scaffolds_of(y) for y in s] if t is SCAFFOLDS
                else expand(t, sample=s))]
        if "mitofinder" in STEPS:
            out += mitofinder_targets(rules.mitofinder.output.fasta)
            srcs = sorted({x for s_ in SAMPLES for x in mitofinder_sources(s_)})
            out += [rules.mitofinder_summary.output.tsv.format(src=x) for x in srcs]
            if config["mitofinder_blast"]:
                out += [rules.mitofinder_blast.output.samples.format(src=x) for x in srcs]
        if "blast" in STEPS:
            out.append(rules.blast.output.tsv)
        if "summary" in STEPS:
            out.append(rules.summary.output.tsv)
        return out + tool_summary_targets()
    out = [f"{OUT}summary/samples.tsv"] + tool_summary_targets()
    if config["mitofinder"]:
        srcs = sorted({x for s_ in SAMPLES for x in mitofinder_sources(s_)})
        out += [rules.mitofinder_summary.output.tsv.format(src=x) for x in srcs]
        if config["mitofinder_blast"]:
            out += [rules.mitofinder_blast.output.samples.format(src=x) for x in srcs]
    if config["mito"]:
        out += [mt_of(x) for x in s]
    if config["assembly"]:
        out += expand([rules.busco_assembly.output.summary.replace("{asm}", m) for m in OTHER_MODES], sample=s)
    if config["nr"]:
        out += expand(f"{OUT}rdna/{{sample}}/{{sample}}.nr.fasta", sample=s)
        out += expand([rules.nr_annotate.output.gb, rules.histone_annotate.output.gb], sample=s)
        if config["assembly"] and config["nuc_markers"]:
            out += expand(rules.nuc_compare.output.tsv, sample=s)
    return out


def mt_seed(wildcards):
    """MEANGS seed, or an empty file (GetOrganelle then uses its animal_mt seeds) when meangs is switched off."""
    if config["meangs"]:
        return rules.meangs.output.seed.format(sample=wildcards.sample)
    return str(RES / "no_seed.fasta")


def mt_select_inputs(wildcards):
    d = {"fasta": rules.mt_iterate.output.fasta, "status": rules.mt_iterate.output.status, "seed": mt_seed(wildcards)}
    if config["assembly"] and config["mt_from_spades"]:
        d["scaffolds"] = scaffolds_of(wildcards.sample)
    return d


def busco_lineage_dir(wildcards):
    lin = str(config["busco_lineage"])
    return lin if Path(lin).expanduser().is_dir() else str(DB / "busco_downloads" / "lineages" / lin)


def blast_inputs(wildcards):
    """mt: COX1 from the MitoM annotation (the barcode; whole mitogenome if not annotated); nr: the rDNA."""
    s = list(SAMPLES)
    d = {"mt": [], "gb": [], "nr": []}
    if config["mito"]:
        d["mt"] = [mt_of(x) for x in s]
        if config["annotate"] and any(x != "full" for x in BLAST_MT):
            d["gb"] = expand(rules.mitom_annotate.output.gb, sample=s)
    if config["nr"]:
        d["nr"] = expand(rules.getorganelle_nr.output.fasta, sample=s)
        if config["annotate"] and any(x != "full" for x in BLAST_NR):
            d["nr_gb"] = expand(rules.nr_annotate.output.gb, sample=s)
    return d


def summary_inputs(wildcards):
    s = list(SAMPLES)
    d = {"fastp": [reads_json(x) for x in s if x not in NO_READS]}
    if config["mito"]:
        d["mt"] = [mt_of(x, "status") for x in s]
        d["mt_fasta"] = [mt_of(x) for x in s]
    if config["nr"]:
        d["nr"] = expand(rules.getorganelle_nr.output.fasta, sample=s)
    if config["assembly"]:
        d["scaffolds"] = [scaffolds_of(x) for x in s]
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


def tool_summaries():
    """{tool: per-sample files} of the tools this run uses: each gets <tool>/<tool>_summary.tsv (rule tool_summary)"""
    s, own = list(SAMPLES), [x for x in SAMPLES if x not in GIVEN_ASM]
    mt_own = [x for x in SAMPLES if x not in GIVEN_MT]
    steps = set(STEPS)
    t = {}
    with_reads = [x for x in s if x not in NO_READS]
    if config["qc"] and with_reads:
        t["fastp"] = expand(rules.fastp.output.json, sample=with_reads)
    if mt_own and ("meangs" in steps or (config["mito"] and config["meangs"])):     # MEANGS alone too
        t["meangs"] = expand(rules.meangs.output.seed, sample=mt_own)
    if config["mito"] and mt_own:
        t["getorganelle"] = expand(rules.mt_iterate.output.status, sample=mt_own)
        t["mitogenome"] = expand(rules.mt_select.output.status, sample=mt_own)
    if config["assembly"]:
        t["spades"] = (expand(rules.spades.output.scaffolds, sample=own) if "genome" in MODES else []) + \
            [f for m in OTHER_MODES for f in expand(rules.spades_mode.output.fasta, sample=s, mode=m)]
        if not STEPS or steps & {"busco", "transcriptome", "nuclear"}:
            t["busco"] = (expand(rules.busco.output.summary, sample=s) if "genome" in MODES else []) + \
                [f for m in OTHER_MODES for f in expand(rules.busco_assembly.output.summary, sample=s, asm=m)]
    if steps & {"trinity", "transcriptome"}:
        t["trinity"] = expand(rules.trinity.output.fasta, sample=s)
        if steps & {"busco", "transcriptome"}:
            t.setdefault("busco", []).extend(expand(rules.busco_assembly.output.summary, sample=s, asm="trinity"))
    if config["genome_size"]:
        t["genomescope"] = expand(rules.genomescope.output.summary, sample=s)
    if config["annotate"] and config["mito"]:
        t["mitom"] = expand(rules.mitom_annotate.output.gb, sample=s)
    if config["nr"]:
        t["rdna"] = expand(rules.getorganelle_nr.output.fasta, sample=s)
        t["histone"] = expand(rules.getorganelle_histone.output.fasta, sample=s)
    return {k: v for k, v in t.items() if v}


def tool_summary_targets():
    return [f"{OUT}{tool}/{tool}_summary.tsv" for tool in tool_summaries()]

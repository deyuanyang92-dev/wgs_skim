"""`wgs_skim <tool>`: one tool on its own, for every sample - the same rules as `wgs_skim run` (Snakemake schedules
the samples by cores and memory, one log per job, resumable). Reads tools take -r / -s and use the reads as they
are (--qc: fastp first); FASTA tools take -i (files, directories, or one merged FASTA with --split-by-id).

Adding a tool = one Tool(...) in TOOLS: its inputs, the workflow step(s) it runs, its own options (`options`) and
how they become workflow options (`configure`); the common options (-r/-s/-i, -o, -j, --mem-gb, -n ...) and the
batch handling are shared (docs/DEVELOPING.md)."""
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

from .paths import HOME, scripts_dir


@dataclass
class Tool:
    name: str
    inputs: str                      # "reads", "fasta" or "both" (-i assemblies or -r reads)
    help: str
    step: str = ""                   # workflow step(s) it runs (configure may choose)
    column: str = "mitogenome"       # sample-sheet column of the -i files: mitogenome or scaffolds
    options: Optional[Callable] = None      # options(argparse parser): the tool's own options
    configure: Optional[Callable] = None    # configure(args, conf, use_fasta) -> steps (or None: `step`)


# ── each tool's own options ────────────────────────────────────────────────
def _getorganelle_opts(p):
    p.add_argument("--target", choices=["mt", "nr"], default="mt",
                   help="mt = mitogenome; nr = nuclear rDNA (18S-ITS-28S) + histone cluster (5 Gb, then 10 Gb seeded "
                        "with round 1)")
    p.add_argument("--meangs", action="store_true", help="MEANGS seed first (= meangs-getorganelle)")
    p.add_argument("--rounds", type=int, help="reseeded rounds when not circular (0 = none; default 3)")


def _getorganelle_conf(a, conf, fasta):
    if a.rounds is not None:
        conf["mt_iter_rounds"] = a.rounds
    return "nr" if a.target == "nr" else ("mito" if a.meangs else "getorganelle")


def _spades_opts(p):
    p.add_argument("--mode", default="genome", metavar="M[,M]",
                   help="genome, isolate, sc, meta, rna, rnaviral, plasmid, metaplasmid, metaviral, bio, sewage, "
                        "corona (several: one assembly each)")


def _spades_conf(a, conf, fasta):
    conf["spades_modes"] = a.mode


def _busco_opts(p):
    p.add_argument("-m", "--mode", choices=["genome", "transcriptome", "proteins"], default="genome")
    p.add_argument("-l", "--lineage", help="lineage name or downloaded lineage directory (default from --taxon: "
                   "metazoa_odb12 for annelids)")


def _busco_conf(a, conf, fasta):
    conf.update(busco_mode=a.mode, busco_lineage=a.lineage)


MF_ASSEMBLERS = ("megahit", "metaspades", "spades", "rnaspades", "trinity")


def _mitofinder_opts(p):
    p.add_argument("--ref", nargs="+", metavar="GB|DIR", help="reference mitogenome(s) of related species, "
                   "GenBank (required): several files / a directory are merged - several related species find more "
                   "genes (sim5: 12-13 -> 15 of 15 with 5 related references)")
    p.add_argument("--assembler", default="megahit", metavar="A[,A]",
                   help="with -r: what assembles the reads for MitoFinder: megahit / metaspades (MitoFinder's own), "
                        "spades (SPAdes genome), rnaspades / trinity (transcripts: mt genes from RNA-seq); several = "
                        "one MitoFinder run each")
    p.add_argument("--trna", choices=["mitfi", "trnascan", "arwen"], default="mitfi")
    p.add_argument("--no-check", action="store_true", help="skip the contamination check (BLAST of every gene of "
                   "every sample against nt, genes hitting another family flagged)")
    p.add_argument("--db", help="BLAST database of the check (default: blast_db of ~/.wgs_skim/config.yaml, else "
                   "NCBI online)")


def _mitofinder_conf(a, conf, fasta):
    if not a.ref:
        sys.exit("[wgs_skim] MitoFinder needs a reference mitogenome GenBank: --ref related_species.gb")
    asm = [x.strip() for x in a.assembler.split(",") if x.strip()]
    bad = [x for x in asm if x not in MF_ASSEMBLERS]
    if bad:
        sys.exit(f"[wgs_skim] unknown --assembler {', '.join(bad)}: {' '.join(MF_ASSEMBLERS)}")
    conf.update(mitofinder_blast=not a.no_check, blast_db=a.db, mitofinder_ref=merged_ref(a.ref, a.outdir_path),
                mitofinder_trna=a.trna,
                mitofinder_assembly="assembly" if fasta else ",".join(asm))
    if not fasta:                                   # the assemblers that are workflow steps run first
        return ",".join(["mitofinder"] + [x for x in asm if x in ("spades", "rnaspades", "trinity")])


def _blast_opts(p):
    p.add_argument("--db", help="local BLAST database prefix (default: NCBI online BLAST)")
    p.add_argument("--genes", default="full", help="full = the whole sequence; or genes cut from a MitoM annotation "
                   "of each input, e.g. cox1,16s")
    p.add_argument("--topn", type=int, default=5)
    p.add_argument("--email", help="e-mail for NCBI Entrez (organism and taxonomy of the hits)")


def _blast_conf(a, conf, fasta):
    conf.update(blast_mt=a.genes, blast_nr="none", blast_db=a.db, email=a.email, blast_topn=a.topn)


TOOLS = {t.name: t for t in [
    Tool("fastp", "reads", "clean reads (adapters, quality) and their report -> fastp/<s>/, fastp/fastp_summary.tsv", "fastp"),
    Tool("meangs", "reads", "MEANGS mitochondrial seed -> meangs/<s>/", "meangs"),
    Tool("getorganelle", "reads", "GetOrganelle: mitogenome (reseeded rounds until circular; getorganelle/<s>/, mitogenome/<s>/) or rDNA + histones (rdna/, histone/)",
         "getorganelle", options=_getorganelle_opts, configure=_getorganelle_conf),
    Tool("spades", "reads", "SPAdes in any mode (genome rna meta plasmid ...) -> spades/<s>/", "spades",
         options=_spades_opts, configure=_spades_conf),
    Tool("trinity", "reads", "Trinity transcriptome -> trinity/<s>/<s>.trinity.fasta", "trinity"),
    Tool("genome-size", "reads", "jellyfish k-mers + GenomeScope2: genome size, heterozygosity -> genomescope/<s>/",
         "genome_size"),
    Tool("busco", "fasta", "BUSCO completeness of assemblies / transcriptomes / proteins -> busco/<s>/, busco/busco_summary.tsv",
         "busco", "scaffolds", _busco_opts, _busco_conf),
    Tool("mitofinder", "both", "MitoFinder: the mitogenome in assemblies (-i) or MitoFinder's own assembly of reads "
         "(-r) -> mitofinder_<source>/<s>/, mitofinder_<source>/mitofinder_summary.tsv", "mitofinder", "scaffolds", _mitofinder_opts, _mitofinder_conf),
    Tool("mitom", "fasta", "MitoM annotation of mitogenomes (genes, rRNA, tRNA; NCBI checks) -> mitom/<s>/",
         "annotate"),
    Tool("blast", "fasta", "BLAST sequences against nt, top hits with GenBank taxonomy -> blast/blast_topn.tsv",
         "blast", options=_blast_opts, configure=_blast_conf),
]}


def fasta_sheet(inputs, out, column, split_by_id=False, only=None):
    """sample sheet of FASTA inputs (no reads): <out>/pipeline_info/inputs/inputs.tsv (sample = file name, or the first word of
    the headers with split_by_id)"""
    if str(scripts_dir()) not in sys.path:
        sys.path.insert(0, str(scripts_dir()))
    from samples import fasta_samples
    meta = Path(out) / "pipeline_info" / "inputs"
    fastas = fasta_samples(inputs, meta / "split", split_by_id)
    if only:
        fastas = {k: v for k, v in fastas.items() if k in only.split(",")}
    meta.mkdir(parents=True, exist_ok=True)
    sheet = meta / "inputs.tsv"
    with open(sheet, "w") as fh:
        fh.write(f"sample\tfastq_1\tfastq_2\t{column}\n")
        fh.writelines(f"{k}\t\t\t{v}\n" for k, v in fastas.items())
    print(f"[wgs_skim] {len(fastas)} sample(s) from {', '.join(map(str, inputs))} -> {sheet}", flush=True)
    return sheet


def merged_ref(refs, out):
    """one GenBank file of every reference given (files, directories of .gb / .gbk / .gbf); a single file as it is"""
    files = [f for r in refs for f in (sorted(x for x in Path(r).iterdir() if x.suffix in (".gb", ".gbk", ".gbf",
                                                                                        ".genbank"))
                                       if Path(r).is_dir() else [Path(r)])]
    missing = [str(f) for f in files if not f.is_file()]
    if missing or not files:
        sys.exit(f"[wgs_skim] reference GenBank not found: {', '.join(missing) or ' '.join(refs)}")
    if len(files) == 1:
        return files[0].resolve()
    dst = Path(out) / "pipeline_info" / "inputs" / "mitofinder_refs.gb"
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("".join(f.read_text().rstrip() + "\n" for f in files))
    print(f"[wgs_skim] {len(files)} reference GenBank files merged -> {dst}", flush=True)
    return dst


def tool_config(a, out):
    """(workflow options, uses FASTA input) of `wgs_skim <tool>` arguments"""
    t = TOOLS[a.tool]
    use_fasta = t.inputs == "fasta" or (t.inputs == "both" and bool(getattr(a, "input", None)))
    if use_fasta and not a.input:
        sys.exit(f"[wgs_skim] give the FASTA files: wgs_skim {t.name} -i FILE|DIR [...]")
    if not use_fasta and not (a.reads_dir or a.samples):
        sys.exit(f"[wgs_skim] give the reads: wgs_skim {t.name} -r READS_DIR (or -s sheet.tsv)"
                 + (" or assemblies: -i FASTA" if t.inputs == "both" else ""))
    a.outdir_path = out
    conf = {"outdir": ".", "db_dir": HOME / "db", "taxon": a.taxon, "qc": True}
    if use_fasta:
        conf["samples"] = fasta_sheet(a.input, out, t.column, a.split_by_id, a.only)
        conf["qc"] = False                            # no reads
    else:
        conf["reads_dir"] = Path(a.reads_dir).resolve() if a.reads_dir else None
        conf["samples"] = Path(a.samples).resolve() if a.samples else None
        conf["qc"] = t.name == "fastp" or a.qc
    steps = t.configure(a, conf, use_fasta) if t.configure else None
    conf["steps"] = steps or t.step
    return conf, use_fasta

"""wgs_skim command line: a thin wrapper that builds the Snakemake (or Nextflow) command for you.

    wgs_skim run -r /path/to/rawdata                  # everything automatic (annelid preset)
    wgs_skim run -r rawdata -o results --taxon fish -j 64
    wgs_skim run -r rawdata -n                        # dry run: show what would run
    wgs_skim run -r rawdata --steps fastp,spades,busco   # only some steps (plus what they need)
    wgs_skim run -r rawdata --engine nextflow         # the same run with Nextflow (-resume, work/ in the results)
    wgs_skim run -r rawdata --steps spades-busco      # a combination: meangs-getorganelle spades-mitofinder ...
    wgs_skim spades -r clean_reads -o asm --mode rna  # one tool, every sample: fastp meangs getorganelle spades
    wgs_skim busco -i assemblies/ -o busco            #   trinity genome-size busco mitofinder mitom blast
    wgs_skim install                                  # create all software environments once
    wgs_skim samples -r rawdata                       # show detected sample names
    wgs_skim circularize -i mitofinder_mt.fasta --spades spades_dir   # circularize MitoFinder / SPAdes contigs
    wgs_skim compare -a getorganelle_dir -b mitofinder_mt.fasta --spades spades_dir   # are the assemblies the same?
    wgs_skim compare --results results/ -o compare/                              # every sample at once
    wgs_skim nuclear --results results/ -o nuclear/   # rDNA + histones: GetOrganelle, MitoM .gbf, vs SPAdes
    wgs_skim config > my.yaml                         # all options with comments

Anything after `--` goes to snakemake (or nextflow) unchanged, e.g. `wgs_skim run -r rawdata -- --set-threads spades=48`.
Personal defaults (local nt, e-mail, BUSCO lineage path ...) go in ~/.wgs_skim/config.yaml and apply to every run.
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from . import __version__
from .engine import ascii_paths, base_cmd, execute, nextflow, snakemake
from .paths import HOME, TAXA, default_conda_prefix, job_mem_gb, preset_file, presets, scripts_dir, workflow_root
from .tools import TOOLS, merged_ref, tool_config

# the modules behind the commands: paths (where things are), engine (Snakemake / Nextflow command), tools (the
# per-tool commands), api (the same from Python), lib (workflow/scripts as a library)


def cmd_run(args, passthrough):
    if not (args.reads_dir or args.samples or args.configfile or args.preset):
        sys.exit("[wgs_skim] give the reads directory: wgs_skim run -r /path/to/rawdata")
    reads = Path(args.reads_dir).resolve() if args.reads_dir else None
    out = Path(args.outdir or (reads.parent / "wgs_skim_results" if reads else "wgs_skim_results")).resolve()
    conf = {"reads_dir": reads, "samples": Path(args.samples).resolve() if args.samples else None,
            "outdir": ".", "db_dir": HOME / "db", "taxon": args.taxon, "blast_db": args.blast_db, "email": args.email,
            "busco_lineage": args.busco_lineage, "steps": args.steps,
            "blast_mt": args.blast_mt, "blast_nr": args.blast_nr, "spades_modes": args.spades_modes,
            "mitofinder_ref": merged_ref(args.mitofinder_ref, out) if args.mitofinder_ref else None}
    if args.no_qc:                                    # the reads are already clean (e.g. fastp/ of an earlier run)
        conf["qc"] = False
    if args.mt_rounds is not None:
        conf["mt_iter_rounds"] = args.mt_rounds
    files = ([preset_file(args.preset)] if args.preset else []) + ([args.configfile] if args.configfile else [])
    run = nextflow if args.engine == "nextflow" else snakemake
    return run(conf, out, args.cores, job_mem_gb(args.mem_gb), args.conda_prefix, files, args.config,
               args.dry_run, passthrough)


def cmd_gui(args, passthrough):
    from .gui import serve
    return serve(args.port, not args.no_browser)


def cmd_presets(args, passthrough):
    """the presets: name, file, first comment line (what it runs)"""
    for name, f in presets().items():
        first = next((x.lstrip("# ").strip() for x in f.read_text().splitlines() if x.startswith("#")), "")
        print(f"{name:16s} {first}\n{'':16s} {f}")
    print(f"\nyour own: {HOME / 'presets'}/<name>.yaml (any workflow options, see `wgs_skim config`); "
          "use: wgs_skim run -r rawdata/ --preset <name>")
    return 0


def cmd_tool(args, passthrough):
    out = Path(args.outdir or f"{args.tool}_results").resolve()
    ascii_paths(out, {"reads_dir": args.reads_dir, "samples": args.samples} if hasattr(args, "reads_dir") else {})
    out.mkdir(parents=True, exist_ok=True)
    conf, _ = tool_config(args, out)
    return snakemake(conf, out, args.cores, job_mem_gb(args.mem_gb), args.conda_prefix, args.configfile,
                     args.config, args.dry_run, passthrough)


def cmd_install(args, passthrough):
    """Create every rule's conda environment (one-time, shared by all projects)."""
    root = workflow_root()
    tmp = HOME / "install"
    tmp.mkdir(parents=True, exist_ok=True)
    (tmp / "reads").mkdir(exist_ok=True)
    for r in ("x_R1.fq.gz", "x_R2.fq.gz"):     # a placeholder sample so that every rule is in the DAG
        (tmp / "reads" / r).touch()
    conf = {"reads_dir": tmp / "reads", "outdir": ".", "db_dir": HOME / "db"}
    base = base_cmd(conf, args.conda_prefix, args.configfile, args.config)
    print(f"[wgs_skim] creating environments in {args.conda_prefix} (first time: 20-60 min)", flush=True)
    rc = execute(base + ["-d", str(tmp), "--cores", "1", "--conda-create-envs-only"] + passthrough, cwd=str(root))
    if rc == 0:                                       # databases every run needs: MitoM, MITOS2 refs, GetOrganelle
        print(f"[wgs_skim] fetching MitoM and the GetOrganelle / MITOS2 databases into {HOME / 'db'}", flush=True)
        rc = execute(base + ["-d", str(tmp), "--cores", "2", "mitom_db", "mitom_nuc_db", "getorganelle_db"]
                     + passthrough, cwd=str(root))
    if rc == 0:
        print("[wgs_skim] all environments and databases ready")
    return rc


def cmd_samples(args, passthrough):
    if not (args.reads_dir or args.samples):
        sys.exit("[wgs_skim] give the reads directory or a sample sheet: wgs_skim samples -r DIR | -s sheet.tsv")
    from .lib.samples import find_pairs, read_sheet
    pairs = read_sheet(args.samples) if args.samples else find_pairs(args.reads_dir)
    for s, (r1, r2) in pairs.items():
        print(s, r1, r2, sep="\t")
    print(f"[wgs_skim] {len(pairs)} samples", file=sys.stderr)
    return 0


def cmd_circularize(args, passthrough):
    """Circularize mitochondrial contigs (MitoFinder / SPAdes) without reads; see workflow/scripts/circularize.py."""
    cmd = [sys.executable, str(scripts_dir() / "circularize.py"), "-i", args.input,
           "-o", args.outdir, "--prefix", args.prefix] + (["--spades", args.spades] if args.spades else [])
    return subprocess.run(cmd).returncode


def compare_tools(prefix):
    """bin directory with mafft + minimap2: <conda_prefix>/wgs_skim_compare, created on first use (69 MB:
    minimap2 without its k8 / sysroot dependencies, which only paftools.js needs - with them 564 MB)."""
    env = prefix / "wgs_skim_compare"
    if not ((env / "bin" / "mafft").exists() and (env / "bin" / "minimap2").exists()):
        conda = shutil.which("mamba") or shutil.which("conda") or next(
            (str(p) for p in (Path.home() / "miniforge3" / "bin" / x for x in ("mamba", "conda")) if p.exists()),
            None)
        if not conda:
            sys.exit("[wgs_skim] compare needs mafft and minimap2: no conda / mamba found to install them")
        print(f"[wgs_skim] installing mafft + minimap2 into {env} (once)", file=sys.stderr, flush=True)
        ch = ["-c", "conda-forge", "-c", "bioconda", "--override-channels"]
        subprocess.run([conda, "create", "-y", "-q", "-p", str(env), *ch, "mafft=7.526", "libzlib"], check=True)
        subprocess.run([conda, "install", "-y", "-q", "-p", str(env), *ch, "--no-deps", "minimap2=2.28"], check=True)
    return env / "bin"


def cmd_compare(args, passthrough):
    """GetOrganelle vs other mitogenome assemblies: one sample directly; a batch (--results / --sheet) through
    workflow/compare_batch.smk (Snakemake: scheduled by cores and memory, resumable). See scripts/compare_mt.py."""
    tools = compare_tools(Path(args.conda_prefix))
    if args.results or args.sheet:
        out = Path(args.outdir).resolve()
        conf = {"outdir": out, "tools": tools, "plot": "false" if args.no_plot else "true", "samples": args.samples,
                "results": Path(args.results).resolve() if args.results else None,
                "sheet": Path(args.sheet).resolve() if args.sheet else None}
        return snakemake(conf, out, args.cores, job_mem_gb(args.mem_gb), args.conda_prefix,
                         passthrough=passthrough, snakefile="compare_batch.smk")
    env = dict(os.environ)
    env["PATH"] = f"{tools}{os.pathsep}{env['PATH']}"
    cmd = [sys.executable, str(scripts_dir() / "compare_mt.py"),
           "-o", args.outdir, "--prefix", args.prefix]
    cmd += [x for b in args.other or [] for x in ("-b", b)]
    cmd += ["--no-plot"] if args.no_plot else []
    cmd += [x for k in ("ref", "spades", "gb", "gff", "gfa") if getattr(args, k) for x in (f"--{k}", getattr(args, k))]
    return subprocess.run(cmd, env=env).returncode


def cmd_nuclear(args, passthrough):
    """rDNA + histones of every sample of an existing results directory: workflow/nuclear_batch.smk (Snakemake:
    jobs scheduled by cores and memory, one log each, resumable); --status reads the output directory."""
    nm = str(scripts_dir() / "nuclear_markers.py")
    if args.status:
        return subprocess.run([sys.executable, nm, "status", "-o", args.outdir, "--markers", args.markers]).returncode
    if not args.results:
        sys.exit("[wgs_skim] give the results directory: wgs_skim nuclear --results DIR -o OUT")
    mitom = Path(args.mitom) if args.mitom else HOME / "db" / "mitom_nuc" / "Mitom"
    if not (mitom / "mito_nuc.py").exists():
        sys.exit(f"[wgs_skim] {mitom}: no MitoM with `mitom nuc` (--mitom, or `wgs_skim install`)")
    out = Path(args.outdir).resolve()
    conf = {"results": Path(args.results).resolve(), "outdir": out, "markers": args.markers, "max_gb": args.max_gb,
            "mitom": mitom.resolve(), "samples": args.samples}
    return snakemake(conf, out, args.cores, job_mem_gb(args.mem_gb), args.conda_prefix, dry_run=args.dry_run,
                     passthrough=passthrough, snakefile="nuclear_batch.smk")


def cmd_build_blastdb(args, passthrough):
    """Small marker BLAST database (one taxon, short records) from a local nt; see build_marker_db.py."""
    cmd = [sys.executable, str(scripts_dir() / "build_marker_db.py"), "--nt", args.nt,
           "--out", args.out, "--taxid", str(args.taxid), "--max-length", str(args.max_length)]
    return subprocess.run(cmd + passthrough).returncode


def cmd_config(args, passthrough):
    sys.stdout.write((workflow_root() / "config" / "config.yaml").read_text())
    return 0


def parser():
    """the wgs_skim argument parser (every command)"""
    p = argparse.ArgumentParser(prog="wgs_skim", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-v", "--version", action="version", version=f"wgs_skim {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp):
        sp.add_argument("--conda-prefix", default=default_conda_prefix(),
                        help="where software environments are kept, shared by all runs (default: %(default)s)")
        sp.add_argument("--configfile", help="YAML file with any workflow options (see `wgs_skim config`)")
        sp.add_argument("--config", nargs="+", metavar="KEY=VALUE", help="set workflow options, as in snakemake")

    r = sub.add_parser("run", help="run the workflow", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    r.add_argument("-r", "--reads_dir", help="directory of paired FASTQ files (sample names detected automatically)")
    r.add_argument("-s", "--samples", help="sample sheet instead of -r: TSV/CSV with sample, fastq_1, fastq_2 "
                   "[, scaffolds: an existing SPAdes assembly to use instead of assembling]")
    r.add_argument("-o", "--outdir", help="results directory (default: <reads_dir>/../wgs_skim_results)")
    r.add_argument("-t", "--taxon", choices=TAXA, help="taxon preset (default annelid)")
    r.add_argument("-j", "--cores", type=int, default=os.cpu_count(), help="CPU cores to use")
    r.add_argument("--mem-gb", type=int, default=0,
                   help="memory the jobs may use together (default: 90%% of the machine); every rule declares "
                        "its need, SPAdes from the data size")
    r.add_argument("-n", "--dry-run", action="store_true", help="only show what would be run")
    r.add_argument("--engine", choices=["snakemake", "nextflow"], default="snakemake",
                   help="workflow manager: Snakemake (workflow/Snakefile) or Nextflow (nextflow/main.nf, run with "
                        "-resume in the results directory)")
    r.add_argument("--steps", metavar="STEP[,STEP]",
                   help="run only these steps (and what they need): fastp, meangs, mito (MEANGS + GetOrganelle), "
                        "getorganelle (without MEANGS), nr (rDNA + histones), spades, rnaspades, busco, genome_size, "
                        "annotate (MitoM), topology, mitofinder, blast, nuclear (rDNA + histones in the SPAdes "
                        "assembly, copy numbers, vs GetOrganelle), trinity, rna (rnaSPAdes + Trinity + BUSCO), "
                        "summary; or combinations: meangs-getorganelle, spades-mitofinder, spades-busco, "
                        "fastp-spades-busco, rnaspades-busco, trinity-busco (default: everything)")
    r.add_argument("--spades-modes", metavar="M[,M]",
                   help="SPAdes modes: genome (default) and/or isolate, sc, meta, rna, rnaviral, plasmid, "
                        "metaplasmid, metaviral, bio, sewage, corona (one assembly + BUSCO each)")
    r.add_argument("--mitofinder-ref", nargs="+", metavar="GB|DIR", help="(several files / a directory: merged) "
                   + "MitoFinder reference GenBank (combination spades-mitofinder; default: "
                   "the sample's own MitoM annotation when the mitogenome step runs too)")
    r.add_argument("--blast-db", help="local nt database prefix (default: NCBI online BLAST)")
    r.add_argument("--blast-mt", metavar="Q[,Q]",
                   help="mitogenome BLAST queries: full and/or genes, e.g. cox1,16s,cytb (default cox1)")
    r.add_argument("--blast-nr", metavar="Q[,Q]",
                   help="rDNA BLAST queries: full and/or 18S, ITS, ITS1, 5.8S, ITS2, 28S (default 18S,ITS,28S)")
    r.add_argument("--email", help="e-mail for NCBI Entrez (BLAST metadata)")
    r.add_argument("--preset", metavar="NAME|FILE", help="a set of steps and options: full (everything incl. "
                   "SPAdes-MitoFinder), mito, nuclear, genome, transcriptome, or your own ~/.wgs_skim/presets/NAME.yaml "
                   "(`wgs_skim presets` lists them); --steps / --config on the command line still win")
    r.add_argument("--mt-rounds", type=int, metavar="N", help="GetOrganelle rounds reseeded with the previous "
                   "result while the mitogenome is not circular (default 3; 0 = one run)")
    r.add_argument("--no-qc", action="store_true", help="the reads are already clean (e.g. the fastp/ directory of an "
                   "earlier run): used as they are, no fastp, no copies")
    r.add_argument("--busco-lineage", help="BUSCO lineage name or downloaded lineage directory")
    common(r)
    r.set_defaults(func=cmd_run)

    for t, tool in TOOLS.items():
        kind, helptext = tool.inputs, tool.help
        tp = sub.add_parser(t, help=helptext, formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                            description=f"{helptext}. Every sample in one batch; resumable (run the same command "
                                        "again). Options after -- go to snakemake.")
        if kind in ("reads", "both"):
            tp.add_argument("-r", "--reads_dir", help="directory of paired FASTQ files (sample names detected)")
            tp.add_argument("-s", "--samples", help="sample sheet instead of -r: sample, fastq_1, fastq_2")
            if t != "fastp":
                tp.add_argument("--qc", action="store_true", help="clean the reads with fastp first "
                                "(default: the reads are used as they are)")
        if kind in ("fasta", "both"):
            tp.add_argument("-i", "--input", nargs="+", default=[], metavar="FASTA|DIR",
                            help="FASTA files or directories of them (sample = file name)")
            tp.add_argument("--split-by-id", action="store_true", help="the input is one FASTA of many samples: "
                            "split it by the first word of the headers (>B2 ..., >B17 ...)")
            tp.add_argument("--only", help="only these samples (comma-separated)")
        tp.add_argument("-o", "--outdir", help=f"results directory (default: ./{t}_results)")
        tp.add_argument("-t", "--taxon", choices=TAXA, help="taxon preset: genetic code, MitoM clade, BUSCO "
                        "lineage, MEANGS class (default annelid)")
        tp.add_argument("-j", "--cores", type=int, default=os.cpu_count(), help="CPU cores for all jobs together")
        tp.add_argument("--mem-gb", type=int, default=0, help="memory for all jobs together (default 90%% of the "
                        "machine)")
        tp.add_argument("-n", "--dry-run", action="store_true", help="only show what would run")
        if tool.options:
            tool.options(tp)
        common(tp)
        tp.set_defaults(func=cmd_tool, tool=t)

    i = sub.add_parser("install", help="create all software environments (once)")
    common(i)
    i.set_defaults(func=cmd_install)

    s = sub.add_parser("samples", help="list the samples detected in a reads directory")
    s.add_argument("-r", "--reads_dir")
    s.add_argument("-s", "--samples")
    s.set_defaults(func=cmd_samples)

    z = sub.add_parser("circularize", help="circularize mitochondrial contigs (e.g. MitoFinder output) without reads")
    z.add_argument("-i", "--input", required=True, help="mt contigs FASTA, e.g. MitoFinder <sample>_all_mt.fasta")
    z.add_argument("-o", "--outdir", default="circularize")
    z.add_argument("--spades", help="SPAdes output directory: use its assembly graph (strongest check, fills gaps)")
    z.add_argument("--prefix", default="sample")
    z.set_defaults(func=cmd_circularize)

    m = sub.add_parser("compare", help="compare GetOrganelle with SPAdes-MitoFinder (or other) mitogenomes")
    m.add_argument("-a", "--ref", help="GetOrganelle output directory or FASTA (reference)")
    m.add_argument("-b", "--other", action="append", help="[label=]FASTA to compare (repeatable)")
    m.add_argument("--spades", help="SPAdes directory of the first -b assembly (graph-based circularization)")
    m.add_argument("--gb", help="GenBank annotation of the reference (gene of each difference)")
    m.add_argument("--gff", help="GFF annotation on the reference sequence, instead of --gb (as QUAST -g)")
    m.add_argument("--gfa", help="reference assembly graph when -a is a FASTA")
    m.add_argument("--results", help="batch: results directory (all samples with GetOrganelle + MitoFinder)")
    m.add_argument("--sheet", help="batch: sample sheet TSV with sample, ref, other, [spades, gb, gfa]")
    m.add_argument("--samples", help="batch: only these samples (comma-separated)")
    m.add_argument("-j", "--cores", type=int, default=os.cpu_count(), help="batch: CPU cores for all jobs together")
    m.add_argument("--mem-gb", type=int, default=0, help="batch: memory for all jobs (default 90%% of the machine; "
                   "~8 GB per sample with a SPAdes graph)")
    m.add_argument("--no-plot", action="store_true", help="no figures (default: one SVG per comparison)")
    m.add_argument("-o", "--outdir", default="compare")
    m.add_argument("--prefix", default="sample")
    m.add_argument("--conda-prefix", default=default_conda_prefix(),
                   help="where software environments are kept (mafft + minimap2 go to wgs_skim_compare there)")
    m.set_defaults(func=cmd_compare)

    n = sub.add_parser("nuclear", help="rDNA + histones of a results directory: GetOrganelle, MitoM .gbf, vs SPAdes")
    n.add_argument("--results", default="", help="results directory (wgs_skim or wgs_skim_pipeline layout)")
    n.add_argument("-o", "--outdir", default="nuclear_markers")
    n.add_argument("--samples", help="only these samples (comma-separated)")
    n.add_argument("--markers", default="rdna,histone")
    n.add_argument("-j", "--cores", type=int, default=os.cpu_count(), help="CPU cores for all jobs together")
    n.add_argument("--mem-gb", type=int, default=0, help="memory for all jobs together (default 90%% of the machine)")
    n.add_argument("-n", "--dry-run", action="store_true", help="only show what would run")
    n.add_argument("--max-gb", default="5,10", help="Gb of reads per GetOrganelle round (round 2, seeded with "
                   "round 1, only when round 1 does not close the unit)")
    n.add_argument("--mitom", help="MitoM checkout (default: the workflow's, <WGS_SKIM_HOME>/db/mitom_nuc/Mitom)")
    n.add_argument("--conda-prefix", default=default_conda_prefix(), help="where the software environments are")
    n.add_argument("--status", action="store_true", help="only show what the batch in -o has done so far")
    n.set_defaults(func=cmd_nuclear)

    b = sub.add_parser("build-blastdb", help="build a fast marker BLAST database (Metazoa, <= 50 kb) from local nt")
    b.add_argument("--nt", required=True, help="local nt prefix, e.g. /db/blast/nt")
    b.add_argument("--out", required=True, help="output prefix, e.g. /db/blast_marker/metazoa_marker")
    b.add_argument("--taxid", type=int, default=33208, help="taxon to keep (default 33208 = Metazoa)")
    b.add_argument("--max-length", type=int, default=50000, help="longest record kept (default 50000 bp)")
    b.set_defaults(func=cmd_build_blastdb)

    pr = sub.add_parser("presets", help="list the presets (wgs_skim run --preset NAME)")
    pr.set_defaults(func=cmd_presets)

    g = sub.add_parser("gui", help="tick functions and options in a web page, save them as a preset")
    g.add_argument("--port", type=int, default=8765)
    g.add_argument("--no-browser", action="store_true", help="only print the address")
    g.set_defaults(func=cmd_gui)

    c = sub.add_parser("config", help="print the default configuration with comments")
    c.set_defaults(func=cmd_config)

    return p


def dispatch(argv):
    """run one wgs_skim command line (list, without the program name); returns the exit code (wgs_skim.api)"""
    argv = [str(x) for x in argv]
    passthrough = argv[argv.index("--") + 1:] if "--" in argv else []
    argv = argv[:argv.index("--")] if "--" in argv else argv
    args = parser().parse_args(argv)
    return args.func(args, passthrough)


def main(argv=None):
    sys.exit(dispatch(sys.argv[1:] if argv is None else argv))


if __name__ == "__main__":
    main()

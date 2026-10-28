"""wgs_skim command line: a thin wrapper that builds the Snakemake command for you.

    wgs_skim run -r /path/to/rawdata                  # everything automatic (annelid preset)
    wgs_skim run -r rawdata -o results --taxon fish -j 64
    wgs_skim run -r rawdata -n                        # dry run: show what would run
    wgs_skim run -r rawdata --steps fastp,spades,busco   # only some steps (plus what they need)
    wgs_skim install                                  # create all software environments once
    wgs_skim samples -r rawdata                       # show detected sample names
    wgs_skim circularize -i mitofinder_mt.fasta --spades spades_dir   # circularize MitoFinder / SPAdes contigs
    wgs_skim compare -a getorganelle_dir -b mitofinder_mt.fasta --spades spades_dir   # are the assemblies the same?
    wgs_skim config > my.yaml                         # all options with comments

Anything after `--` goes to snakemake unchanged, e.g. `wgs_skim run -r rawdata -- --set-threads spades=48`.
Personal defaults (local nt, e-mail, BUSCO lineage path ...) go in ~/.wgs_skim/config.yaml and apply to every run.
"""
import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from . import __version__

HOME = Path(os.environ.get("WGS_SKIM_HOME", Path.home() / ".wgs_skim"))
USER_CONFIG = HOME / "config.yaml"


def default_conda_prefix():
    """WGS_SKIM_CONDA_PREFIX, else `conda_prefix:` in ~/.wgs_skim/config.yaml, else ~/.wgs_skim/conda."""
    if os.environ.get("WGS_SKIM_CONDA_PREFIX"):
        return os.environ["WGS_SKIM_CONDA_PREFIX"]
    if USER_CONFIG.exists():
        for line in USER_CONFIG.read_text().splitlines():
            if line.startswith("conda_prefix:"):
                return line.split(":", 1)[1].split("#")[0].strip()
    return HOME / "conda"
TAXA = ["annelid", "mollusc", "fish", "crustacean"]


def workflow_root():
    """Directory holding workflow/ config/ resources/: the installed package data, or the git checkout."""
    here = Path(__file__).resolve().parent
    for root in (here, here.parent.parent):
        if (root / "workflow" / "Snakefile").is_file():
            return root
    sys.exit("[wgs_skim] cannot find workflow/Snakefile (broken installation?)")


def snakemake_exe():
    """Snakemake >= 8 from the same environment as wgs_skim, else from PATH."""
    local = Path(sys.executable).parent / "snakemake"
    exe = str(local) if local.exists() else shutil.which("snakemake")
    if not exe:
        sys.exit("[wgs_skim] snakemake not found: conda install -c conda-forge -c bioconda 'snakemake>=8'")
    return exe


def base_cmd(args, extra_config):
    root = workflow_root()
    cmd = [snakemake_exe(), "-s", str(root / "workflow" / "Snakefile"),
           "--sdm", "conda", "--conda-prefix", str(args.conda_prefix)]
    if USER_CONFIG.is_file():
        cmd += ["--configfile", str(USER_CONFIG)]
    if getattr(args, "configfile", None):
        cmd += ["--configfile", str(Path(args.configfile).resolve())]
    cfg = [f"{k}={v}" for k, v in extra_config.items() if v not in (None, "")] + (getattr(args, "config", None) or [])
    if cfg:
        cmd += ["--config", *cfg]
    return cmd


def cmd_run(args, passthrough):
    if not (args.reads_dir or args.samples or args.configfile):
        sys.exit("[wgs_skim] give the reads directory: wgs_skim run -r /path/to/rawdata")
    reads = Path(args.reads_dir).resolve() if args.reads_dir else None
    out = Path(args.outdir or (reads.parent / "wgs_skim_results" if reads else "wgs_skim_results")).resolve()
    out.mkdir(parents=True, exist_ok=True)
    conf = {"reads_dir": reads, "samples": Path(args.samples).resolve() if args.samples else None,
            "outdir": ".", "db_dir": HOME / "db", "taxon": args.taxon, "blast_db": args.blast_db, "email": args.email,
            "busco_lineage": args.busco_lineage, "steps": args.steps,
            "blast_mt": args.blast_mt, "blast_nr": args.blast_nr}
    cmd = base_cmd(args, conf) + ["-d", str(out), "--cores", str(args.cores)]
    if args.dry_run:
        cmd.append("-n")
    cmd += passthrough
    print(f"[wgs_skim] results: {out}\n[wgs_skim] {shlex.join(cmd)}", flush=True)
    return subprocess.run(cmd).returncode


def cmd_install(args, passthrough):
    """Create every rule's conda environment (one-time, shared by all projects)."""
    root = workflow_root()
    tmp = HOME / "install"
    tmp.mkdir(parents=True, exist_ok=True)
    (tmp / "reads").mkdir(exist_ok=True)
    for r in ("x_R1.fq.gz", "x_R2.fq.gz"):     # a placeholder sample so that every rule is in the DAG
        (tmp / "reads" / r).touch()
    cmd = base_cmd(args, {"reads_dir": tmp / "reads", "outdir": ".", "db_dir": HOME / "db"}) + \
        ["-d", str(tmp), "--cores", "1", "--conda-create-envs-only"] + passthrough
    print(f"[wgs_skim] creating environments in {args.conda_prefix} (first time: 20-60 min)\n"
          f"[wgs_skim] {shlex.join(cmd)}", flush=True)
    rc = subprocess.run(cmd, cwd=str(root)).returncode
    if rc == 0:
        print("[wgs_skim] all environments ready")
    return rc


def cmd_samples(args, passthrough):
    sys.path.insert(0, str(workflow_root() / "workflow" / "scripts"))
    from samples import find_pairs, read_sheet
    pairs = read_sheet(args.samples) if args.samples else find_pairs(args.reads_dir)
    for s, (r1, r2) in pairs.items():
        print(s, r1, r2, sep="\t")
    print(f"[wgs_skim] {len(pairs)} samples", file=sys.stderr)
    return 0


def cmd_circularize(args, passthrough):
    """Circularize mitochondrial contigs (MitoFinder / SPAdes) without reads; see workflow/scripts/circularize.py."""
    cmd = [sys.executable, str(workflow_root() / "workflow" / "scripts" / "circularize.py"), "-i", args.input,
           "-o", args.outdir, "--prefix", args.prefix] + (["--spades", args.spades] if args.spades else [])
    return subprocess.run(cmd).returncode


def cmd_compare(args, passthrough):
    """GetOrganelle vs other mitogenome assemblies of one sample; see workflow/scripts/compare_mt.py."""
    env = dict(os.environ)
    if not shutil.which("mafft"):               # MAFFT from the workflow's own environments
        found = sorted(Path(args.conda_prefix).glob("*/bin/mafft"))
        if not found:
            sys.exit("[wgs_skim] mafft not found: put it on PATH or run `wgs_skim install` first")
        env["PATH"] = f"{found[0].parent}{os.pathsep}{env['PATH']}"
    cmd = [sys.executable, str(workflow_root() / "workflow" / "scripts" / "compare_mt.py"), "-a", args.ref,
           "-o", args.outdir, "--prefix", args.prefix]
    cmd += [x for b in args.other for x in ("-b", b)]
    cmd += [x for k in ("spades", "gb", "gfa") if getattr(args, k) for x in (f"--{k}", getattr(args, k))]
    return subprocess.run(cmd, env=env).returncode


def cmd_build_blastdb(args, passthrough):
    """Small marker BLAST database (one taxon, short records) from a local nt; see build_marker_db.py."""
    cmd = [sys.executable, str(workflow_root() / "workflow" / "scripts" / "build_marker_db.py"), "--nt", args.nt,
           "--out", args.out, "--taxid", str(args.taxid), "--max-length", str(args.max_length)]
    return subprocess.run(cmd + passthrough).returncode


def cmd_config(args, passthrough):
    sys.stdout.write((workflow_root() / "config" / "config.yaml").read_text())
    return 0


def main(argv=None):
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
    r.add_argument("-s", "--samples", help="sample sheet instead of -r: TSV/CSV with sample, fastq_1, fastq_2")
    r.add_argument("-o", "--outdir", help="results directory (default: <reads_dir>/../wgs_skim_results)")
    r.add_argument("-t", "--taxon", choices=TAXA, help="taxon preset (default annelid)")
    r.add_argument("-j", "--cores", type=int, default=os.cpu_count(), help="CPU cores to use")
    r.add_argument("-n", "--dry-run", action="store_true", help="only show what would be run")
    r.add_argument("--steps", metavar="STEP[,STEP]",
                   help="run only these steps (and what they need): fastp, meangs, mito (MEANGS + GetOrganelle), "
                        "getorganelle (without MEANGS), nr, spades, busco, genome_size, annotate, topology, blast, "
                        "summary (default: everything)")
    r.add_argument("--blast-db", help="local nt database prefix (default: NCBI online BLAST)")
    r.add_argument("--blast-mt", metavar="Q[,Q]",
                   help="mitogenome BLAST queries: full and/or genes, e.g. cox1,16s,cytb (default cox1)")
    r.add_argument("--blast-nr", metavar="Q[,Q]",
                   help="rDNA BLAST queries: full and/or 18S, ITS, ITS1, 5.8S, ITS2, 28S (default 18S,ITS,28S)")
    r.add_argument("--email", help="e-mail for NCBI Entrez (BLAST metadata)")
    r.add_argument("--busco-lineage", help="BUSCO lineage name or downloaded lineage directory")
    common(r)
    r.set_defaults(func=cmd_run)

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
    m.add_argument("-a", "--ref", required=True, help="GetOrganelle output directory or FASTA (reference)")
    m.add_argument("-b", "--other", action="append", required=True, help="[label=]FASTA to compare (repeatable)")
    m.add_argument("--spades", help="SPAdes directory of the other assembly (graph-based circularization)")
    m.add_argument("--gb", help="GenBank annotation of the reference (gene of each difference)")
    m.add_argument("--gfa", help="reference assembly graph when -a is a FASTA")
    m.add_argument("-o", "--outdir", default="compare")
    m.add_argument("--prefix", default="sample")
    m.add_argument("--conda-prefix", default=default_conda_prefix(),
                   help="where `wgs_skim install` put the environments (MAFFT is taken from there)")
    m.set_defaults(func=cmd_compare)

    b = sub.add_parser("build-blastdb", help="build a fast marker BLAST database (Metazoa, <= 50 kb) from local nt")
    b.add_argument("--nt", required=True, help="local nt prefix, e.g. /db/blast/nt")
    b.add_argument("--out", required=True, help="output prefix, e.g. /db/blast_marker/metazoa_marker")
    b.add_argument("--taxid", type=int, default=33208, help="taxon to keep (default 33208 = Metazoa)")
    b.add_argument("--max-length", type=int, default=50000, help="longest record kept (default 50000 bp)")
    b.set_defaults(func=cmd_build_blastdb)

    c = sub.add_parser("config", help="print the default configuration with comments")
    c.set_defaults(func=cmd_config)

    argv = list(sys.argv[1:] if argv is None else argv)
    passthrough = argv[argv.index("--") + 1:] if "--" in argv else []
    argv = argv[:argv.index("--")] if "--" in argv else argv
    args = p.parse_args(argv)
    sys.exit(args.func(args, passthrough))


if __name__ == "__main__":
    main()

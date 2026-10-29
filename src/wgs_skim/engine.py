"""Workflow engines: the Snakemake (or Nextflow) command for a run, built from a dict of workflow options.

    from wgs_skim.engine import snakemake, nextflow
    snakemake({"reads_dir": "/data/raw", "steps": "spades-busco"}, out="results", cores=32, mem_gb=200)
Options set in ~/.wgs_skim/config.yaml or a --configfile are not overridden by the built-in defaults (db_dir), and
KEY=VALUE pairs given by the user (`config=[...]`) win over everything, as with `snakemake --config`."""
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from .paths import USER_CONFIG, workflow_root


def snakemake_exe():
    """Snakemake >= 8 from the same environment as wgs_skim, else from PATH."""
    local = Path(sys.executable).parent / "snakemake"
    exe = str(local) if local.exists() else shutil.which("snakemake")
    if not exe:
        sys.exit("[wgs_skim] snakemake not found: conda install -c conda-forge -c bioconda 'snakemake>=8'")
    return exe


def user_files(configfile=None):
    """~/.wgs_skim/config.yaml and the --configfile, those that exist"""
    return [Path(f) for f in [USER_CONFIG] + ([configfile] if configfile else []) if Path(f).is_file()]


def user_config(configfile=None):
    """the options of ~/.wgs_skim/config.yaml and the --configfile (the latter wins)"""
    import yaml
    params = {}
    for f in user_files(configfile):
        params.update(yaml.safe_load(f.read_text()) or {})
    return params


def config_args(conf, configfile=None, config=()):
    """--config KEY=VALUE list: `conf` minus empty values (and db_dir when the user set it), then `config`"""
    given = set(user_config(configfile))
    cfg = [f"{k}={v}" for k, v in conf.items() if v not in (None, "") and not (k == "db_dir" and k in given)]
    return cfg + list(config or [])


def base_cmd(conf, conda_prefix, configfile=None, config=(), snakefile="Snakefile"):
    """snakemake -s workflow/<snakefile> --sdm conda ... --configfile ... --config ..."""
    cmd = [snakemake_exe(), "-s", str(workflow_root() / "workflow" / snakefile),
           "--sdm", "conda", "--conda-prefix", str(conda_prefix)]
    for f in user_files(configfile) if snakefile == "Snakefile" else []:
        cmd += ["--configfile", str(f.resolve())]
    cfg = config_args(conf, configfile, config) if snakefile == "Snakefile" else \
        [f"{k}={v}" for k, v in conf.items() if v not in (None, "")]
    return cmd + (["--config", *cfg] if cfg else [])


def resources(cores, mem_gb):
    """--cores and, when known, the memory every job's mem_mb is scheduled against"""
    return ["--cores", str(cores)] + (["--resources", f"mem_mb={int(mem_gb * 1024)}"] if mem_gb else [])


def execute(cmd, cwd=None, env=None, show=None):
    print((f"[wgs_skim] results: {show}\n" if show else "") + f"[wgs_skim] {shlex.join(map(str, cmd))}", flush=True)
    return subprocess.run([str(c) for c in cmd], cwd=cwd, env=env).returncode


def snakemake(conf, out, cores, mem_gb=0, conda_prefix=None, configfile=None, config=(), dry_run=False,
              passthrough=(), snakefile="Snakefile"):
    """snakemake on workflow/<snakefile> with results in `out`, scheduled by cores and memory (GB); exit code"""
    from .paths import default_conda_prefix
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    cmd = base_cmd(conf, conda_prefix or default_conda_prefix(), configfile, config, snakefile)
    cmd += ["-d", str(out)] + resources(cores, mem_gb) + (["-n"] if dry_run else []) + list(passthrough)
    return execute(cmd, show=out)


# steps / SPAdes modes the Nextflow version does not have (it refuses instead of running half a workflow)
NEXTFLOW_MISSING = {"rna", "rnaspades", "trinity", "transcriptome", "mitofinder", "spades-mitofinder",
                    "rnaspades-mitofinder", "trinity-mitofinder", "rnaspades-busco", "trinity-busco"}


def nextflow(conf, out, cores, mem_gb=0, conda_prefix=None, configfile=None, config=(), dry_run=False,
             passthrough=()):
    """The same run with Nextflow (nextflow/main.nf): the same options as params, always -resume, work/ in the
    results directory. ~/.wgs_skim/config.yaml and --configfile are merged into one -params-file (Nextflow takes
    one); command-line options win over both, as with Snakemake."""
    import yaml
    from .paths import default_conda_prefix
    exe = shutil.which("nextflow")
    if not exe:
        sys.exit("[wgs_skim] nextflow not found: conda install -c bioconda nextflow (or use --engine snakemake)")
    other = [m for m in str(conf.get("spades_modes") or "genome").split(",") if m.strip() not in ("", "genome")]
    if other or set(str(conf.get("steps") or "").split(",")) & NEXTFLOW_MISSING:
        sys.exit("[wgs_skim] other SPAdes modes, Trinity and MitoFinder run with the Snakemake engine only "
                 "(drop --engine nextflow)")
    conda_prefix = conda_prefix or default_conda_prefix()
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    params = user_config(configfile)
    params_file = out / "nextflow_params.yaml"
    params_file.write_text(yaml.safe_dump(params, sort_keys=False))
    cmd = [exe, "run", str(workflow_root() / "nextflow" / "main.nf"), "-profile", "conda", "-resume",
           "-work-dir", str(out / "work"), "-params-file", str(params_file)]
    opts = {k: v for k, v in dict(conf, outdir=out).items()
            if v not in (None, "") and not (k == "db_dir" and k in params)}
    opts |= {"max_cpus": cores, "python": sys.executable, "conda_prefix": conda_prefix}
    if mem_gb:
        opts["max_memory"] = f"{mem_gb} GB"
    for kv in config or []:
        k, _, v = kv.partition("=")
        opts[k] = v
    cmd += [x for k, v in opts.items() if v not in (None, "") for x in (f"--{k}", str(v))]
    cmd += (["-preview"] if dry_run else []) + list(passthrough)
    # Nextflow >= 26.04 parses with the strict v2 syntax by default, which nextflow/ does not use yet
    env = dict(os.environ, WGS_SKIM_CONDA_PREFIX=str(conda_prefix), NXF_SYNTAX_PARSER="v1")
    return execute(cmd, cwd=str(out), env=env, show=out)   # .nextflow/ (the -resume cache) lives there

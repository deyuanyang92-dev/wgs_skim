"""Running external programs, shared by the workflow scripts: commands with a log, bash pipelines, BLAST.

    from runtools import run, sh, blastn
    run(["get_organelle_from_reads.py", ...], log="x.log")          # stdout + stderr to the log, returns the exit code
    hits = blastn("query.fasta", "assembly.fasta", fields="sseqid pident length")   # [{'sseqid': ..., ...}]
Standard library only: also usable from the Python 2-era tool environments' callers (rule mitofinder runs its
script with Snakemake's Python)."""
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# a target FASTA bigger than this is searched through a BLAST database (-subject keeps it all in memory and is
# single-threaded: a 2 Gb genome assembly does not finish)
SUBJECT_MAX_BYTES = 50_000_000


def tool_env(exe=None, prepend=()):
    """os.environ with the directory of `exe` (a tool's own environment: GetOrganelle finds bowtie2 / SPAdes /
    BLAST beside itself) and `prepend` directories first on PATH"""
    env = dict(os.environ)
    dirs = [str(d) for d in prepend]
    if exe and Path(exe).parent != Path("."):
        dirs.append(str(Path(exe).parent))
    if dirs:
        env["PATH"] = os.pathsep.join(dirs + [env.get("PATH", "")])
    return env


def run(cmd, log=None, env=None, check=False, cwd=None, append=False):
    """run a command (list); stdout + stderr to `log` (path or open file; default: this process's stderr).
    Returns the exit code; check=True raises on failure."""
    cmd = [str(c) for c in cmd]
    if log is None:
        rc = subprocess.run(cmd, stdout=sys.stderr, stderr=subprocess.STDOUT, env=env, cwd=cwd).returncode
    elif hasattr(log, "write"):
        log.flush()
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=cwd).returncode
    else:
        with open(log, "a" if append else "w") as fh:
            rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=cwd).returncode
    if check and rc:
        raise RuntimeError(f"{Path(cmd[0]).name} failed (exit {rc})" + (f"; see {log}" if isinstance(log, (str, Path))
                                                                         else ""))
    return rc


def sh(cmd, check=True, env=None):
    """a bash pipeline with pipefail (stderr to this process's stderr)"""
    return subprocess.run(f"set -o pipefail; {cmd}", shell=True, check=check, executable="/bin/bash",
                          stderr=sys.stderr, env=env).returncode


def output(cmd, env=None, check=True):
    """stdout of a command as text"""
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, env=env, check=check).stdout


def blastn(query, target, fields="qseqid sseqid pident length", evalue=1e-10, threads=1, exe="blastn",
           extra=(), db=False):
    """blastn of `query` against the FASTA `target` -> [{field: value}] (text values, in BLAST's order).
    A large target (or db=True) is searched through a temporary BLAST database, removed afterwards."""
    fields = fields.split() if isinstance(fields, str) else list(fields)
    cmd = [exe, "-query", str(query), "-evalue", str(evalue), "-outfmt", "6 " + " ".join(fields), *map(str, extra)]
    big = db or Path(target).stat().st_size > SUBJECT_MAX_BYTES
    tmp = tempfile.mkdtemp(prefix="wgs_skim_blast_", dir=os.environ.get("TMPDIR")) if big else None
    try:
        if big:
            mk = str(Path(exe).parent / "makeblastdb") if Path(exe).parent != Path(".") else "makeblastdb"
            subprocess.run([mk, "-in", str(target), "-dbtype", "nucl", "-out", f"{tmp}/db"], check=True,
                           stdout=subprocess.DEVNULL, stderr=sys.stderr)
            cmd += ["-db", f"{tmp}/db", "-num_threads", str(threads)]
        else:
            cmd += ["-subject", str(target)]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)
    return [dict(zip(fields, line.split("\t"))) for line in out.splitlines() if line.strip()]


def aligned_bp(hits, key="sseqid", min_pident=0):
    """{key value: aligned bp} summed over the HSPs with pident >= min_pident (hits from blastn(... pident length))"""
    cov = {}
    for h in hits:
        if float(h["pident"]) >= min_pident:
            cov[h[key]] = cov.get(h[key], 0) + int(h["length"])
    return cov

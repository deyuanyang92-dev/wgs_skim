"""Where things are: the workflow files, the user's settings (~/.wgs_skim), the software environments, the machine."""
import os
import sys
from pathlib import Path

HOME = Path(os.environ.get("WGS_SKIM_HOME", Path.home() / ".wgs_skim"))
USER_CONFIG = HOME / "config.yaml"
TAXA = ["annelid", "mollusc", "fish", "crustacean"]


def default_conda_prefix():
    """WGS_SKIM_CONDA_PREFIX, else `conda_prefix:` in ~/.wgs_skim/config.yaml, else ~/.wgs_skim/conda."""
    if os.environ.get("WGS_SKIM_CONDA_PREFIX"):
        return os.environ["WGS_SKIM_CONDA_PREFIX"]
    if USER_CONFIG.exists():
        for line in USER_CONFIG.read_text().splitlines():
            if line.startswith("conda_prefix:"):
                return line.split(":", 1)[1].split("#")[0].strip()
    return HOME / "conda"


def workflow_root():
    """Directory holding workflow/ config/ resources/: the installed package data, or the git checkout."""
    here = Path(__file__).resolve().parent
    for root in (here, here.parent.parent):
        if (root / "workflow" / "Snakefile").is_file():
            return root
    sys.exit("[wgs_skim] cannot find workflow/Snakefile (broken installation?)")


def scripts_dir():
    """workflow/scripts: the library every rule uses (see wgs_skim.lib)"""
    return workflow_root() / "workflow" / "scripts"


def machine_mem_gb():
    """total memory of this machine (GB), 0 if unknown"""
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) / 1e6
    except OSError:
        pass
    return 0


def job_mem_gb(asked=0):
    """memory all jobs may use together: what was asked, else 90 % of the machine"""
    return asked or int(machine_mem_gb() * 0.9)

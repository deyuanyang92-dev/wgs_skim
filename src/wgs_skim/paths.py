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


def presets():
    """{name: preset file}: the built-in ones (config/presets) and the user's (~/.wgs_skim/presets, which win)"""
    found = {}
    for d in (workflow_root() / "config" / "presets", HOME / "presets"):
        found |= {p.stem: p for p in sorted(Path(d).glob("*.yaml"))} if Path(d).is_dir() else {}
    return found


def preset_file(name):
    """a preset by name, or a YAML file path"""
    if Path(name).is_file():
        return Path(name).resolve()
    known = presets()
    if name not in known:
        sys.exit(f"[wgs_skim] unknown preset '{name}': {', '.join(known)} (or a YAML file; your own presets go in "
                 f"{HOME / 'presets'}/<name>.yaml)")
    return known[name]


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

"""The workflow's Python library (workflow/scripts) as modules of the package:

    from wgs_skim.lib import seqtools, samples, blast_taxonomy
    recs = seqtools.read_fasta("mt.fasta")
    pairs = samples.find_pairs("rawdata/")
    from wgs_skim.lib.mitofinder_check import check

Every module there imports without side effects (docs/DEVELOPING.md); they import each other by plain name
(`from seqtools import ...`, as Snakemake runs them), so the scripts directory is put on sys.path here."""
import importlib
import sys

from .paths import scripts_dir

_DIR = str(scripts_dir())
if _DIR not in sys.path:
    sys.path.insert(0, _DIR)
__path__ = [_DIR]                     # `import wgs_skim.lib.<module>` / `from wgs_skim.lib.<module> import ...`


def modules():
    """names of the library modules"""
    return sorted(p.stem for p in scripts_dir().glob("*.py") if not p.stem.startswith("_"))


def __getattr__(name):
    if name in modules():
        mod = importlib.import_module(name)
        sys.modules[f"{__name__}.{name}"] = mod
        return mod
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

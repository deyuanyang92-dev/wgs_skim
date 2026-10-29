"""wgs_skim: genome-skimming workflow (Snakemake) for animal mitogenomes, nuclear rDNA and genome surveys.

From Python: wgs_skim.run(...), wgs_skim.tool(...), wgs_skim.command(...) (wgs_skim.api); the library of the
workflow scripts: wgs_skim.lib (e.g. `from wgs_skim.lib import seqtools`)."""
__version__ = "0.1.0"


def __getattr__(name):            # wgs_skim.run / tool / command without importing argparse on `import wgs_skim`
    if name in ("run", "tool", "command", "argv"):
        from . import api
        return getattr(api, name)
    raise AttributeError(f"module 'wgs_skim' has no attribute {name!r}")

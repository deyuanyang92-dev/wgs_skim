#!/usr/bin/env python3
"""Run a script written for Snakemake's `script:` directive (workflow/scripts/*.py) from a Nextflow process.

    smk_shim.py snakemake.json SCRIPT.py

The JSON holds input / output / params / log / threads / wildcards / scriptdir; the script gets them as the global
`snakemake` object, with the namedlist access it uses in Snakemake (sm.input.fasta, sm.input.get("x"), sm.log[0]),
so one copy of every script serves both workflows.
"""
import json
import runpy
import sys
from pathlib import Path


class NamedList(list):
    """Snakemake's namedlist: the values in order ([0], iteration) and by name (.name, .get(name))."""

    def __init__(self, data=None):
        self._names = dict(data) if isinstance(data, dict) else {}
        values = data.values() if isinstance(data, dict) else (data or [])
        # positional items are the files of all names, flattened, as in Snakemake
        super().__init__(x for v in values for x in (v if isinstance(v, list) else [v]))

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            return self._names[name]
        except KeyError:
            raise AttributeError(f"{name!r} not given (keys: {', '.join(self._names)})") from None

    def get(self, name, default=None):
        return self._names.get(name, default)

    def keys(self):
        return self._names.keys()

    def items(self):
        return self._names.items()


class Snakemake:
    def __init__(self, d):
        for k in ("input", "output", "params", "log", "wildcards", "resources"):
            setattr(self, k, NamedList(d.get(k)))
        self.threads = int(d.get("threads", 1))
        self.scriptdir = d["scriptdir"]
        self.rule = d.get("rule", "")
        self.config = d.get("config", {})


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    spec, script = sys.argv[1:]
    sm = Snakemake(json.load(open(spec)))
    # as Snakemake: output / log directories exist; a log the script never writes (summary.py) is still there
    for f in list(sm.output) + list(sm.log):
        Path(f).parent.mkdir(parents=True, exist_ok=True)
    for f in sm.log:
        open(f, "a").close()
    sys.argv = [script]
    runpy.run_path(script, init_globals={"snakemake": sm}, run_name="__main__")


if __name__ == "__main__":
    main()

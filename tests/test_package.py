"""The Python package: the API builds the same command lines as the CLI; the workflow library is importable."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import wgs_skim  # noqa: E402
from wgs_skim import tools  # noqa: E402


def test_api_argv_matches_the_parser():
    assert wgs_skim.argv("busco", {"input": ["a/", "b/"], "mode": "transcriptome", "mem_gb": 50, "dry_run": True,
                                   "lineage": None}) == \
        ["busco", "--input", "a/", "b/", "--mode", "transcriptome", "--mem-gb", "50", "--dry-run"]
    assert wgs_skim.argv("run", {"reads_dir": "r"}, snakemake=["--rerun-incomplete"]) == \
        ["run", "--reads_dir", "r", "--", "--rerun-incomplete"]
    with pytest.raises(TypeError):
        wgs_skim.argv("busco", {"lineag": "x"})
    with pytest.raises(ValueError):
        wgs_skim.tool("nosuchtool")


def test_lib_exposes_the_workflow_scripts():
    from wgs_skim.lib import seqtools
    from wgs_skim.lib.mitofinder_check import family
    import wgs_skim.lib as lib
    assert seqtools.rc("AAC") == "GTT" and family("A; Bidae; C") == "Bidae"
    assert {"seqtools", "samples", "runtools", "tables", "collect_fasta", "blast_taxonomy"} <= set(lib.modules())


def test_every_tool_has_a_step_and_parses():
    p = wgs_skim.cli.parser() if hasattr(wgs_skim, "cli") else __import__("wgs_skim.cli").cli.parser()
    for name, t in tools.TOOLS.items():
        assert t.inputs in ("reads", "fasta", "both") and t.step
        flag = ["-i", "x.fasta"] if t.inputs == "fasta" else ["-r", "reads/"]
        args = p.parse_args([name, *flag] + (["--ref", "r.gb"] if name == "mitofinder" else []))
        assert args.tool == name

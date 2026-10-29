"""wgs_skim from Python: the same commands, options and behaviour as the command line.

    import wgs_skim
    wgs_skim.run(reads_dir="rawdata/", outdir="res/", steps="spades-busco", cores=32)
    wgs_skim.tool("busco", input=["asm/"], outdir="busco/", mode="transcriptome")
    wgs_skim.tool("mitofinder", reads_dir="rnaseq/", ref="related.gb", assembler="rnaspades")
    wgs_skim.command("compare", results="res/", outdir="compare/")
Keyword arguments are the command's options by their names in `wgs_skim <command> -h` (dashes as underscores:
mem_gb, dry_run, split_by_id); True switches a flag on, a list gives several values. `snakemake=[...]` is passed
on to Snakemake (what follows `--` on the command line). Returns the exit code (0 = success).
The library functions themselves: wgs_skim.lib (workflow/scripts)."""
from .cli import dispatch, parser


def _options(command):
    """{dest: option string} of a command's parser"""
    sub = next(a for a in parser()._actions if a.dest == "command").choices[command]
    return {a.dest: max(a.option_strings, key=len) for a in sub._actions if a.option_strings}


def argv(command, options, snakemake=()):
    """the command line of `command` with these options (see the module doc)"""
    known, out = _options(command), [command]
    for k, v in options.items():
        if k not in known:
            raise TypeError(f"wgs_skim {command}: no option {k!r} (options: {', '.join(sorted(known))})")
        if v is None or v is False:
            continue
        if v is True:
            out.append(known[k])
        elif isinstance(v, (list, tuple)):
            out += [known[k], *map(str, v)]
        else:
            out += [known[k], str(v)]
    return out + (["--", *map(str, snakemake)] if snakemake else [])


def command(name, snakemake=(), **options):
    """any wgs_skim command (run, a tool, compare, nuclear, circularize, samples ...)"""
    return dispatch(argv(name, options, snakemake))


def run(snakemake=(), **options):
    """the workflow: wgs_skim run"""
    return command("run", snakemake, **options)


def tool(name, snakemake=(), **options):
    """one tool for every sample: wgs_skim <tool> (wgs_skim.tools.TOOLS)"""
    from .tools import TOOLS
    if name not in TOOLS:
        raise ValueError(f"unknown tool {name!r}: {', '.join(TOOLS)}")
    return command(name, snakemake, **options)

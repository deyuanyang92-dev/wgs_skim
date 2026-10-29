"""GetOrganelle runs, shared by every script that starts or reads one: the command (getorganelle_cmd), the final
sequence(s) (result, best), renaming to `<sample>|<locus>|<n> <topology> length=<bp>`, cleaning (clean).

GetOrganelle writes one *.path_sequence.fasta per alternative path (graph1.1, graph1.2 ... differ only in the
orientation of flip-flop regions); the first one is used. Circular sequences have "circular" in the header.
`clean()` removes a run's bulky intermediates (the user's remove_getorganelle_templefile*.py): the extended reads,
the extension's SPAdes run, the seed index and, after a failed run, the unzipped reads and their index - the result,
graphs (gfa / fastg) and log stay.   Usage: collect_fasta.py DIR SAMPLE LOCUS > out.fasta | collect_fasta.py --clean DIR
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from seqtools import is_circular, read_fasta, write_fasta  # noqa: E402,F401

EXE = "get_organelle_from_reads.py"


def result_fastas(getorganelle_dir):
    """the *.path_sequence.fasta files: complete graph first, scaffolds only as fallback"""
    d = Path(getorganelle_dir)
    fas = sorted(d.glob("*.path_sequence.fasta")) if d.is_dir() else []
    complete = [f for f in fas if "scaffolds" not in f.name]
    return complete or fas


def result(path):
    """(fasta, gfa or None, number of alternative paths) of a GetOrganelle output directory (or a FASTA itself);
    (None, None, 0) without a result. Alternative paths (graph1.1, graph1.2 ...) differ only in the orientation of
    flip-flop regions; the first is used."""
    p = Path(path)
    if p.is_file():
        return p, None, 1
    fas = result_fastas(p)
    if not fas:
        return None, None, 0
    graph = fas[0].name.split(".path_sequence")[0].rsplit(".", 1)[0]       # <...>.graph1
    gfa = p / f"{graph}.selected_graph.gfa"
    return fas[0], (gfa if gfa.exists() else None), sum(1 for f in fas if f.name.startswith(graph + "."))


def result_fasta(getorganelle_dir):
    """the final FASTA of a GetOrganelle run, or None"""
    return result(getorganelle_dir)[0] if Path(getorganelle_dir).is_dir() else None


def best(getorganelle_dir):
    """(records, circular) of the first result file, or ([], False)."""
    fasta = result_fasta(getorganelle_dir)
    if not fasta:
        return [], False
    recs = read_fasta(fasta)
    return recs, any(is_circular(n) for n, _ in recs)


def getorganelle_cmd(r1, r2, out, target="animal_mt", threads=8, seed=None, rounds=None, kmers=None,
                     max_reads=None, config_dir=None, genes=None, extra=(), exe=EXE):
    """get_organelle_from_reads.py command (list): -F target, optional seed (-s), label genes (--genes), -R rounds,
    -k k-mers, --max-reads per file (None = all reads), --config-dir; always --overwrite and no read reduction
    for coverage (the data cap is --max-reads)"""
    cmd = [exe, "-1", r1, "-2", r2, "-o", out, "-F", target, "-t", threads, "--overwrite",
           "--max-reads", max_reads or "inf", "--reduce-reads-for-coverage", "inf"]
    for flag, v in (("-s", seed), ("--genes", genes), ("-R", rounds), ("-k", kmers), ("--config-dir", config_dir)):
        if v not in (None, ""):
            cmd += [flag, v]
    return [str(c) for c in cmd + list(extra)]


def clean(getorganelle_dir):
    """remove the intermediates of a GetOrganelle run (see module doc); returns the bytes freed"""
    d = Path(getorganelle_dir)
    if not d.is_dir():
        return 0
    gone = list(d.glob("extended_*.fq")) + list(d.glob("[0-9]-*.fastq")) + list(d.glob("*temp.indices*")) + \
        [d / "extended_spades", d / "seed"]
    freed = 0
    for x in gone:
        if x.is_dir():
            freed += sum(f.stat().st_size for f in x.rglob("*") if f.is_file())
            shutil.rmtree(x, ignore_errors=True)
        elif x.exists():
            freed += x.stat().st_size
            x.unlink()
    return freed


def renamed(recs, sample, locus):
    out = []
    for i, (name, seq) in enumerate(recs, 1):
        topo = "circular" if is_circular(name) else "linear"
        out.append((f"{sample}|{locus}|{i} {topo} length={len(seq)}", seq))
    return out


def main():
    if sys.argv[1] == "--clean":
        print(f"{clean(sys.argv[2]) / 1e6:.0f} MB removed from {sys.argv[2]}", file=sys.stderr)
        return
    getorganelle_dir, sample, locus = sys.argv[1:4]
    recs, _ = best(getorganelle_dir)
    write_fasta(renamed(recs, sample, locus), sys.stdout)


if __name__ == "__main__":
    main()

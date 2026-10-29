"""Pick the final sequence(s) of a GetOrganelle run and rename them `<sample>|<locus>|<n> <topology> length=<bp>`.

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
from seqtools import read_fasta, write_fasta  # noqa: E402,F401


def result_fastas(getorganelle_dir):
    d = Path(getorganelle_dir)
    fas = sorted(d.glob("*.path_sequence.fasta")) if d.is_dir() else []
    complete = [f for f in fas if "scaffolds" not in f.name]   # complete graph first, scaffolds only as fallback
    return complete or fas


def best(getorganelle_dir):
    """(records, circular) of the first result file, or ([], False)."""
    fas = result_fastas(getorganelle_dir)
    if not fas:
        return [], False
    recs = read_fasta(fas[0])
    return recs, any("circular" in n.lower() for n, _ in recs)


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
        topo = "circular" if "circular" in name.lower() else "linear"
        out.append((f"{sample}|{locus}|{i} {topo} length={len(seq)}", seq))
    return out


if __name__ == "__main__":
    if sys.argv[1] == "--clean":
        print(f"{clean(sys.argv[2]) / 1e6:.0f} MB removed from {sys.argv[2]}", file=sys.stderr)
        sys.exit(0)
    getorganelle_dir, sample, locus = sys.argv[1:4]
    recs, _ = best(getorganelle_dir)
    write_fasta(renamed(recs, sample, locus), sys.stdout)

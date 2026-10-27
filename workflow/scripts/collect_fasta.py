"""Pick the final sequence(s) of a GetOrganelle run and rename them `<sample>|<locus>|<n> <topology> length=<bp>`.

GetOrganelle writes one *.path_sequence.fasta per alternative path (graph1.1, graph1.2 ... differ only in the
orientation of flip-flop regions); the first one is used. Circular sequences have "circular" in the header.
"""
import sys
from pathlib import Path


def read_fasta(path):
    recs, name, seq = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(seq)))
            name, seq = line[1:], []
        elif line:
            seq.append(line)
    if name is not None:
        recs.append((name, "".join(seq)))
    return recs


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


def renamed(recs, sample, locus):
    out = []
    for i, (name, seq) in enumerate(recs, 1):
        topo = "circular" if "circular" in name.lower() else "linear"
        out.append((f"{sample}|{locus}|{i} {topo} length={len(seq)}", seq))
    return out


def write_fasta(recs, fh):
    for name, seq in recs:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")


if __name__ == "__main__":
    getorganelle_dir, sample, locus = sys.argv[1:4]
    recs, _ = best(getorganelle_dir)
    write_fasta(renamed(recs, sample, locus), sys.stdout)

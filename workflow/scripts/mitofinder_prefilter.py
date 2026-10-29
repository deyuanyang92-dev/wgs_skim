"""The scaffolds MitoFinder needs to see from a whole-genome assembly: >= min_len bp and k-mer coverage
>= factor x the length-weighted median coverage (mitochondrial scaffolds are tens to hundreds of times deeper than
nuclear ones). As the user's spades_mitofinder.py: MJ20260820041 CHY006 2 Gb SPAdes assembly -> 3,871 scaffolds,
MitoFinder 5 min and one mitochondrial contig (so that its circularization check runs). Assemblies without SPAdes
coverage in the headers (NODE_x_length_L_cov_C) keep every record >= min_len.
    mitofinder_prefilter.py IN.fasta OUT.fasta [--min-len 1000] [--factor 5]   (factor 0: length filter only)
"""
import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from seqtools import iter_fasta  # noqa: E402

COV = re.compile(r"_length_(\d+)_cov_([\d.]+)")


def weighted_median(pairs):
    """median coverage weighted by length: half of the assembled bases are at or below it"""
    total, acc = sum(n for n, _ in pairs), 0
    for n, c in sorted(pairs, key=lambda x: x[1]):
        acc += n
        if acc >= total / 2:
            return c
    return 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta")
    ap.add_argument("out")
    ap.add_argument("--min-len", type=int, default=1000)
    ap.add_argument("--factor", type=float, default=5)
    a = ap.parse_args()
    covs = []
    for name, seq in iter_fasta(a.fasta, upper=False):
        m = COV.search(name)
        if len(seq) >= a.min_len:
            covs.append((len(seq), float(m.group(2)) if m else None))
    with_cov = [(n, c) for n, c in covs if c is not None]
    use_cov = a.factor > 0 and with_cov and len(with_cov) == len(covs)
    cut = a.factor * weighted_median(with_cov) if use_cov else 0
    kept = 0
    with open(a.out, "w") as out:
        for name, seq in iter_fasta(a.fasta, upper=False):
            m = COV.search(name)
            if len(seq) >= a.min_len and (not use_cov or float(m.group(2)) >= cut):
                out.write(f">{name}\n{seq}\n")
                kept += 1
    how = f"coverage >= {cut:.1f} ({a.factor:g} x weighted median)" if use_cov else "no coverage filter"
    print(f"{len(covs)} records >= {a.min_len} bp, {how}: {kept} kept")


if __name__ == "__main__":
    main()

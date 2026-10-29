"""The scaffolds MitoFinder needs to see from a whole-genome assembly: >= min_len bp and k-mer coverage
>= factor x the length-weighted median coverage (mitochondrial scaffolds are usually tens to hundreds of times
deeper than nuclear ones), plus every scaffold the reference mitogenome hits (--ref: blastn dc-megablast,
e <= 1e-10), whatever its coverage. As the user's spades_mitofinder.py: MJ20260820041 CHY006 2 Gb SPAdes
assembly -> 3,871 scaffolds, MitoFinder 5 min and one mitochondrial contig (so that its circularization check runs).
The reference rescue: a sub-sampled genome skim (wgs_skim_test cwsc010, 0.3 Gb) has its mitochondrial scaffolds at
1.5-2.9x, under 5 x the median (6.3x) - coverage alone dropped all of them and MitoFinder found nothing.
Assemblies without SPAdes coverage in the headers (NODE_x_length_L_cov_C) keep every record >= min_len.
    mitofinder_prefilter.py IN.fasta OUT.fasta [--min-len 1000] [--factor 5] [--ref ref.gb] [--threads 4]
    (factor 0: length filter only)
"""
import argparse
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from runtools import blastn  # noqa: E402
from seqtools import genbank_seqs, iter_fasta, write_fasta  # noqa: E402

COV = re.compile(r"_length_(\d+)_cov_([\d.]+)")


def weighted_median(pairs):
    """median coverage weighted by length: half of the assembled bases are at or below it"""
    total, acc = sum(n for n, _ in pairs), 0
    for n, c in sorted(pairs, key=lambda x: x[1]):
        acc += n
        if acc >= total / 2:
            return c
    return 0.0


def coverage(name):
    m = COV.search(name)
    return float(m.group(2)) if m else None


def ref_hits(fasta, ref, threads=4):
    """record names (first word) of `fasta` that the reference mitogenome(s) (GenBank or FASTA) hit"""
    recs = genbank_seqs(ref) if Path(ref).read_text(errors="replace").startswith("LOCUS") else None
    with tempfile.NamedTemporaryFile("w", suffix=".fasta") as q:
        if recs is None:
            q.write(Path(ref).read_text())
        else:
            write_fasta(recs, q)
        q.flush()
        hits = blastn(q.name, fasta, "sseqid", 1e-10, threads, extra=["-task", "dc-megablast"])
    return {h["sseqid"] for h in hits}


def select(fasta, min_len=1000, factor=5, ref=None, threads=4):
    """(names kept, description): the module doc's rule"""
    covs = [(len(seq), coverage(name)) for name, seq in iter_fasta(fasta, upper=False) if len(seq) >= min_len]
    with_cov = [(n, c) for n, c in covs if c is not None]
    use_cov = bool(factor > 0 and with_cov and len(with_cov) == len(covs))
    cut = factor * weighted_median(with_cov) if use_cov else 0
    hits = ref_hits(fasta, ref, threads) if ref and use_cov else set()
    keep, rescued = set(), 0
    for name, seq in iter_fasta(fasta, upper=False):
        if len(seq) < min_len:
            continue
        if not use_cov or coverage(name) >= cut:
            keep.add(name)
        elif name.split()[0] in hits:
            keep.add(name)
            rescued += 1
    how = f"coverage >= {cut:.1f} ({factor:g} x weighted median)" if use_cov else "no coverage filter"
    if ref and use_cov:
        how += f" or hit by the reference ({rescued} below the coverage cut)"
    return keep, f"{len(covs)} records >= {min_len} bp, {how}: {len(keep)} kept"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta")
    ap.add_argument("out")
    ap.add_argument("--min-len", type=int, default=1000)
    ap.add_argument("--factor", type=float, default=5)
    ap.add_argument("--ref", help="reference mitogenome (GenBank or FASTA): scaffolds it hits are kept")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    keep, how = select(a.fasta, a.min_len, a.factor, a.ref, a.threads)
    with open(a.out, "w") as out:
        out.writelines(f">{n}\n{s}\n" for n, s in iter_fasta(a.fasta, upper=False) if n in keep)
    print(how)


if __name__ == "__main__":
    main()

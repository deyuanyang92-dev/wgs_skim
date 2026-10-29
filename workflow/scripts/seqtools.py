"""Sequences: FASTA I/O, reverse complement, lengths / N50, GenBank sequence, end-overlap detection - shared by
the workflow scripts (import them from here, not again).

Every script in this directory imports with no side effects (its work is in main() / main(sm)), so any function
can be reused: `sys.path.insert(0, <workflow/scripts>)` - a Snakemake `script:` adds snakemake.scriptdir."""
import gzip
from collections import defaultdict
from pathlib import Path


def iter_fasta(path, first_word=False, upper=True):
    """(header, sequence) one record at a time (large assemblies), plain or .gz; header = the whole line after '>'
    (first_word: up to the first space)"""
    opener = gzip.open if str(path).endswith(".gz") else open
    name, seq = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = (line[1:].split()[0] if first_word and line[1:].strip() else line[1:]), []
            elif name is not None and line:
                seq.append(line.upper() if upper else line)
    if name is not None:
        yield name, "".join(seq)


def read_fasta(path, first_word=False, upper=True):
    """[(header, SEQUENCE)] (see iter_fasta)"""
    return list(iter_fasta(path, first_word, upper))


def parse_fasta(text, upper=True):
    """[(header, sequence)] of FASTA text (a program's output)"""
    recs = []
    for line in text.splitlines():
        if line.startswith(">"):
            recs.append([line[1:].strip(), []])
        elif recs and line.strip():
            recs[-1][1].append(line.strip().upper() if upper else line.strip())
    return [(n, "".join(s)) for n, s in recs]


def write_fasta(recs, fh, width=80):
    """write (header, sequence) records to an open file or a path; width 0 = one line per sequence"""
    if isinstance(fh, (str, Path)):
        with open(fh, "w") as out:
            return write_fasta(recs, out, width)
    for name, seq in recs:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width or max(len(seq), 1)):
            fh.write(seq[i:i + (width or len(seq))] + "\n")


COMP = str.maketrans("ACGTNRYKMSWBDHVacgtnrykmswbdhv", "TGCANYRMKSWVHDBtgcanyrmkswvhdb")


def rc(s):
    """reverse complement (IUPAC codes, case kept)"""
    return s.translate(COMP)[::-1]


def fasta_lengths(path):
    """[length of every non-empty record] (streamed)"""
    return [len(s) for _, s in iter_fasta(path, upper=False) if s]


def n50(lens):
    """N50 of a list of lengths ('' when empty)"""
    if not lens:
        return ""
    half, acc = sum(lens) / 2, 0
    for x in sorted(lens, reverse=True):
        acc += x
        if acc >= half:
            return x


def is_circular(header):
    """GetOrganelle / wgs_skim headers mark a circular sequence with 'circular'"""
    return "circular" in header.lower()


def genbank_seqs(path):
    """[(accession or LOCUS name, SEQUENCE)] of a GenBank flat file, without Biopython (any Python 3 environment)"""
    recs, name, seq, in_seq = [], None, [], False
    for line in open(path):
        if line.startswith("LOCUS"):
            name, seq, in_seq = line.split()[1], [], False
        elif line.startswith("ORIGIN"):
            in_seq = True
        elif line.startswith("//"):
            if name is not None:
                recs.append((name, "".join(seq).upper()))
            name, in_seq = None, False
        elif in_seq:
            seq.append("".join(c for c in line if c.isalpha()))
    return recs


def find_circle(s, k=21, span=2000, end_window=5000, min_run=25, min_len=1000):
    """End overlap of a circular molecule assembled as a line (improved MITObim circules.py): k-mers of the first
    `span` bp that recur in the last `end_window` bp vote for the circle length P; the longest stretch with
    s[i] == s[i+P] (>= min_run, ragged ends tolerated) confirms it. Returns dict(P, overlap, run, votes, anchor)
    with the circle = s[anchor:anchor+P], or None."""
    L = len(s)
    if L < min_len:
        return None
    tail = defaultdict(list)
    for j in range(max(span, L - end_window), L - k + 1):
        tail[s[j:j + k]].append(j)
    votes = defaultdict(int)
    for i in range(min(span, L - k + 1)):
        for j in tail.get(s[i:i + k], ()):
            votes[j - i] += 1
    best = None
    for P, v in sorted(votes.items(), key=lambda x: -x[1])[:10]:
        if P < min_len:
            continue
        run = start = cur = 0
        for i in range(L - P):
            cur = cur + 1 if s[i] == s[i + P] else 0
            if cur > run:
                run, start = cur, i - cur + 1
        if run >= min_run and (best is None or (run, v) > (best["run"], best["votes"])):
            best = dict(P=P, overlap=L - P, run=run, votes=v, anchor=start + run // 2)
    return best

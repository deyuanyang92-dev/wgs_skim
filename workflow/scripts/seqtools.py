"""FASTA I/O and end-overlap detection shared by the workflow scripts (import them from here, not again).

Every script in this directory imports with no side effects (its work is in main() / main(sm)), so any function
can be reused: `sys.path.insert(0, <workflow/scripts>)` - a Snakemake `script:` adds snakemake.scriptdir."""
import gzip
from collections import defaultdict


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


def write_fasta(recs, fh):
    for name, seq in recs:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")


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

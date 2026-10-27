"""FASTA I/O and end-overlap detection shared by the workflow scripts."""
from collections import defaultdict


def read_fasta(path):
    """[(header, SEQUENCE)]; the header is the whole line after '>'."""
    recs = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            recs.append([line[1:], []])
        elif recs and line:
            recs[-1][1].append(line.upper())
    return [(n, "".join(q)) for n, q in recs]


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

"""Circularize mitochondrial contigs from a de novo assembly (MitoFinder / SPAdes) without reads.

Why: MitoFinder and MitoHiFi decide circularity by a self-BLAST hit >= 220 bp within 40 bp of both ends
(circularizationCheck.py, Schomaker 2014). Short-read assemblers leave an overlap of only ~k-1 bp, often with a few
ragged bases (MJ20260820041 CHY006: 58 bp, of which the outer 24 bp disagree), so circular contigs are reported as
"not found". Four checks, strongest first:

  graph    the contig's path in the SPAdes assembly graph (scaffolds.paths + assembly_graph_with_scaffolds.gfa)
           links its last segment back to its first: the molecule closes in the graph itself. The circle is the
           contig minus the link overlap. Segments used more than once in the path are tandem repeat units whose
           copy number the graph cannot resolve (CHY007: control-region unit used 4-5 times) - reported.
  walk     the path end does not link to its start, but a short walk through the graph does (<= 12 segments,
           each >= 0.2x the mt depth): the gap SPAdes left in a repeat region is filled from the graph. Of the walks found, the one with the highest
           minimum depth is used; its segments and depths are reported. MitoFinder-trimmed contigs are accepted
           when they lie inside the scaffold the path spells; the full scaffold is then circularized.
           CHY008: 3 control-region segments (depth 404-601x, contig 998x) close a 17048 bp circle, identical
           (99.99%, one 2-bp indel) to the circle GetOrganelle assembled and reads confirmed.
           Tandem repeats: a segment whose depth says it occurs more often than the path uses it gets the missing
           turns of its loop (copies = round(segment depth / depth of the contig's longest segment)). CHY006: a
           94-bp control-region unit at 2500x vs 858x -> 2 more turns -> 17052 bp, identical (100%) to GetOrganelle.
           Copy numbers that disagree with depth are reported ("repeat copies uncertain"): CHY007's scaffold uses
           a unit 4-5 times at ~1.0-1.2x depth (17228 bp). Short reads cannot count tandem copies exactly.
  overlap  the ends overlap (k-mer vote + longest agreeing stretch >= 25 bp, ragged ends tolerated; seqtools).
           Tried after the graph: an overlap inside a tandem repeat collapses repeat copies (CHY006: 58-bp overlap
           gave 16770 bp; the walk and GetOrganelle 17052 bp).
  join     several contigs (MitoFinder often returns pieces): joined end to end through graph links or exact
           end overlaps (>= 20 bp); circular if the chain closes.

Circular by sequence is not proof of a circular molecule: check junction-spanning reads (rule topology) for that.

Usage:  circularize.py -i mt_contigs.fasta -o out/ [--spades spades_dir] [--prefix sample]
        -> out/<prefix>.circularized.fasta, out/<prefix>.circularize.tsv
"""
import argparse
import csv
import itertools
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from seqtools import find_circle, read_fasta, write_fasta  # noqa: E402

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
FLIP = {"+": "-", "-": "+"}


def rc(s):
    return s.translate(COMP)[::-1]


def read_paths(path_file, names):
    """{contig: [[(segment, orientation), ...] per gap-free chunk]} for the wanted SPAdes contig names."""
    out, cur = {}, None
    for line in open(path_file):
        line = line.strip()
        if line.startswith("NODE_"):
            cur = line if line in names else None
            if cur:
                out[cur] = [[]]
        elif cur:
            out[cur][-1] += [(x[:-1], x[-1]) for x in line.rstrip(";").split(",") if x]
            if line.endswith(";"):                     # scaffold gap after this chunk
                out[cur].append([])
    return {k: [c for c in v if c] for k, v in out.items()}


def read_links(gfa, segs):
    """Lengths of the wanted segments and every link touching them: {(a, ao, b, bo): overlap} in both directions.
    Streams the GFA (several GB for a genome)."""
    length, links = {}, {}
    with open(gfa) as fh:
        for line in fh:
            if line[0] == "S":
                t = line.split("\t", 3)
                if t[1] in segs:
                    length[t[1]] = len(t[2])
            elif line[0] == "L":
                _, a, ao, b, bo, cig = line.rstrip("\n").split("\t")[:6]
                if a in segs or b in segs:
                    ov = int(re.match(r"(\d+)", cig).group(1)) if cig[0].isdigit() else 0
                    links[(a, ao, b, bo)] = ov
                    links[(b, FLIP[bo], a, FLIP[ao])] = ov
    return length, links


def read_graph(gfa):
    """Whole graph: ({(segment, orientation): [(next, orientation, overlap)]}, {segment: depth}). ~1 min and ~6 GB
    for a 2.5 GB GFA (9 M links); read only when a contig needs a walk."""
    adj, depth = {}, {}
    with open(gfa) as fh:
        for line in fh:
            if line[0] == "L":
                _, a, ao, b, bo, cig = line.rstrip("\n").split("\t")[:6]
                ov = int(re.match(r"(\d+)", cig).group(1)) if cig[0].isdigit() else 0
                adj.setdefault((a, ao), []).append((b, bo, ov))
                adj.setdefault((b, FLIP[bo]), []).append((a, FLIP[ao], ov))
            elif line[0] == "S":
                t = line.split("\t", 2)
                m = re.search(r"\tDP:f:([\d.eE+-]+)", t[2])
                depth[t[1]] = float(m.group(1)) if m else 0.0
    return adj, depth


def read_sequences(gfa, names):
    with open(gfa) as fh:
        return {t[1]: t[2] for t in (line.split("\t", 3) for line in fh if line[0] == "S") if t[1] in names}


def walks(adj, ok, start, goal, max_steps=12, budget=200000):
    """Paths start -> ... -> goal through segments ok() accepts, without revisits: [[(seg, ori), ...]] (start and
    goal excluded)."""
    found, todo = [], [[start]]
    while todo and budget:
        budget -= 1
        p = todo.pop()
        for b, bo, _ in adj.get(p[-1], []):
            if (b, bo) == goal:
                found.append(p[1:])
            elif len(p) <= max_steps and ok(b) and (b, bo) not in p:
                todo.append(p + [(b, bo)])
    return found


def add_repeat_copies(circ, adj, copies):
    """Tandem repeats: a segment used fewer times in the circle than its depth says (copies(s) = round(depth / mt))
    gets the missing turns of its loop inserted (a loop through segments that all miss copies). Returns
    (circle, note)."""
    notes = []
    for _ in range(20):
        used = {}
        for s, _o in circ:
            used[s] = used.get(s, 0) + 1

        def miss(s, used=used):
            return copies(s) - used.get(s, 0)

        for i, x in enumerate(circ):
            if miss(x[0]) > 0:
                loops = walks(adj, lambda s, x=x, miss=miss: miss(s) > 0 and s != x[0], x, x, max_steps=10)
                if loops:
                    loop = min(loops, key=len)
                    n = min(miss(s) for s, _ in loop + [x])
                    circ = circ[:i + 1] + (loop + [x]) * n + circ[i + 1:]
                    notes.append(f"repeat unit {'-'.join(s for s, _ in [x] + loop)} +{n} copies (depth)")
                    break
        else:
            break
    return circ, "; ".join(notes)


def path_sequence(circ, seq, adj, circular=True):
    """Sequence of a segment path; circular: each segment minus its overlap with the one before (wrapping)."""
    def part(s, o):
        return seq[s] if o == "+" else rc(seq[s])
    ov = {(a, ao, b, bo): o for a, ao in set(circ) for b, bo, o in adj.get((a, ao), [])}
    out = [] if circular else [part(*circ[0])]
    for k in range(0 if circular else 1, len(circ)):
        (a, ao), (b, bo) = circ[k - 1], circ[k]
        out.append(part(b, bo)[ov[(a, ao, b, bo)]:])
    return "".join(out)

def graph_check(path, links, length):
    """(closes, overlap, repeat note, branch note) for one contig path."""
    flat = [x for chunk in path for x in chunk]
    first, last = flat[0], flat[-1]
    ov = links.get((last[0], last[1], first[0], first[1]))
    counts = {}
    for s, _ in flat:
        counts[s] = counts.get(s, 0) + 1
    reps = [f"{s} ({length.get(s, '?')} bp) x{n}" for s, n in counts.items() if n > 1]
    outs = {(a, ao) for (a, ao, b, bo) in links if (a, ao) == last and (b, bo) != first}
    branch = f"end has {len(outs) + (ov is not None)} continuations" if outs else ""
    return ov is not None, ov or 0, "; ".join(reps), branch


def overlap_len(a, b, lo=20, hi=300):
    """Longest exact suffix(a) == prefix(b) of lo..hi bp, or 0."""
    for n in range(min(hi, len(a), len(b)), lo - 1, -1):
        if a[-n:] == b[:n]:
            return n
    return 0


def try_join(contigs, links=None, paths=None):
    """Order / orient up to 6 contigs into a closed chain. Returns (circle, description) or None."""
    if not 1 < len(contigs) <= 6:
        return None
    names = list(contigs)
    for perm in itertools.permutations(names[1:]):
        for orients in itertools.product("+-", repeat=len(perm)):
            chain = [(names[0], "+")] + list(zip(perm, orients))
            seqs = [contigs[n] if o == "+" else rc(contigs[n]) for n, o in chain]
            ovs = []
            for i in range(len(chain)):
                (a, ao), (b, bo) = chain[i], chain[(i + 1) % len(chain)]
                ov = None
                if links is not None and a in paths and b in paths:     # graph link between the contigs' ends
                    fa = [x for c in paths[a] for x in c]
                    fb = [x for c in paths[b] for x in c]
                    end = fa[-1] if ao == "+" else (fa[0][0], FLIP[fa[0][1]])
                    start = fb[0] if bo == "+" else (fb[-1][0], FLIP[fb[-1][1]])
                    ov = links.get((end[0], end[1], start[0], start[1]))
                if ov is None:
                    ov = overlap_len(seqs[i], seqs[(i + 1) % len(seqs)]) or None
                if ov is None:
                    break
                ovs.append(ov)
            else:
                circle = "".join(s[:len(s) - ov] for s, ov in zip(seqs, ovs))
                return circle, " > ".join(f"{n}{o}" for n, o in chain) + f" (overlaps {ovs})"
    return None


def circularize(recs, spades=None):
    """recs: [(name, seq)]. Returns (rows, [(name, circle or linear seq)])."""
    # MitoFinder writes some contigs reverse-complemented ("NODE_x (reverse)"): back to the SPAdes orientation
    seqs = {n.split()[0]: (rc(s) if "(reverse)" in n else s) for n, s in recs}
    paths, links, length = {}, None, {}
    if spades and Path(spades, "scaffolds.paths").exists() and Path(spades, "assembly_graph_with_scaffolds.gfa").exists():
        paths = read_paths(Path(spades, "scaffolds.paths"), set(seqs))
        segs = {s for p in paths.values() for c in p for s, _ in c}
        length, links = read_links(Path(spades, "assembly_graph_with_scaffolds.gfa"), segs)
    rows, out, unclosed = [], [], {}
    for cid, s in seqs.items():
        row = dict(contig=cid, length=len(s), method="not_circular", circle_length="", overlap_bp="",
                   gaps=s.upper().count("N"), repeats="", branches="", note="")
        circle = None
        if cid in paths:
            closes, ov, reps, branch = graph_check(paths[cid], links, length)
            row.update(repeats=reps, branches=branch)
            if len(paths[cid]) > 1:
                row["note"] = f"{len(paths[cid]) - 1} scaffold gap(s) in the path"
            if closes and (ov == 0 or s[-ov:] == s[:ov]):
                circle = s[:len(s) - ov]
                row.update(method="graph", overlap_bp=ov)
        if circle is None and cid in paths:
            unclosed[cid] = len(rows)
        if circle is not None:
            row["circle_length"] = len(circle)
        rows.append(row)
        out.append((cid, circle or s))
    if unclosed:                                 # graph walk from each open end back to the start
        gfa = Path(spades, "assembly_graph_with_scaffolds.gfa")
        adj, depth = read_graph(gfa)
        for cid, i in unclosed.items():
            if len(paths[cid]) > 1:              # chunks separated by a scaffold gap: no link to walk from
                rows[i]["note"] = "; ".join(x for x in (rows[i]["note"], "scaffold gap: no graph walk") if x)
                continue
            flat = [x for c in paths[cid] for x in c]
            mt = depth[max(flat, key=lambda x: length.get(x[0], 0))[0]]   # depth of the contig's longest segment
            ws = walks(adj, lambda s, mt=mt: depth.get(s, 0) >= 0.2 * mt, flat[-1], flat[0])
            if not ws:
                rows[i]["note"] = "; ".join(x for x in (rows[i]["note"], "no mt-depth graph walk closes the ends") if x)
                continue
            best = max(ws, key=lambda w: (min((depth[s] for s, _ in w), default=float("inf")), -len(w)))
            circ, rep = add_repeat_copies(flat + best, adj, lambda s, mt=mt: round(depth.get(s, 0) / mt))
            seq = read_sequences(gfa, {s for s, _ in circ})
            if seqs[cid] not in path_sequence(flat, seq, adj, circular=False):   # MitoFinder may trim the scaffold
                rows[i]["note"] = "; ".join(x for x in (rows[i]["note"], "graph path does not rebuild the contig") if x)
                continue
            circle = path_sequence(circ, seq, adj)
            used = {}
            for x, _ in circ:
                used[x] = used.get(x, 0) + 1
            odd = [f"{x} x{n} in circle, depth ~x{depth[x] / mt:.1f}" for x, n in used.items()
                   if max(n, round(depth[x] / mt)) > 1 and n != round(depth[x] / mt)]
            if odd:                              # short reads cannot count tandem copies exactly: say so
                rep = "; ".join(x for x in (rep, "repeat copies uncertain: " + ", ".join(odd)) if x)
            filled = ", ".join(f"{s}{o} ({len(seq[s])} bp, {depth[s]:.0f}x)" for s, o in best)
            rows[i].update(method="walk", circle_length=len(circle), overlap_bp="",
                           note="; ".join(x for x in (f"gap filled from the graph: {filled or 'direct link'} "
                                                      f"(mt depth {mt:.0f}x, {len(ws)} walk(s))", rep) if x))
            out[i] = (cid, circle)
        del adj, depth

    for i, row in enumerate(rows):               # no graph (or no walk): overlapping ends
        if row["method"] == "not_circular":
            s = seqs[row["contig"]]
            c = find_circle(s)
            if c:
                circle = s[c["anchor"]:c["anchor"] + c["P"]]
                row.update(method="overlap", overlap_bp=c["overlap"], circle_length=len(circle))
                out[i] = (row["contig"], circle)
    if len(recs) > 1 and not any(r["method"] != "not_circular" for r in rows):
        joined = try_join(seqs, links, paths)
        if joined:
            circle, desc = joined
            rows.append(dict(contig="+".join(seqs), length=sum(map(len, seqs.values())), method="join",
                             circle_length=len(circle), overlap_bp="", gaps=circle.count("N"), repeats="",
                             branches="", note=desc))
            out = [("joined", circle)]
    return rows, out


FIELDS = ["sample", "contig", "length", "method", "circle_length", "overlap_bp", "gaps", "repeats", "branches", "note"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="mitochondrial contigs (FASTA), e.g. MitoFinder *_all_mt.fasta")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--spades", help="SPAdes output directory (scaffolds.paths + assembly_graph_with_scaffolds.gfa)")
    ap.add_argument("--prefix", default="sample")
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rows, seqs = circularize(read_fasta(a.input), a.spades)
    circ = {r["contig"] for r in rows if r["method"] != "not_circular"}
    with open(out / f"{a.prefix}.circularized.fasta", "w") as fh:
        write_fasta([(f"{a.prefix}|{n} {'circular' if n in circ or n == 'joined' else 'linear'} length={len(s)}", s)
                     for n, s in seqs], fh)
    with open(out / f"{a.prefix}.circularize.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows({"sample": a.prefix, **r} for r in rows)
    for r in rows:
        print("\t".join(str(r[k]) for k in FIELDS[1:]))


if __name__ == "__main__":
    main()

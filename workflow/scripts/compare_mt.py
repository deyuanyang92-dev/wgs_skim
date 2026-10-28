"""Compare mitogenome assemblies of one sample: GetOrganelle (reference) against SPAdes-MitoFinder or any other.

Why: different assemblers should give the same mitogenome; where they do not, it is almost always a region short
reads cannot resolve - a tandem repeat in the control region, a multi-copy segment, a gap. GetOrganelle's assembly
graph shows those regions (segments used several times in the path, segments with several times the depth), so its
result is the reference and every difference is reported with the graph region and gene it falls in.

Steps: the other assembly is circularized first if needed (circularize.py: graph / walk / overlap / join, with
--spades), turned to the reference strand and rotated to its start, then aligned to it (MAFFT, global). Differences
are merged into blocks (mismatch, insertion, deletion) and placed on the reference: gene (--gb) and GetOrganelle
problem region (repeat segment, high depth, Ns).

Verdict: identical | consistent (every difference is inside a flagged repeat / gap region, at a linear end, or a
<= 2-bp length difference of a homopolymer >= 6 bp - all things short reads do not settle) | differ.

Usage: compare_mt.py -a getorganelle_out_dir_or_fasta -b mitofinder_mt.fasta [-b label=other.fasta]
                     [--spades spades_dir] [--gb reference.gb] [--gfa reference_graph.gfa] -o out/ --prefix sample
       -> out/<prefix>.compare.tsv (one row per comparison), <prefix>.diffs.tsv, <prefix>.<label>.aligned.fasta
Needs mafft on PATH.
"""
import argparse
import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from circularize import circularize, rc  # noqa: E402
from seqtools import read_fasta, write_fasta  # noqa: E402


def getorganelle_result(path):
    """(fasta, gfa or None, number of alternative paths) from a GetOrganelle output directory or a FASTA."""
    p = Path(path)
    if p.is_file():
        return p, None, 1
    for kind in ("complete", "scaffolds", "contigs"):
        fs = sorted(p.glob(f"*.{kind}.graph1.*.path_sequence.fasta"))
        if fs:
            gfa = next(iter(p.glob(f"*.{kind}.graph1.selected_graph.gfa")), None)
            return fs[0], gfa, len(fs)
    sys.exit(f"no GetOrganelle path_sequence.fasta in {p}")


def graph_regions(header, gfa, ref):
    """Problem regions of the reference from its GetOrganelle path (FASTA header '163-,1174+,...(circular)') and
    graph: [(start, end, label)], 0-based on ref. Segments used more than once in the path or with >= 1.5x the depth
    of the longest segment are repeats short reads cannot count."""
    path = [(x[:-1], x[-1]) for x in re.sub(r"\(.*", "", header.split()[0]).split(",") if x and x[-1] in "+-"]
    if not gfa or not path:
        return []
    seq, depth, ov = {}, {}, {}
    for line in open(gfa):
        t = line.rstrip("\n").split("\t")
        if t[0] == "S":
            tags = dict(x.split(":", 2)[::2] for x in t[3:] if x.count(":") >= 2)
            seq[t[1]] = t[2]
            depth[t[1]] = float(tags["DP"]) if "DP" in tags else float(tags.get("RC", 0)) / max(len(t[2]), 1)
        elif t[0] == "L":
            o = int(re.match(r"\d+", t[5]).group()) if t[5][0].isdigit() else 0
            ov[(t[1], t[2], t[3], t[4])] = ov[(t[3], "+-"[t[4] == "+"], t[1], "+-"[t[2] == "+"])] = o
    if any(s not in seq for s, _ in path):
        return []
    spans, built = [], ""
    for k, (s, o) in enumerate(path):               # circular: each segment minus its overlap with the previous
        prev = path[k - 1]
        cut = ov.get((prev[0], prev[1], s, o), 0) if len(path) > 1 else 0
        part = seq[s] if o == "+" else rc(seq[s])
        spans.append((len(built) - cut, len(built) + len(part) - cut, s))
        built += part[cut:]
    shift = (built + built).find(ref[:60])
    if shift < 0:
        return []
    main = depth[max(seq, key=lambda s: len(seq[s]))]
    used = {s: sum(x == s for x, _ in path) for s, _ in path}
    n = len(built)
    out = []
    for a, b, s in spans:
        why = []
        if used[s] > 1:
            why.append(f"x{used[s]} in path")
        if depth[s] >= 1.5 * main:
            why.append(f"depth {depth[s] / main:.1f}x")
        if why:
            out.append(((a - shift) % n, (a - shift) % n + (b - a), f"repeat segment {s} ({len(seq[s])} bp, {', '.join(why)})"))
    return out


def genes(gb, ref):
    """[(start, end, name)] 0-based on ref from a GenBank file (gene / CDS / tRNA / rRNA / D-loop / misc_feature);
    the annotated sequence may start elsewhere on the circle (shifted) but must be on the same strand."""
    out, cur, origin = [], None, []
    lines = open(gb).read().split("\n")
    if "ORIGIN" in [x.strip() for x in lines]:
        origin = lines[[x.strip() for x in lines].index("ORIGIN") + 1:]
    gseq = re.sub(r"[^a-z]", "", "".join(origin).lower()).upper()
    shift = (gseq + gseq).find(ref[:60]) if gseq else 0
    if shift < 0:
        print(f"{gb}: annotated sequence is not the reference (or other strand): genes not used", file=sys.stderr)
        return []
    for line in lines:
        m = re.match(r" {5}(\S+)\s+(\S+)", line)
        if m:
            cur = None
            if m.group(1) in ("gene", "CDS", "tRNA", "rRNA", "D-loop", "misc_feature", "rep_origin"):
                pos = [int(x) for x in re.findall(r"\d+", m.group(2))]
                if pos:
                    cur = [min(pos) - 1, max(pos), m.group(1)]
                    out.append(cur)
        elif cur and re.search(r'/(gene|product|note)="', line) and cur[2] in ("gene", "CDS", "tRNA", "rRNA", "D-loop",
                                                                           "misc_feature", "rep_origin"):
            cur[2] = re.search(r'="([^"]*)', line).group(1)
    n = len(gseq) or len(ref)
    return [((a - shift) % n, (a - shift) % n + b - a, name) for a, b, name in out]


def orient(ref, seq, circular=True, k=31):
    """seq on the reference strand, rotated to the reference start when circular: (seq, strand, rotation)."""
    for strand, s in (("+", seq), ("-", rc(seq))):
        for i in range(0, min(len(ref) - k, 3000), 100):          # first unique k-mer of ref found in seq
            j = (s + s).find(ref[i:i + k])
            if j >= 0:
                r = (j - i) % len(s) if circular else 0
                return s[r:] + s[:r], strand, r
    return seq, "?", 0


def mafft(a, b):
    with tempfile.NamedTemporaryFile("w", suffix=".fa") as fh:
        fh.write(f">ref\n{a}\n>other\n{b}\n")
        fh.flush()
        res = subprocess.run(["mafft", "--auto", "--quiet", "--thread", "4", fh.name], capture_output=True, text=True,
                             check=True)
    return [s for _, s in read_fasta_text(res.stdout)]


def read_fasta_text(text):
    recs = []
    for line in text.splitlines():
        if line.startswith(">"):
            recs.append([line[1:], ""])
        elif recs:
            recs[-1][1] += line.strip().upper()
    return recs


def diff_blocks(a, b, merge=10):
    """Differences of an alignment as blocks on a: [(start, end, a_bp, b_bp, a_seq, b_seq)] (0-based, a coords)."""
    blocks, pos, cur = [], 0, None
    for x, y in zip(a, b):
        if x != y:
            if cur and pos - cur["last"] <= merge:
                cur.update(last=pos)
            else:
                cur = {"start": pos, "last": pos, "a": "", "b": ""}
                blocks.append(cur)
            cur["a"] += x if x != "-" else ""
            cur["b"] += y if y != "-" else ""
        elif cur and pos - cur["last"] <= merge:
            cur["a"] += x
            cur["b"] += y
        if x != "-":
            pos += 1
    out = []
    for c in blocks:                                # trim the identical tail picked up while merging
        a_s, b_s = c["a"], c["b"]
        while a_s and b_s and a_s[-1] == b_s[-1]:
            a_s, b_s = a_s[:-1], b_s[:-1]
        out.append((c["start"], c["start"] + len(a_s), a_s, b_s))
    return out


def kind(a_s, b_s):
    if not a_s:
        return "insertion"
    if not b_s:
        return "deletion"
    return "mismatch" if len(a_s) == len(b_s) and len(a_s) <= 3 else "complex"


def compare(ref, ref_circular, regions, feats, label, other, other_circular, method):
    seq, strand, rot = orient(ref, other, other_circular)
    if strand == "?":
        return dict(other=label, verdict="differ", note="no shared 31-mer with the reference"), [], None
    a, b = mafft(ref, seq)
    blocks = diff_blocks(a, b)
    n = len(ref)
    rows, explained = [], 0
    for s, e, a_s, b_s in blocks:
        hit = [lab for rs, re_, lab in regions if s < re_ and e + 1 > rs or s < re_ - n and e + 1 > rs - n]
        if "N" in a_s or "N" in b_s:
            hit.append("Ns (gap)")
        indel = a_s + b_s
        if indel and len(set(indel)) == 1 and abs(len(a_s) - len(b_s)) <= 2:
            i = j = min(s, n - 1)                     # the run of that base around the difference
            while i > 0 and ref[i - 1] == indel[0]:
                i -= 1
            while j < n and ref[j] == indel[0]:
                j += 1
            run = j - i
            if run >= 6:
                hit.append(f"homopolymer {indel[0]}x{run}")
        end = not ref_circular and (s == 0 or e >= n)
        if hit or end:
            explained += 1
        g = sorted({name for gs, ge, name in feats if gs <= e and ge >= s})
        if feats and not g:                          # intergenic: name the flanking genes
            before = max((f for f in feats if f[1] <= s), key=lambda f: f[1], default=None)
            after = min((f for f in feats if f[0] > e), key=lambda f: f[0], default=None)
            g = [f"intergenic {before[2] if before else 'start'}..{after[2] if after else 'end'}"]
        rows.append({"other": label, "ref_start": s + 1, "ref_end": max(e, s + 1), "type": kind(a_s, b_s),
                     "ref_bp": len(a_s), "other_bp": len(b_s), "ref_seq": a_s[:40], "other_seq": b_s[:40],
                     "gene": ", ".join(g) or "intergenic", "graph_region": "; ".join(hit) or ("linear end" if end else "")})
    mis = sum(r["ref_bp"] for r in rows if r["type"] == "mismatch")
    indel = sum(abs(r["ref_bp"] - r["other_bp"]) for r in rows)
    same = sum(x == y for x, y in zip(a, b))
    verdict = "identical" if not rows else "consistent" if explained == len(rows) else "differ"
    summary = {"other": label, "other_method": method, "ref_length": n, "other_length": len(other),
               "ref_topology": "circular" if ref_circular else "linear",
               "other_topology": "circular" if other_circular else "linear", "strand": strand, "rotation": rot,
               "identity_pct": f"{100 * same / len(a):.3f}", "diff_blocks": len(rows), "mismatch_bp": mis,
               "indel_bp": indel, "unexplained_blocks": len(rows) - explained, "verdict": verdict}
    return summary, rows, (a, b)


SUMMARY = ["sample", "other", "other_method", "ref_length", "other_length", "ref_topology", "other_topology",
           "strand", "rotation", "identity_pct", "diff_blocks", "mismatch_bp", "indel_bp", "unexplained_blocks",
           "verdict", "ref_alternative_paths", "ref_problem_regions", "note"]
DIFFS = ["sample", "other", "ref_start", "ref_end", "type", "ref_bp", "other_bp", "gene", "graph_region", "ref_seq",
         "other_seq"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-a", "--ref", required=True, help="GetOrganelle output directory or FASTA (reference)")
    ap.add_argument("-b", "--other", action="append", required=True, help="[label=]FASTA to compare (repeatable)")
    ap.add_argument("--spades", help="SPAdes directory of the other assembly (graph-based circularization)")
    ap.add_argument("--gb", help="GenBank annotation of the reference (genes for each difference)")
    ap.add_argument("--gfa", help="reference assembly graph, if -a is a FASTA")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--prefix", default="sample")
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    fasta, gfa, alts = getorganelle_result(a.ref)
    (header, ref), = read_fasta(fasta)[:1]
    ref_circ = "circular" in header and "linear" not in header
    regions = graph_regions(header, a.gfa or gfa, ref)
    feats = genes(a.gb, ref) if a.gb else []
    summaries, diffs = [], []
    for item in a.other:
        label, path = item.split("=", 1) if "=" in item and not Path(item).exists() else (Path(item).stem, item)
        recs = read_fasta(path)
        if any("circular" in n and "linear" not in n for n, _ in recs) and len(recs) == 1:
            seq, circ, method = recs[0][1], True, "as given (circular)"
        else:
            rows, seqs = circularize(recs, a.spades)
            done = [(r, s) for r, (_, s) in zip(rows, seqs) if r["method"] != "not_circular"]
            if rows and rows[-1]["method"] == "join":
                done = [(rows[-1], seqs[0][1])]
            if done:
                r, seq = max(done, key=lambda x: len(x[1]))
                circ, method = True, f"circularized: {r['method']} ({r['note'] or r['contig']})"[:300]
            else:
                seq, circ, method = max((s for _, s in recs), key=len), False, "not circular (longest contig)"
        with open(out / f"{a.prefix}.{label}.oriented.fasta", "w") as fh:
            write_fasta([(f"{a.prefix}|{label} {'circular' if circ else 'linear'}", orient(ref, seq, circ)[0])], fh)
        s, rows, aln = compare(ref, ref_circ, regions, feats, label, seq, circ, method)
        s.update(sample=a.prefix, ref_alternative_paths=alts,
                 ref_problem_regions="; ".join(f"{x + 1}-{y}: {lab}" for x, y, lab in regions))
        summaries.append(s)
        diffs += [{"sample": a.prefix, **r} for r in rows]
        if aln:
            with open(out / f"{a.prefix}.{label}.aligned.fasta", "w") as fh:
                write_fasta([(f"{a.prefix}|getorganelle", aln[0]), (f"{a.prefix}|{label}", aln[1])], fh)
        print(f"{a.prefix} getorganelle {len(ref)} bp vs {label} {len(seq)} bp: {s['verdict']} "
              f"(identity {s.get('identity_pct')}%, {s.get('diff_blocks')} difference block(s), "
              f"{s.get('unexplained_blocks')} outside GetOrganelle-flagged regions)")
        for r in rows:
            print(f"  {r['ref_start']}-{r['ref_end']} {r['type']} {r['ref_bp']}->{r['other_bp']} bp  {r['gene']}  "
                  f"{r['graph_region'] or '-'}")
    for name, fields, data in (("compare", SUMMARY, summaries), ("diffs", DIFFS, diffs)):
        with open(out / f"{a.prefix}.{name}.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore", restval="")
            w.writeheader()
            w.writerows(data)


if __name__ == "__main__":
    main()

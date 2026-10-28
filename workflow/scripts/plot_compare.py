"""Figure of one assembly comparison as SVG: GetOrganelle on top, the other assembly below, genes on both,
alignment links between them; differences in red, GetOrganelle repeat regions (from its graph) in purple.

Why: a table of differences says where; the figure shows at a glance whether they sit in the control region, a gene
or a repeat the graph already flagged. Genes come from the reference annotation and are carried to the other assembly
through the alignment, so one annotation is enough. Plain SVG from the standard library (pyGenomeViz drew the same
figure but needs matplotlib / numpy / Biopython, a 389 MB environment): opens in any browser with hover details
(<title>), and is vector for papers (Inkscape / Illustrator).

Used by compare_mt.py (draw); standalone:
plot_compare.py --aligned P.L.aligned.fasta --genes P.genes.tsv --regions P.regions.tsv --diffs P.diffs.tsv
                --label L -o P.L.compare.svg
"""
import argparse
import csv
from html import escape

COLORS = {"CDS": "#e8a33d", "rRNA": "#3b82c4", "tRNA": "#4caf50", "other": "#9e9e9e"}
DIFF, REPEAT, LINK = "#d62728", "#8e44ad", "#c8c8c8"


def gene_type(name, typ):
    """tRNA / rRNA / CDS from the feature type, or from the name when every feature is a 'gene' (MitoZ, Mitom)."""
    n = name.lower()
    if typ in ("tRNA", "rRNA", "CDS"):
        return typ
    if n.startswith(("trn", "trna")):
        return "tRNA"
    if n.startswith(("rrn", "rns", "rnl", "12s", "16s", "s-rrna", "l-rrna")):
        return "rRNA"
    return "CDS" if typ == "gene" else "other"


def column_map(a, b):
    """before[i] = bases of b before the alignment column of a's base i; before[len(a)] = len(b)."""
    before, j = [], 0
    for x, y in zip(a, b):
        if x != "-":
            before.append(j)
        j += y != "-"
    return before + [j]


def matched_blocks(a, b, min_len=20):
    """Gap-free stretches of the alignment: (a_start, a_end, b_start, b_end)."""
    out, i, j, cur = [], 0, 0, None
    for x, y in zip(a + "-", b + "-"):
        if x != "-" and y != "-":
            cur = cur or [i, j, 0]
            cur[2] += 1
        elif cur:
            if cur[2] >= min_len:
                out.append((cur[0], cur[0] + cur[2], cur[1], cur[1] + cur[2]))
            cur = None
        i += x != "-"
        j += y != "-"
    return out


def draw(aligned, genes, regions, diffs, label, title=""):
    """SVG text. aligned = (ref with gaps, other with gaps); genes = [(start, end, name, strand, type)] 0-based on ref;
    regions = [(start, end, label)]; diffs = rows of compare_mt.compare (ref_start 1-based)."""
    ra, ob = aligned
    rlen, olen = len(ra.replace("-", "")), len(ob.replace("-", ""))
    before = column_map(ra, ob)
    W, left, right = 1500, 190, 30
    k = (W - left - right) / max(rlen, olen)
    X = lambda p: left + p * k  # noqa: E731
    y1, y2, h = 110, 250, 9                          # track baselines, half height of a big arrow
    el = []

    def arrow(s, e, y, strand, color, tip, half):
        x1, x2 = X(s), X(e)
        head = min(8, (x2 - x1) * 0.6)
        if strand < 0:
            pts = [(x2, y - half), (x1 + head, y - half), (x1, y), (x1 + head, y + half), (x2, y + half)]
        else:
            pts = [(x1, y - half), (x2 - head, y - half), (x2, y), (x2 - head, y + half), (x1, y + half)]
        el.append(f'<polygon points="{" ".join(f"{px:.1f},{py:.1f}" for px, py in pts)}" fill="{color}" '
                  f'stroke="#333" stroke-width="0.4"><title>{escape(tip)}</title></polygon>')

    def track(y, length, name, feats, label_above):
        el.append(f'<text x="{left - 12}" y="{y + 5}" text-anchor="end" font-size="15">{escape(name)}</text>')
        el.append(f'<line x1="{X(0)}" y1="{y}" x2="{X(length)}" y2="{y}" stroke="#555" stroke-width="1"/>')
        for s, e, gname, strand, typ in feats:
            s, e = max(0, s), min(length, e)
            if e <= s:
                continue
            t = gene_type(gname, typ)
            big = t in ("CDS", "rRNA")
            arrow(s, e, y, strand, COLORS[t], f"{gname} ({t}) {s + 1}-{e}", h if big else h * 0.7)
            if big:
                ty = y - h - 6 if label_above else y + h + 14
                el.append(f'<text x="{X((s + e) / 2):.1f}" y="{ty}" text-anchor="middle" font-size="11">'
                          f'{escape(gname)}</text>')

    for s, e, os_, oe in matched_blocks(ra, ob):     # alignment links first, under everything else
        el.append(f'<polygon points="{X(s):.1f},{y1 + h} {X(e):.1f},{y1 + h} {X(oe):.1f},{y2 - h} {X(os_):.1f},'
                  f'{y2 - h}" fill="{LINK}" opacity="0.6"><title>aligned {s + 1}-{e} / {os_ + 1}-{oe}</title>'
                  '</polygon>')
    for s, e, lab in regions:                        # GetOrganelle graph: repeats short reads cannot count
        s, e = max(0, s), min(rlen, e)
        el.append(f'<rect x="{X(s):.1f}" y="{y1 - h - 4}" width="{max(X(e) - X(s), 2):.1f}" height="{2 * h + 8}" '
                  f'fill="{REPEAT}" opacity="0.35"><title>{escape(f"{s + 1}-{e}: {lab}")}</title></rect>')
    other_genes = [(before[s], before[max(e - 1, s)] + 1, n, st, t) for s, e, n, st, t in genes if s < rlen]
    track(y1, rlen, f"GetOrganelle ({rlen:,} bp)", genes, True)
    track(y2, olen, f"{label} ({olen:,} bp)", other_genes, False)
    for d in diffs:
        s, e, obp = int(d["ref_start"]) - 1, int(d["ref_end"]), int(d["other_bp"])
        if d["type"] == "insertion":
            e = s
            os_, oe = before[s] - obp, before[s]
        else:
            os_ = before[s]
            oe = os_ + obp
        tip = (f"{d['type']} {d['ref_bp']} -> {d['other_bp']} bp at {d['ref_start']}-{d['ref_end']} | {d['gene']} | "
               f"{d['graph_region'] or 'outside flagged regions'}")
        xs, xe, xos, xoe = X(s), max(X(e), X(s) + 2), X(os_), max(X(oe), X(os_) + 2)
        el.append(f'<polygon points="{xs:.1f},{y1 + h} {xe:.1f},{y1 + h} {xoe:.1f},{y2 - h} {xos:.1f},{y2 - h}" '
                  f'fill="{DIFF}" opacity="0.8"><title>{escape(tip)}</title></polygon>')
        for y, a_, b_ in ((y1, xs, xe), (y2, xos, xoe)):
            el.append(f'<rect x="{a_:.1f}" y="{y - h - 3}" width="{b_ - a_:.1f}" height="{2 * h + 6}" fill="{DIFF}">'
                      f'<title>{escape(tip)}</title></rect>')
    ya = y2 + 45                                     # scale
    step = 2000 if max(rlen, olen) > 8000 else 1000
    el.append(f'<line x1="{X(0)}" y1="{ya}" x2="{X(max(rlen, olen))}" y2="{ya}" stroke="#555"/>')
    for p in range(0, max(rlen, olen) + 1, step):
        el.append(f'<line x1="{X(p):.1f}" y1="{ya}" x2="{X(p):.1f}" y2="{ya + 5}" stroke="#555"/>'
                  f'<text x="{X(p):.1f}" y="{ya + 18}" text-anchor="middle" font-size="11">{p // 1000} kb</text>')
    legend = [(DIFF, "difference"), (REPEAT, "GetOrganelle repeat region (graph)"), (COLORS["CDS"], "CDS"),
              (COLORS["rRNA"], "rRNA"), (COLORS["tRNA"], "tRNA")]
    x = left
    for color, text in legend:
        el.append(f'<rect x="{x}" y="{ya + 32}" width="12" height="12" fill="{color}"/>'
                  f'<text x="{x + 17}" y="{ya + 43}" font-size="12">{text}</text>')
        x += 30 + 7 * len(text)
    head = f'<text x="{left}" y="28" font-size="16" font-weight="bold">{escape(title)}</text>' if title else ""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{ya + 60}" viewBox="0 0 {W} {ya + 60}" '
            f'font-family="Arial, Helvetica, sans-serif"><rect width="100%" height="100%" fill="white"/>{head}'
            + "".join(el) + "</svg>\n")


def read_fasta(path):
    recs = []
    for line in open(path):
        if line.startswith(">"):
            recs.append("")
        elif recs:
            recs[-1] += line.strip()
    return recs


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aligned", required=True)
    ap.add_argument("--genes")
    ap.add_argument("--regions")
    ap.add_argument("--diffs", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("-o", "--out", required=True, help="output .svg")
    a = ap.parse_args()
    genes = [(int(r["start"]), int(r["end"]), r["name"], int(r["strand"]), r["type"])
             for r in tsv(a.genes)] if a.genes else []
    regions = [(int(r["start"]), int(r["end"]), r["label"]) for r in tsv(a.regions)] if a.regions else []
    diffs = [r for r in tsv(a.diffs) if r["other"] == a.label]
    with open(a.out, "w") as fh:
        fh.write(draw(read_fasta(a.aligned)[:2], genes, regions, diffs, a.label, a.title))


if __name__ == "__main__":
    main()

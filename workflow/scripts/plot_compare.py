"""Figure of one assembly comparison as SVG: GetOrganelle on top, the other assembly below (its contigs ordered
along the reference and turned to its strand), genes on both, alignments as links between them (grey; orange when
inverted); differences red, GetOrganelle repeat regions (from its graph) purple, reference parts missing from the
assembly light red, misassembly junctions marked with a red triangle.

Why: the tables say where; the figure shows at a glance whether problems sit in the control region, a gene or a
repeat the graph already flagged. Genes come from the reference annotation and are carried to the other assembly
through the alignments, so one annotation is enough. Plain SVG from the standard library (pyGenomeViz draws similar
figures but needs matplotlib / numpy / Biopython, a 389 MB environment): opens in any browser with hover details
(<title>), and is vector for papers (Inkscape / Illustrator). Used by compare_mt.py.
"""
from html import escape

from assembly_check import gene_type

COLORS = {"CDS": "#e8a33d", "rRNA": "#3b82c4", "tRNA": "#4caf50", "other": "#9e9e9e"}
DIFF, REPEAT, LINK, INVERTED, MISSING = "#d62728", "#8e44ad", "#c8c8c8", "#f39c12", "#f4b6b6"


def layout(contigs, blocks):
    """Contigs in reference order, each turned to the reference strand: {name: (offset, length, flipped)}, total."""
    info = {}
    for name, seq, _ in contigs:
        bs = [b for b in blocks if b["contig"] == name]
        big = max(bs, key=lambda b: b["re"] - b["rs"], default=None)
        info[name] = (big["rs"] if big else 10 ** 9, len(seq), bool(big) and big["strand"] == "-")
    gap = max(1, int(0.008 * sum(v[1] for v in info.values()))) if len(info) > 1 else 0
    pos, off = {}, 0
    for name in sorted(info, key=lambda c: info[c][0]):
        pos[name] = (off, info[name][1], info[name][2])
        off += info[name][1] + gap
    return pos, off - gap


def draw(rlen, contigs, blocks, genes, regions, diffs, events, label, title=""):
    """SVG text. contigs = [(name, seq, circular)]; blocks = assembly_check.minimap2 alignments; genes =
    [(start, end, name, strand, type)] 0-based on the reference; regions = [(start, end, label)]; diffs = rows of
    compare_mt.describe; events = assembly_check.structure events."""
    pos, olen = layout(contigs, blocks)
    W, left, right = 1560, 260, 30
    k = (W - left - right) / max(rlen, olen)
    X = lambda p: left + p * k  # noqa: E731
    y1, y2, h = 110, 250, 9                          # track baselines, half height of a big arrow
    el = []

    def qx(contig, q):                               # contig position -> x on the lower track
        off, ln, flip = pos[contig]
        return X(off + (ln - q if flip else q))

    def pieces(rs, re_):                             # a reference interval split at the origin of the circle
        return [(rs, min(re_, rlen), 0)] + ([(rlen, re_, rlen)] if re_ > rlen else [])

    def to_query(b, p):
        return b["qs"] + (p - b["rs"]) if b["strand"] == "+" else b["qe"] - (p - b["rs"])

    def arrow(x1, x2, y, strand, color, tip, half):
        x1, x2 = min(x1, x2), max(x1, x2)
        head = min(8, (x2 - x1) * 0.6)
        if strand < 0:
            pts = [(x2, y - half), (x1 + head, y - half), (x1, y), (x1 + head, y + half), (x2, y + half)]
        else:
            pts = [(x1, y - half), (x2 - head, y - half), (x2, y), (x2 - head, y + half), (x1, y + half)]
        el.append(f'<polygon points="{" ".join(f"{px:.1f},{py:.1f}" for px, py in pts)}" fill="{color}" '
                  f'stroke="#333" stroke-width="0.4"><title>{escape(tip)}</title></polygon>')

    def gene(x1, x2, y, strand, gname, typ, above, tip):
        t = gene_type(gname, typ)
        big = t in ("CDS", "rRNA")
        arrow(x1, x2, y, strand, COLORS[t], tip, h if big else h * 0.7)
        if big and abs(x2 - x1) > 12:
            ty = y - h - 6 if above else y + h + 14
            el.append(f'<text x="{(x1 + x2) / 2:.1f}" y="{ty}" text-anchor="middle" font-size="11">'
                      f'{escape(gname)}</text>')

    for b in blocks:                                 # links first, under everything else
        flip = pos[b["contig"]][2]
        inv = (b["strand"] == "-") != flip
        for rs, re_, sh in pieces(b["rs"], b["re"]):     # rs, re_ in alignment coordinates; drawn at - sh
            a1, a2 = qx(b["contig"], to_query(b, rs)), qx(b["contig"], to_query(b, re_))
            el.append(f'<polygon points="{X(rs - sh):.1f},{y1 + h} {X(re_ - sh):.1f},{y1 + h} {a2:.1f},{y2 - h} '
                      f'{a1:.1f},{y2 - h}" fill="{INVERTED if inv else LINK}" opacity="0.6"><title>'
                      f'{escape(b["contig"])} {b["qs"] + 1}-{b["qe"]} ({b["strand"]}) = reference '
                      f'{b["rs"] % rlen + 1}-{(b["re"] - 1) % rlen + 1}</title></polygon>')
    for e in events:                                 # reference parts the assembly lacks
        if e["event"] == "not in assembly":
            s, t = (int(x) for x in e["ref_pos"].split("-"))
            el.append(f'<rect x="{X(s - 1):.1f}" y="{y1 - h - 6}" width="{max(X(t) - X(s - 1), 2):.1f}" '
                      f'height="{2 * h + 12}" fill="{MISSING}"><title>not in {escape(label)}: {s}-{t}'
                      f'{" | " + escape(e["graph_region"]) if e["graph_region"] else ""}</title></rect>')
    for s, e, lab in regions:                        # GetOrganelle graph: repeats short reads cannot count
        s, e = max(0, s), min(rlen, e)
        el.append(f'<rect x="{X(s):.1f}" y="{y1 - h - 4}" width="{max(X(e) - X(s), 2):.1f}" height="{2 * h + 8}" '
                  f'fill="{REPEAT}" opacity="0.35"><title>{escape(f"{s + 1}-{e}: {lab}")}</title></rect>')
    el.append(f'<text x="{left - 12}" y="{y1 + 5}" text-anchor="end" font-size="15">GetOrganelle ({rlen:,} bp)'
              '</text>')
    el.append(f'<line x1="{X(0)}" y1="{y1}" x2="{X(rlen)}" y2="{y1}" stroke="#555"/>')
    n_contigs = len(pos)
    el.append(f'<text x="{left - 12}" y="{y2 + 5}" text-anchor="end" font-size="15">{escape(label)} '
              f'({olen if n_contigs == 1 else sum(v[1] for v in pos.values()):,} bp'
              f'{f", {n_contigs} contigs" if n_contigs > 1 else ""})</text>')
    for name, (off, ln, flip) in pos.items():
        el.append(f'<line x1="{X(off):.1f}" y1="{y2}" x2="{X(off + ln):.1f}" y2="{y2}" stroke="#555" stroke-width="2">'
                  f'<title>{escape(name)} ({ln:,} bp{", reverse complement" if flip else ""})</title></line>')
    for s, e, gname, strand, typ in genes:
        gene(X(s % rlen), X(s % rlen + e - s), y1, strand, gname, typ, True, f"{gname} {s + 1}-{e}")
        for b in blocks:                             # carried to the other assembly through each alignment
            for sh in (0, rlen, -rlen):
                a, z = max(s + sh, b["rs"]), min(e + sh, b["re"])
                if z - a >= min(30, e - s):
                    flip = pos[b["contig"]][2]
                    st = strand * (-1 if b["strand"] == "-" else 1) * (-1 if flip else 1)
                    gene(qx(b["contig"], to_query(b, a)), qx(b["contig"], to_query(b, z)), y2, st, gname, typ,
                         False, f"{gname} on {b['contig']} ({100 * (z - a) / (e - s):.0f}% of the gene)")
    for d in diffs:
        s, e = int(d["ref_start"]) - 1, int(d["ref_end"])
        tip = (f"{d['type']} {d['ref_bp']} -> {d['other_bp']} bp at {d['ref_start']}-{d['ref_end']} | {d['gene']} | "
               f"{d['graph_region'] or 'outside flagged regions'}")
        xs, xe = X(s), max(X(e), X(s) + 2)
        el.append(f'<rect x="{xs:.1f}" y="{y1 - h - 3}" width="{xe - xs:.1f}" height="{2 * h + 6}" fill="{DIFF}">'
                  f'<title>{escape(tip)}</title></rect>')
        b = next((b for b in blocks for sh in (0, rlen) if b["rs"] <= s + sh <= b["re"]), None)
        if b:
            p = s if b["rs"] <= s <= b["re"] else s + rlen
            q = qx(b["contig"], to_query(b, p))
            w = max(2.0, int(d["other_bp"]) * k)
            el.append(f'<polygon points="{xs:.1f},{y1 + h} {xe:.1f},{y1 + h} {q + w:.1f},{y2 - h} {q:.1f},{y2 - h}" '
                      f'fill="{DIFF}" opacity="0.8"><title>{escape(tip)}</title></polygon>'
                      f'<rect x="{q:.1f}" y="{y2 - h - 3}" width="{w:.1f}" height="{2 * h + 6}" fill="{DIFF}">'
                      f'<title>{escape(tip)}</title></rect>')
    for e in events:                                 # misassembly junctions on the reference
        if e["event"] in ("inversion", "relocation", "local misassembly"):
            for p in e["ref_pos"].split("->"):
                x = X(int(p) - 1)
                el.append(f'<polygon points="{x - 6:.1f},{y1 - h - 16} {x + 6:.1f},{y1 - h - 16} {x:.1f},'
                          f'{y1 - h - 5}" fill="{DIFF}"><title>{escape(e["event"])} {e["ref_pos"]} '
                          f'({e["size"]} bp) {escape(e["graph_region"])}</title></polygon>')
    ya = y2 + 45                                     # scale
    step = 2000 if max(rlen, olen) > 8000 else 1000
    el.append(f'<line x1="{X(0)}" y1="{ya}" x2="{X(max(rlen, olen))}" y2="{ya}" stroke="#555"/>')
    for p in range(0, max(rlen, olen) + 1, step):
        el.append(f'<line x1="{X(p):.1f}" y1="{ya}" x2="{X(p):.1f}" y2="{ya + 5}" stroke="#555"/>'
                  f'<text x="{X(p):.1f}" y="{ya + 18}" text-anchor="middle" font-size="11">{p // 1000} kb</text>')
    legend = [(DIFF, "difference / misassembly"), (REPEAT, "GetOrganelle repeat region (graph)"),
              (MISSING, "not in this assembly"), (INVERTED, "inverted alignment"), (COLORS["CDS"], "CDS"),
              (COLORS["rRNA"], "rRNA"), (COLORS["tRNA"], "tRNA")]
    x = left
    for color, text in legend:
        el.append(f'<rect x="{x}" y="{ya + 32}" width="12" height="12" fill="{color}"/>'
                  f'<text x="{x + 17}" y="{ya + 43}" font-size="12">{text}</text>')
        x += 30 + 6.5 * len(text)
    head = f'<text x="{left}" y="28" font-size="16" font-weight="bold">{escape(title)}</text>' if title else ""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{ya + 60}" viewBox="0 0 {W} {ya + 60}" '
            f'font-family="Arial, Helvetica, sans-serif"><rect width="100%" height="100%" fill="white"/>{head}'
            + "".join(el) + "</svg>\n")

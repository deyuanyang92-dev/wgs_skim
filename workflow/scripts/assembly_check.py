"""QUAST-style check of an assembly against a reference mitogenome: the parts of QUAST 5.2 that matter for one
small circular genome, re-implemented without QUAST (119 MB source, GeneMark / BUSCO / Icarus / plots).

What is taken from QUAST (quast_libs, read 2026-10-28):
- alignment: minimap2 -c -x asm5 --cs with QUAST's own options (ca_utils/align_contigs.py, run_minimap): -r 200
  ends an alignment at any jump over 200 bp, so a moved or deleted block becomes a junction, not a long gap;
- misassemblies between consecutive alignments of one contig (ca_utils/analyze_misassemblies.py, is_misassembly):
  other strand = inversion; inconsistency = distance on reference - distance on contig, |inconsistency| > 1000 =
  relocation, > 200 = local misassembly, smaller = indel; circular reference: the shorter way round the circle;
- genes (genome_analyzer.py): complete = inside one alignment, partial = overlapping one (QUAST: >= 100 bp);
  here judged on the bases present (a gene can sit in one alignment with most of it deleted): complete >= 90 %,
  partial >= 10 %, missing below;
- genome fraction, duplication ratio, unaligned length.
Changed for mitogenomes: the reference is doubled when circular, so alignments and genes across its origin stay
whole; partial = any overlap (a tRNA is ~65 bp, QUAST's 100 bp would hide it); every junction is placed against
GetOrganelle's repeat regions (repeat copies short reads cannot count); genes also get the base differences inside
them (minimap2 cs tag) and, for protein genes, whether the indels break the reading frame.

Needs minimap2 on PATH.
"""
import re
import subprocess
import tempfile
from pathlib import Path
from urllib.parse import unquote

EXTENSIVE, LOCAL = 1000, 200        # QUAST defaults (qconfig: DEFAULT_EXT_MIS_SIZE, local_misassembly_min_length)
MIN_UNALIGNED = 50
BIG_INDEL = 50                      # indels inside an alignment listed with the junctions
QUAST_MINIMAP2 = ["-B5", "-O4,16", "--no-long-join", "-r", "200", "-s", "65", "-z", "200", "--mask-level", "0.9",
                  "--min-occ", "200", "-g", "2500", "--score-N", "2"]


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


def read_gff(path):
    """Genes from a GFF on the reference sequence (as QUAST -g): [(start, end, name, strand, type)], 0-based.
    Cleaned as in the user's modify_gff2.py (Quast/, 2024): region / source / gene rows are dropped when the file
    has CDS / rRNA / tRNA rows (Geneious writes both), Name= keeps its first word ('COX1 gene' -> COX1) without the
    URL-encoded anticodon ('trnN%28guu%29' -> trnN)."""
    rows = []
    for line in open(path):
        if line.startswith("##FASTA"):
            break
        t = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(t) < 9:
            continue
        attrs = dict(x.split("=", 1) for x in t[8].split(";") if "=" in x)
        name = unquote(next((attrs[k] for k in ("Name", "gene", "product", "ID") if attrs.get(k)), t[2]))
        name = name.split(" ")[0].split("(")[0]      # 'trnN%28guu%29' -> trnN (modify_gff2 cut at '%')
        rows.append((int(t[3]) - 1, int(t[4]), name, -1 if t[6] == "-" else 1, t[2]))
    typed = any(r[4] in ("CDS", "rRNA", "tRNA") for r in rows)
    out, seen = [], set()
    for r in rows:
        if r[4] in ("region", "source") or (typed and r[4] == "gene") or (r[0], r[1]) in seen:
            continue
        seen.add((r[0], r[1]))
        out.append(r)
    return out


def minimap2(ref, contigs, ref_circular, threads=4):
    """Alignments of every contig: dicts contig, qlen, qs, qe, strand, rs, re (reference coordinates, may run past
    the reference end on a circle), ops [(ref_pos, kind, length, bases)] from the cs tag (reference strand)."""
    n = len(ref)
    with tempfile.TemporaryDirectory() as d:
        r, q = Path(d) / "ref.fa", Path(d) / "q.fa"
        r.write_text(f">ref\n{ref + ref if ref_circular else ref}\n")
        q.write_text("".join(f">{name}\n{seq}\n" for name, seq, _ in contigs))
        res = subprocess.run(["minimap2", "-c", "-x", "asm5", *QUAST_MINIMAP2, "--cs", "--secondary=no",
                              "-t", str(threads), str(r), str(q)], capture_output=True, text=True, check=True)
    blocks = []
    for line in res.stdout.splitlines():
        t = line.split("\t")
        cs = next((x[5:] for x in t[12:] if x.startswith("cs:Z:")), "")
        rs, re_ = int(t[7]), int(t[8])
        if ref_circular and rs >= n:                 # second copy of the doubled reference: back to the first
            rs, re_ = rs - n, re_ - n
        ops, pos = [], rs
        for op, val in re.findall(r"([:*+\-~])([0-9a-z]+)", cs):
            if op == ":":
                pos += int(val)
            elif op == "*":
                ops.append((pos, "mismatch", 1, val))
                pos += 1
            elif op == "+":
                ops.append((pos, "insertion", len(val), val))
            elif op == "-":
                ops.append((pos, "deletion", len(val), val))
                pos += len(val)
        blocks.append(dict(contig=t[0], qlen=int(t[1]), qs=int(t[2]), qe=int(t[3]), strand=t[4], rs=rs, re=re_,
                           ops=ops))
    return blocks


def in_region(pos, regions, n, pad=50):
    pos %= n
    return "; ".join(lab for s, e, lab in regions if s - pad <= pos < e + pad or s - pad <= pos + n < e + pad)


def structure(blocks, contigs, n, ref_circular, regions):
    """Junctions between consecutive alignments of each contig (QUAST misassembly classes) and unaligned parts."""
    circ = {name: c for name, _, c in contigs}
    events = []
    by = {}
    for b in blocks:
        by.setdefault(b["contig"], []).append(b)
    for name, bs in by.items():
        bs.sort(key=lambda b: b["qs"])
        pairs = [(b1, b2, 0) for b1, b2 in zip(bs, bs[1:])]
        if circ.get(name) and len(bs) > 1:           # circular contig: last alignment continues into the first
            pairs.append((bs[-1], bs[0], bs[0]["qlen"]))
        for b1, b2, wrap in pairs:
            dq = b2["qs"] + wrap - b1["qe"]
            at = b1["re"] if b1["strand"] == "+" else b1["rs"]       # where the contig leaves the reference ...
            to = b2["rs"] if b2["strand"] == "+" else b2["re"]       # ... and where it comes back
            if b1["strand"] != b2["strand"]:
                kind, inc = "inversion", None
            else:
                dr = to - at if b1["strand"] == "+" else at - to
                if ref_circular:
                    dr = (dr + n // 2) % n - n // 2
                inc = dr - dq
                kind = "relocation" if abs(inc) > EXTENSIVE else "local misassembly" if abs(inc) > LOCAL else \
                    "indel" if inc else None
            if kind:
                events.append(dict(contig=name, event=kind, ref_pos=f"{at % n + 1}->{to % n + 1}",
                                   contig_pos=b1["qe"], size=inc if inc is not None else "",
                                   graph_region=in_region(at, regions, n) or in_region(to, regions, n)))
        for b in bs:                                 # long indels inside one alignment
            for p, kind, ln, _ in b["ops"]:
                if kind != "mismatch" and ln >= BIG_INDEL:
                    events.append(dict(contig=name, event=kind, ref_pos=str(p % n + 1), contig_pos="",
                                       size=ln if kind == "insertion" else -ln, graph_region=in_region(p, regions, n)))
    cov = coverage(blocks, n)                        # reference stretches no contig covers
    i = 0
    while i < n:
        if cov[i]:
            i += 1
            continue
        j = i
        while j < n and not cov[j]:
            j += 1
        if j - i >= MIN_UNALIGNED:
            events.append(dict(contig="", event="not in assembly", ref_pos=f"{i + 1}-{j}", contig_pos="",
                               size=-(j - i), graph_region=in_region(i, regions, n) or in_region(j - 1, regions, n)))
        i = j
    unaligned = []
    for name, seq, _ in contigs:
        cov = [0] * len(seq)
        for b in by.get(name, []):
            for i in range(b["qs"], b["qe"]):
                cov[i] = 1
        i = 0
        while i < len(seq):
            if not cov[i]:
                j = i
                while j < len(seq) and not cov[j]:
                    j += 1
                if j - i >= MIN_UNALIGNED:
                    unaligned.append((name, i, j))
                i = j
            else:
                i += 1
    return events, unaligned


def coverage(blocks, n):
    """Per reference base: indices of the alignments covering it."""
    cov = [[] for _ in range(n)]
    for k, b in enumerate(blocks):
        for p in range(b["rs"], b["re"]):
            cov[p % n].append(k)
    return cov


def genes_check(feats, blocks, n, regions):
    """One row per gene: complete / partial / split / missing (QUAST classes, see module doc), strand, base
    differences inside the gene, frameshift for protein genes."""
    cov = coverage(blocks, n)
    main = {}                                        # each contig's own direction (fragments come either way)
    for c in {b["contig"] for b in blocks}:
        main[c] = max(("+", "-"), key=lambda s: sum(b["re"] - b["rs"] for b in blocks
                                                    if b["contig"] == c and b["strand"] == s))
    rows = []
    for s, e, name, _, typ in feats:
        pos = [p % n for p in range(s, e)]
        hit = [k for p in pos for k in cov[p]]
        covered = sum(1 for p in pos if cov[p])
        whole = [k for k in set(hit) if all(k in cov[p] for p in pos)]
        ks = whole[:1] if whole else sorted(set(hit))
        mis = ins = dele = 0
        span = set(pos)
        for k in ks:
            for p, kind, ln, _ in blocks[k]["ops"]:
                if kind == "mismatch":
                    mis += p % n in span
                elif kind == "insertion":                # between p - 1 and p: inside when both are
                    ins += ln if (p - 1) % n in span and p % n in span else 0
                else:
                    dele += sum((p + i) % n in span for i in range(ln))
        t = gene_type(name, typ)
        present = (covered - dele) / max(len(pos), 1)
        status = "missing" if present < 0.1 else "partial" if present < 0.9 else "complete" if whole else \
            "split" if covered == len(pos) else "partial"
        rows.append(dict(gene=name, type=t, ref_start=s + 1, ref_end=e, status=status,
                         present_pct=f"{100 * present:.1f}",
                         strand="reverse" if ks and all(blocks[k]["strand"] != main[blocks[k]["contig"]]
                                                        for k in ks) else "",
                         mismatches=mis, inserted_bp=ins, deleted_bp=dele,
                         frameshift="yes" if t == "CDS" and (ins - dele) % 3 else "",
                         graph_region=in_region(s, regions, n) or in_region(e - 1, regions, n)))
    return rows


def check(ref, ref_circular, contigs, feats, regions, threads=4):
    """Everything for one assembly: (stats dict, structure events, unaligned parts, gene rows, alignments).
    contigs = [(name, sequence, circular)]."""
    n = len(ref)
    blocks = minimap2(ref, contigs, ref_circular, threads)
    events, unaligned = structure(blocks, contigs, n, ref_circular, regions)
    cov = coverage(blocks, n)
    covered = sum(1 for c in cov if c)
    genes = genes_check(feats, blocks, n, regions) if feats else []
    count = lambda *kinds: sum(e["event"] in kinds for e in events)  # noqa: E731
    stats = dict(contigs=len(contigs), genome_fraction_pct=f"{100 * covered / n:.2f}",
                 duplication_ratio=f"{sum(len(c) for c in cov) / max(covered, 1):.3f}",
                 unaligned_bp=sum(j - i for _, i, j in unaligned),
                 misassemblies=count("inversion", "relocation"), local_misassemblies=count("local misassembly"),
                 misassemblies_outside_repeats=sum(e["event"] in ("inversion", "relocation", "local misassembly")
                                                   and not e["graph_region"] for e in events),
                 genes_complete=sum(g["status"] == "complete" for g in genes),
                 genes_partial=sum(g["status"] in ("partial", "split") for g in genes),
                 genes_missing=sum(g["status"] == "missing" for g in genes),
                 genes_with_differences=sum(bool(g["mismatches"] or g["inserted_bp"] or g["deleted_bp"])
                                            for g in genes),
                 frameshift_genes=",".join(g["gene"] for g in genes if g["frameshift"]))
    return stats, events, unaligned, genes, blocks

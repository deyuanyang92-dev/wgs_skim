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
Batch:  compare_mt.py --results results_dir [--samples A,B] [-j 2] -o out/    (layouts in LAYOUTS)
        compare_mt.py --sheet samples.tsv -o out/       (columns sample, ref, other, [spades, gb, gfa])
       -> per-sample files as above, plus out/all.compare.tsv and out/all.diffs.tsv (skipped / error rows included)
Figures (plot_compare.py, plain SVG): <prefix>.<label>.compare.svg per comparison - genes of the reference
annotation on both assemblies, alignment links, differences red, GetOrganelle repeat regions purple; batch also
writes out/all.compare.html (table + every figure). --no-plot to skip.
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
from plot_compare import draw  # noqa: E402
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
    """[(start, end, name, strand, type)] 0-based on ref from a GenBank file (gene / CDS / tRNA / rRNA / D-loop / misc_feature);
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
                    cur = [min(pos) - 1, max(pos), m.group(1), -1 if "complement" in m.group(2) else 1, m.group(1)]
                    out.append(cur)
        elif cur and re.search(r'/(gene|product|note)="', line) and cur[2] in ("gene", "CDS", "tRNA", "rRNA", "D-loop",
                                                                           "misc_feature", "rep_origin"):
            cur[2] = re.search(r'="([^"]*)', line).group(1)
    n = len(gseq) or len(ref)
    seen, res = set(), []
    for a, b, name, strand, typ in out:              # gene + CDS of the same gene: keep the first
        if (a, b) not in seen:
            seen.add((a, b))
            res.append(((a - shift) % n, (a - shift) % n + b - a, name, strand, typ))
    return res


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
        g = sorted({f[2] for f in feats if f[0] <= e and f[1] >= s})
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


def compare_sample(ref_path, others, out, prefix, spades=None, gb=None, gfa=None, plot=True):
    """Compare one sample's assemblies with its GetOrganelle result; writes the per-sample files.
    Returns (summary rows, difference rows, printed lines)."""
    out.mkdir(parents=True, exist_ok=True)
    fasta, ref_gfa, alts = getorganelle_result(ref_path)
    (header, ref), = read_fasta(fasta)[:1]
    ref_circ = "circular" in header and "linear" not in header
    regions = graph_regions(header, gfa or ref_gfa, ref)
    feats = genes(gb, ref) if gb else []
    with open(out / f"{prefix}.genes.tsv", "w") as fh:          # for the figure
        fh.write("start\tend\tname\tstrand\ttype\n" + "".join("\t".join(map(str, f)) + "\n" for f in feats))
    with open(out / f"{prefix}.regions.tsv", "w") as fh:
        fh.write("start\tend\tlabel\n" + "".join(f"{x}\t{y}\t{lab}\n" for x, y, lab in regions))
    summaries, diffs, lines = [], [], []
    for item in others:
        label, path = item.split("=", 1) if "=" in item and not Path(item).exists() else (Path(item).stem, item)
        recs = read_fasta(path)
        if any("circular" in n and "linear" not in n for n, _ in recs) and len(recs) == 1:
            seq, circ, method = recs[0][1], True, "as given (circular)"
        else:
            rows, seqs = circularize(recs, spades)
            done = [(r, s) for r, (_, s) in zip(rows, seqs) if r["method"] != "not_circular"]
            if rows and rows[-1]["method"] == "join":
                done = [(rows[-1], seqs[0][1])]
            if done:
                r, seq = max(done, key=lambda x: len(x[1]))
                circ, method = True, f"circularized: {r['method']} ({r['note'] or r['contig']})"[:300]
            else:
                seq, circ, method = max((s for _, s in recs), key=len), False, "not circular (longest contig)"
        with open(out / f"{prefix}.{label}.oriented.fasta", "w") as fh:
            write_fasta([(f"{prefix}|{label} {'circular' if circ else 'linear'}", orient(ref, seq, circ)[0])], fh)
        s, rows, aln = compare(ref, ref_circ, regions, feats, label, seq, circ, method)
        s.update(sample=prefix, ref_alternative_paths=alts,
                 ref_problem_regions="; ".join(f"{x + 1}-{y}: {lab}" for x, y, lab in regions))
        summaries.append(s)
        diffs += [{"sample": prefix, **r} for r in rows]
        if aln:
            with open(out / f"{prefix}.{label}.aligned.fasta", "w") as fh:
                write_fasta([(f"{prefix}|getorganelle", aln[0]), (f"{prefix}|{label}", aln[1])], fh)
        lines.append(f"{prefix} getorganelle {len(ref)} bp vs {label} {len(seq)} bp: {s['verdict']} "
                     f"(identity {s.get('identity_pct')}%, {s.get('diff_blocks')} difference block(s), "
                     f"{s.get('unexplained_blocks')} outside GetOrganelle-flagged regions)")
        lines += [f"  {r['ref_start']}-{r['ref_end']} {r['type']} {r['ref_bp']}->{r['other_bp']} bp  {r['gene']}  "
                  f"{r['graph_region'] or '-'}" for r in rows]
    write_tables(out, prefix, summaries, diffs)
    if plot:
        for r in summaries:
            if (out / f"{prefix}.{r['other']}.aligned.fasta").exists():
                lines.append(figure(out, prefix, r, feats, regions, diffs))
    return summaries, diffs, lines


def figure(out, prefix, row, feats, regions, diffs):
    """SVG of one comparison (plot_compare.draw) -> <prefix>.<label>.compare.svg."""
    base = out / f"{prefix}.{row['other']}"
    aligned = [s for _, s in read_fasta(f"{base}.aligned.fasta")][:2]
    svg = draw(aligned, feats, regions, [d for d in diffs if d["other"] == row["other"]], row["other"],
               f"{prefix}: {row['verdict']} (identity {row.get('identity_pct')}%)")
    (out / f"{base.name}.compare.svg").write_text(svg)
    return f"  figure: {base}.compare.svg"


def write_tables(out, prefix, summaries, diffs):
    for name, fields, data in (("compare", SUMMARY, summaries), ("diffs", DIFFS, diffs)):
        with open(out / f"{prefix}.{name}.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore", restval="")
            w.writeheader()
            w.writerows(data)


# Result layouts found by --results: paths relative to the results directory, {s} = sample.
# wgs_skim_pipeline (Temp_scripts): GetOrganelle 02.organelle/02_mt, SPAdes 03.spades, MitoFinder 10.mitofinder.
LAYOUTS = {
    "wgs_skim_pipeline": {
        "samples": "02.organelle/02_mt/*",
        "ref": "02.organelle/02_mt/{s}",
        "other": ["mitofinder=10.mitofinder/{s}/mitofinder/mitofinder_summary/mt_contigs/{s}_all_mt.fasta"],
        "spades": "03.spades/{s}",
        "gb": "06.mito_annotation/final_annotation/{s}_locus_mt_1.gb",
    },
}


def jobs_from_results(results):
    """One job per sample found in a results directory of a known layout."""
    res = Path(results)
    for name, lay in LAYOUTS.items():
        dirs = sorted(d for d in res.glob(lay["samples"]) if d.is_dir())
        if not dirs:
            continue
        print(f"[compare] {res}: layout {name}, {len(dirs)} sample(s)", file=sys.stderr)
        jobs = []
        for d in dirs:
            s = d.name
            f = lambda t: str(res / t.format(s=s))  # noqa: E731
            jobs.append({"sample": s, "ref": f(lay["ref"]), "spades": f(lay["spades"]), "gb": f(lay["gb"]),
                         "other": [o.split("=", 1)[0] + "=" + f(o.split("=", 1)[1]) for o in lay["other"]]})
        return jobs
    sys.exit(f"[compare] {res}: no known result layout ({', '.join(LAYOUTS)}); use --sheet")


def jobs_from_sheet(sheet):
    """Sample sheet (TSV): sample, ref, other ([label=]fasta, several separated by ','), optional spades, gb, gfa."""
    jobs = []
    with open(sheet) as fh:
        for row in csv.DictReader((x for x in fh if x.strip() and not x.startswith("#")), delimiter="\t"):
            jobs.append({"sample": row["sample"], "ref": row["ref"], "other": row["other"].split(","),
                         **{k: row.get(k) or None for k in ("spades", "gb", "gfa")}})
    return jobs


def missing_inputs(job):
    """Why a job cannot run (required files absent); optional files that are absent are dropped."""
    for k in ("spades", "gb", "gfa"):
        if job.get(k) and not Path(job[k]).exists():
            job[k] = None
    miss = [] if Path(job["ref"]).exists() else ["GetOrganelle result"]
    miss += [o.split("=", 1)[0] for o in job["other"] if not Path(o.split("=", 1)[-1]).exists()]
    return miss


def run_job(job, out):
    try:
        s, d, lines = compare_sample(job["ref"], job["other"], out, job["sample"], job.get("spades"), job.get("gb"),
                                     job.get("gfa"), job.get("plot", True))
        return s, d, lines
    except Exception as e:  # one sample failing must not stop the batch; reported in the table
        return [{"sample": job["sample"], "verdict": "error", "note": f"{type(e).__name__}: {e}"[:300]}], [], \
            [f"{job['sample']}: error: {type(e).__name__}: {e}"]


def batch(jobs, out, threads, only=None, plot=True):
    """All samples, `threads` at a time; per-sample files plus all.compare.tsv / all.diffs.tsv."""
    from concurrent.futures import ProcessPoolExecutor
    if only:
        keep = set(only.split(","))
        jobs = [j for j in jobs if j["sample"] in keep]
    summaries, diffs, todo = [], [], []
    for j in jobs:
        j["plot"] = plot
        miss = missing_inputs(j)
        if miss:
            summaries.append({"sample": j["sample"], "verdict": "skipped", "note": "missing: " + ", ".join(miss)})
            print(f"{j['sample']}: skipped (missing: {', '.join(miss)})")
        else:
            todo.append(j)
    with ProcessPoolExecutor(max_workers=max(1, threads)) as ex:
        for s, d, lines in ex.map(run_job, todo, [out] * len(todo)):
            summaries += s
            diffs += d
            print("\n".join(lines), flush=True)
    order = {j["sample"]: i for i, j in enumerate(jobs)}
    summaries.sort(key=lambda r: order.get(r["sample"], 0))
    diffs.sort(key=lambda r: order.get(r["sample"], 0))
    write_tables(out, "all", summaries, diffs)
    overview(out, summaries)
    n = {}
    for r in summaries:
        n[r["verdict"]] = n.get(r["verdict"], 0) + 1
    print(f"[compare] {len(jobs)} sample(s): " + ", ".join(f"{k} {v}" for k, v in n.items())
          + f" -> {out / 'all.compare.tsv'}, {out / 'all.diffs.tsv'}, {out / 'all.compare.html'}")


def overview(out, summaries):
    """One page: summary table, then each sample's figure (inline SVG)."""
    from html import escape
    cols = ["sample", "other", "ref_length", "other_length", "ref_topology", "other_topology", "identity_pct",
            "diff_blocks", "unexplained_blocks", "verdict", "note"]
    color = {"identical": "#2e7d32", "consistent": "#1565c0", "differ": "#c62828"}
    rows = "".join("<tr>" + "".join(
        f'<td style="color:{color.get(r.get(c), "inherit")}">{escape(str(r.get(c, "")))}</td>' for c in cols)
        + "</tr>" for r in summaries)
    figs = ""
    for r in summaries:
        f = out / f"{r['sample']}.{r.get('other', '')}.compare.svg"
        if r.get("other") and f.exists():
            figs += f'<h3>{escape(r["sample"])} - {escape(r["other"])}</h3>' + f.read_text()
    (out / "all.compare.html").write_text(
        '<!doctype html><meta charset="utf-8"><title>Assembly comparison</title><style>body{font-family:Arial,'
        'sans-serif;margin:20px}table{border-collapse:collapse;font-size:13px}td,th{border:1px solid #ccc;'
        'padding:3px 8px}svg{max-width:100%;height:auto}</style><h2>GetOrganelle vs other assemblies</h2>'
        '<p>identical: no difference; consistent: every difference is in a GetOrganelle repeat region, at a linear '
        'end or a homopolymer length; differ: check by hand. Hover the figures for details.</p>'
        f'<table><tr>{"".join(f"<th>{c}</th>" for c in cols)}</tr>{rows}</table>{figs}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-a", "--ref", help="GetOrganelle output directory or FASTA (reference)")
    ap.add_argument("-b", "--other", action="append", help="[label=]FASTA to compare (repeatable)")
    ap.add_argument("--spades", help="SPAdes directory of the other assembly (graph-based circularization)")
    ap.add_argument("--gb", help="GenBank annotation of the reference (genes for each difference)")
    ap.add_argument("--gfa", help="reference assembly graph, if -a is a FASTA")
    ap.add_argument("--results", help="batch: results directory of a known layout (all samples found there)")
    ap.add_argument("--sheet", help="batch: sample sheet TSV (sample, ref, other, [spades, gb, gfa])")
    ap.add_argument("--samples", help="batch: only these samples (comma-separated)")
    ap.add_argument("-j", "--jobs", type=int, default=2,
                    help="batch: samples at a time (each needs ~7 GB memory with --spades graphs; default 2)")
    ap.add_argument("--no-plot", action="store_true", help="no figures (default: one SVG per comparison)")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--prefix", default="sample")
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    if a.results or a.sheet:
        batch(jobs_from_results(a.results) if a.results else jobs_from_sheet(a.sheet), out, a.jobs, a.samples,
              not a.no_plot)
        return
    if not (a.ref and a.other):
        ap.error("give -a and -b (one sample), or --results / --sheet (batch)")
    print("\n".join(compare_sample(a.ref, a.other, out, a.prefix, a.spades, a.gb, a.gfa, not a.no_plot)[2]))


if __name__ == "__main__":
    main()

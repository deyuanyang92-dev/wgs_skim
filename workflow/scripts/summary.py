"""Collect per-sample results into summary/samples.tsv and summary/report.html (Snakemake script for rule summary)."""
import csv
import html
import json
import re
from datetime import datetime
from pathlib import Path

sm = snakemake  # noqa: F821
inp = sm.input
samples = sm.params.samples


def by_sample(files):
    """Map sample -> file using the per-sample path pattern <...>/<sample>/<sample>.<ext>."""
    out = {}
    for f in files or []:
        p = Path(f)
        out[p.parent.name] = p
    return out


def fasta_lengths(path):
    lens, n = [], 0
    for line in open(path):
        if line.startswith(">"):
            if n:
                lens.append(n)
            n = 0
        else:
            n += len(line.strip())
    if n:
        lens.append(n)
    return lens


def n50(lens):
    if not lens:
        return ""
    lens, half, acc = sorted(lens, reverse=True), sum(lens) / 2, 0
    for x in lens:
        acc += x
        if acc >= half:
            return x


def busco_line(path):
    m = re.search(r"C:[\d.]+%\[S:[\d.]+%,D:[\d.]+%\],F:[\d.]+%,M:[\d.]+%,n:\d+", Path(path).read_text())
    return m.group(0) if m else ""


def genomescope(path):
    t = Path(path).read_text()
    size = re.search(r"Genome Haploid Length\s+([\d,]+) bp\s+([\d,]+) bp", t)
    het = re.search(r"Heterozygous \(ab\)\s+([\d.]+)%\s+([\d.]+)%", t)
    size_mb = f"{int(size.group(2).replace(',', '')) / 1e6:.1f}" if size else ""
    return size_mb, (f"{het.group(2)}" if het else "")


def gene_count(path):
    t = Path(path).read_text() if Path(path).exists() else ""
    return {k: len(re.findall(rf"^ {{5}}{k} ", t, re.M)) for k in ("CDS", "tRNA", "rRNA")}


cols = ["sample", "clean_gb", "q30", "mt_length", "mt_circular", "mt_source", "mt_genes", "mt_topology", "mt_end_diagnosis",
        "nr_length", "nr_topology", "spades_n50", "spades_total_mb", "busco", "genome_size_mb",
        "heterozygosity_pct", "blast_mt", "blast_nr"]
rows = {s: {c: "" for c in cols} | {"sample": s} for s in samples}

for s, f in by_sample(inp.fastp).items():
    j = json.load(open(f))["summary"]["after_filtering"]
    rows[s].update(clean_gb=f"{j['total_bases'] / 1e9:.2f}", q30=f"{j['q30_rate'] * 100:.1f}")
for s, f in by_sample(inp.get("mt")).items():
    j = json.load(open(f))
    rows[s].update(mt_length=j["length"] or "", mt_circular="yes" if j["circular"] else ("no" if j["length"] else "failed"),
                   mt_source=j["note"] if j.get("source") == "spades" else "; ".join(x for x in (
                       f"GetOrganelle {j['round']}" + (f" ({j['stop']})" if j.get("stop") else ""), j.get("note")) if x))
for s, f in by_sample(inp.get("nr")).items():
    lens = fasta_lengths(f)
    rows[s]["nr_length"] = ",".join(map(str, lens))
for s, f in by_sample(inp.get("scaffolds")).items():
    lens = fasta_lengths(f) if Path(f).stat().st_size else []
    rows[s].update(spades_n50=n50(lens), spades_total_mb=f"{sum(lens) / 1e6:.1f}" if lens else "")
for s, f in by_sample(inp.get("busco")).items():
    rows[s]["busco"] = busco_line(f)
for s, f in by_sample(inp.get("genomescope")).items():
    rows[s]["genome_size_mb"], rows[s]["heterozygosity_pct"] = genomescope(f)
for s, f in by_sample(inp.get("annotation")).items():
    g = gene_count(f)
    rows[s]["mt_genes"] = f"{g['CDS']} PCG / {g['tRNA']} tRNA / {g['rRNA']} rRNA" if any(g.values()) else ""
for s, f in by_sample(inp.get("topology")).items():
    if Path(f).stat().st_size:
        for r in csv.DictReader(open(f), delimiter="\t"):
            key = "mt_topology" if "|mt|" in r["contig"] else "nr_topology"
            rows[s][key] = ",".join(x for x in (rows[s][key], r["verdict"]) if x)
            if key == "mt_topology" and r.get("end_diagnosis"):
                rows[s]["mt_end_diagnosis"] = r["end_diagnosis"]
if inp.get("blast") and Path(inp.blast).stat().st_size:
    best = {}
    for r in csv.DictReader(open(inp.blast, encoding="utf-8-sig"), delimiter="\t"):
        best.setdefault(r["qseqid"], r)       # rows are sorted best first
    for q, r in best.items():
        s, locus, what = (q.split("|") + ["", ""])[:3]
        if s in rows and locus in ("mt", "nr"):
            hit = f"{what}: {r.get('organism', '')} ({float(r['pident']):.1f}%, qcov {r.get('qcovs', '')}%)"
            rows[s][f"blast_{locus}"] = "; ".join(x for x in (rows[s][f"blast_{locus}"], hit) if x)

with open(sm.output.tsv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(rows[s] for s in samples)

n_circ = sum(r["mt_circular"] == "yes" for r in rows.values())
head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
body = "".join("<tr>" + "".join(f"<td>{html.escape(str(r[c]))}</td>" for c in cols) + "</tr>" for r in rows.values())
Path(sm.output.html).write_text(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>wgs_skim report</title>
<style>
:root{{--bg:#fff;--fg:#1d2433;--muted:#5b6475;--line:#e3e7ee;--head:#f3f5f9}}
@media (prefers-color-scheme:dark){{:root{{--bg:#14171c;--fg:#e6e9ef;--muted:#a3abb9;--line:#2b313b;--head:#1d2128}}}}
body{{margin:0;padding:24px 16px;background:var(--bg);color:var(--fg);font:14px/1.5 system-ui,sans-serif}}
h1{{font-size:20px;margin:0 0 4px}} p{{color:var(--muted);margin:0 0 16px}}
.wrap{{overflow-x:auto;border:1px solid var(--line);border-radius:8px}}
table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}}
th,td{{padding:6px 10px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}}
th{{background:var(--head);position:sticky;top:0}}
</style></head><body>
<h1>wgs_skim report</h1>
<p>{len(rows)} samples · taxon preset {html.escape(sm.params.taxon)} · mitogenome circular {n_circ}/{len(rows)} ·
generated {datetime.now():%Y-%m-%d %H:%M}</p>
<div class="wrap"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>
</body></html>
""", encoding="utf-8")

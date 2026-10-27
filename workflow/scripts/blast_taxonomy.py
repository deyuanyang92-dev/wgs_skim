"""Species check: BLAST mitogenome and rDNA queries of every sample against nt, with GenBank organism and taxonomy.

Snakemake script (rule blast): local nt if params.db is set, otherwise NCBI online BLAST (blastn -remote).
Queries (params.mt / params.nr, several allowed): "full" = the whole assembled sequence, or genes cut from the
annotations - mitogenome genes from MitoM (cox1 = COI barcode, 16s, 12s, cytb, nad1 ...), rDNA regions from
`mitom nuc` (18S, ITS1, 5.8S, ITS2, 28S; ITS = ITS1-5.8S-ITS2). A sample without the requested gene gets the
whole sequence instead (logged). Query IDs: <sample>|mt|cox1, <sample>|nr|ITS, <sample>|mt|full ...
Output: one row per (query, hit), best first: the TopN distinct accessions per query by bitscore.
"""
import csv
import subprocess
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO

sm = snakemake  # noqa: F821
sys.stderr = open(sm.log[0], "w")
out_dir = Path(sm.output.tsv).parent
FIELDS = "qseqid sacc stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs".split()
GB_FIELDS = ["organism", "taxonomy", "definition", "seq_length", "topology", "reference"]


ALIAS = {"coi": "cox1", "co1": "cox1", "cob": "cytb", "16s": "rrnl", "lrrna": "rrnl", "12s": "rrns", "srrna": "rrns",
         "nd1": "nad1", "nd2": "nad2", "nd3": "nad3", "nd4": "nad4", "nd4l": "nad4l", "nd5": "nad5", "nd6": "nad6"}
NR_PRODUCT = {"18S ribosomal RNA": "18S", "5.8S ribosomal RNA": "5.8S", "28S ribosomal RNA": "28S",
              "internal transcribed spacer 1": "ITS1", "internal transcribed spacer 2": "ITS2"}


def norm(name):
    n = name.strip().lower().replace("-", "")
    return ALIAS.get(n, n)


def mt_genes(gb):
    """{gene: longest sequence} of the CDS / rRNA in a MitoM GenBank file."""
    genes = {}
    for rec in SeqIO.parse(gb, "genbank"):
        for f in rec.features:
            if f.type in ("CDS", "rRNA") and "gene" in f.qualifiers:
                g, seq = norm(f.qualifiers["gene"][0]), str(f.extract(rec.seq))
                if len(seq) > len(genes.get(g, "")):
                    genes[g] = seq
    return genes


def nr_regions(gb):
    """{18S, ITS1, 5.8S, ITS2, 28S, ITS: longest copy} from the `mitom nuc` GenBank file."""
    regions = {}
    for rec in SeqIO.parse(gb, "genbank"):
        feats = [(NR_PRODUCT[f.qualifiers["product"][0]], f) for f in rec.features
                 if f.qualifiers.get("product", [""])[0] in NR_PRODUCT]
        for g, f in feats:
            seq = str(f.extract(rec.seq))
            if len(seq) > len(regions.get(g.lower(), "")):
                regions[g.lower()] = seq
        its = [f for g, f in feats if g in ("ITS1", "5.8S", "ITS2")]
        if its:                                  # ITS1-5.8S-ITS2 as one region (one copy: features are adjacent)
            lo, hi = min(int(f.location.start) for f in its), max(int(f.location.end) for f in its)
            seq = rec.seq[lo:hi]
            seq = str(seq.reverse_complement() if its[0].location.strand == -1 else seq)
            if len(seq) > len(regions.get("its", "")) and hi - lo < 3000:
                regions["its"] = seq
    return regions


def by_sample(files):
    return {s: f for f in files for s in sm.params.samples if Path(f).name.startswith(f"{s}.")}


def queries(locus, fastas, gbs, wanted, extract):
    """[(id, seq)] for one locus of every sample; the whole sequence (once) when a wanted gene is missing."""
    out = []
    for s, f in by_sample(fastas).items():
        recs = read_fasta(f)
        whole = [(f"{s}|{locus}|full" + (str(i) if len(recs) > 1 else ""), q) for i, (_, q) in enumerate(recs, 1)]
        genes = extract(gbs[s]) if s in gbs and Path(gbs[s]).stat().st_size else {}
        missing = [w for w in wanted if w != "full" and norm(w) not in genes]
        out += [(f"{s}|{locus}|{w}", genes[norm(w)]) for w in wanted if w != "full" and norm(w) in genes]
        if "full" in wanted or missing:
            out += whole
        if missing and whole:
            print(f"{s}: {', '.join(missing)} not in the {locus} annotation: whole sequence used", file=sys.stderr)
    return out


def read_fasta(path):
    recs = []
    for line in open(path):
        if line.startswith(">"):
            recs.append([line[1:].split()[0], []])
        elif recs:
            recs[-1][1].append(line.strip())
    return [(n, "".join(q)) for n, q in recs]


query = out_dir / "query.fasta"
qs = queries("mt", sm.input.get("mt", []), by_sample(sm.input.get("gb", [])), sm.params.mt, mt_genes) + \
    queries("nr", sm.input.get("nr", []), by_sample(sm.input.get("nr_gb", [])), sm.params.nr, nr_regions)
with open(query, "w") as fh:
    fh.writelines(f">{n}\n{q}\n" for n, q in qs if q)


rows = []
if query.stat().st_size:
    raw = out_dir / "blast.out"
    cmd = ["blastn", "-query", str(query), "-outfmt", "6 " + " ".join(FIELDS), "-evalue", "1e-10",
           "-max_target_seqs", str(max(50, sm.params.topn * 10)), "-out", str(raw)]
    cmd += (["-db", sm.params.db, "-num_threads", str(sm.threads)] if sm.params.db else ["-db", "nt", "-remote"])
    print(" ".join(cmd), file=sys.stderr, flush=True)
    subprocess.run(cmd, check=True, stderr=sys.stderr)

    hits = {}
    for line in open(raw):
        r = dict(zip(FIELDS, line.rstrip("\n").split("\t")))
        best = hits.setdefault(r["qseqid"], {})
        old = best.get(r["sacc"])
        if old is None or float(r["bitscore"]) > float(old["bitscore"]):   # best HSP per accession
            best[r["sacc"]] = r
    for q in hits:
        rows += sorted(hits[q].values(), key=lambda r: -float(r["bitscore"]))[: sm.params.topn]

    # GenBank records of the hits, 100 per request
    Entrez.email = sm.params.email or "wgs_skim@example.org"
    meta, accs = {}, sorted({r["sacc"] for r in rows})
    for i in range(0, len(accs), 100):
        for attempt in range(3):
            try:
                with Entrez.efetch(db="nucleotide", id=",".join(accs[i:i + 100]), rettype="gb", retmode="text") as h:
                    for rec in SeqIO.parse(h, "genbank"):
                        a = rec.annotations
                        refs = a.get("references") or []
                        meta[rec.id.split(".")[0]] = {
                            "organism": a.get("organism", ""), "taxonomy": "; ".join(a.get("taxonomy", [])),
                            "definition": rec.description, "seq_length": len(rec.seq),
                            "topology": a.get("topology", ""), "reference": refs[0].title if refs else ""}
                break
            except Exception as e:                       # network: retry, then leave the metadata empty
                print(f"efetch failed ({e}), attempt {attempt + 1}/3", file=sys.stderr, flush=True)
                time.sleep(10 * (attempt + 1))
    for r in rows:
        r.update(meta.get(r["sacc"].split(".")[0], dict.fromkeys(GB_FIELDS, "")))

with open(sm.output.tsv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS + GB_FIELDS, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)
print(f"{len(rows)} hits for {len({r['qseqid'] for r in rows})} sequences", file=sys.stderr)

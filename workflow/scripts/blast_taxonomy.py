"""Species check: BLAST mitogenome and rDNA queries of every sample against nt, with GenBank organism and taxonomy.

Snakemake script (rule blast): local nt if params.db is set, otherwise NCBI online BLAST (blastn -remote).
Queries (params.mt / params.nr, several allowed): "full" = the whole assembled sequence, or genes cut from the
annotations - mitogenome genes from MitoM (cox1 = COI barcode, 16s, 12s, cytb, nad1 ...), rDNA regions from
`mitom nuc` (18S, ITS1, 5.8S, ITS2, 28S; ITS = ITS1-5.8S-ITS2). A sample without the requested gene gets the
whole sequence instead (logged). Query IDs: <sample>|mt|cox1, <sample>|nr|ITS, <sample>|mt|full ...
Output: one row per (query, hit), best first: the TopN distinct accessions per query by bitscore.
The search + metadata part (blast_topn) is also used by mitofinder_check.py, and runs on any FASTA:
    blast_taxonomy.py --query q.fasta -o hits.tsv [--db LOCAL_DB] [--topn 5] [--threads 8] [--email E]
"""
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))   # the shared modules next to this script
if "snakemake" in globals():                                  # a Snakemake script: run from a temporary copy
    sys.path.insert(0, snakemake.scriptdir)  # noqa: F821
from samples import by_sample  # noqa: E402
from seqtools import read_fasta  # noqa: E402
from tables import write_tsv  # noqa: E402

FIELDS = "qseqid sacc stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs slen".split()
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
    from Bio import SeqIO                     # Biopython only where needed: the module imports without it
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
    from Bio import SeqIO                     # Biopython only where needed: the module imports without it
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


def queries(locus, fastas, gbs, wanted, extract):
    """[(id, seq)] for one locus of every sample; the whole sequence (once) when a wanted gene is missing."""
    out = []
    for s, f in by_sample(fastas).items():
        recs = read_fasta(f, first_word=True)
        whole = [(f"{s}|{locus}|full" + (str(i) if len(recs) > 1 else ""), q) for i, (_, q) in enumerate(recs, 1)]
        genes = extract(gbs[s]) if s in gbs and Path(gbs[s]).stat().st_size else {}
        missing = [w for w in wanted if w != "full" and norm(w) not in genes]
        out += [(f"{s}|{locus}|{w}", genes[norm(w)]) for w in wanted if w != "full" and norm(w) in genes]
        if "full" in wanted or missing:
            out += whole
        if missing and whole:
            print(f"{s}: {', '.join(missing)} not in the {locus} annotation: whole sequence used", file=sys.stderr)
    return out


def blast_topn(query, out_tsv, db="", topn=5, threads=8, email=""):
    """blastn of every query against a local database (or NCBI nt online), the TopN distinct accessions per query
    by bitscore with their GenBank organism / taxonomy / title; writes out_tsv, returns the rows"""
    query, out_tsv = Path(query), Path(out_tsv)
    rows = []
    if query.stat().st_size:
        raw = out_tsv.parent / (out_tsv.stem + ".blast.out")
        cmd = ["blastn", "-query", str(query), "-outfmt", "6 " + " ".join(FIELDS), "-evalue", "1e-10",
               "-max_target_seqs", str(max(50, topn * 10)), "-out", str(raw)]
        cmd += (["-db", db, "-num_threads", str(threads)] if db else ["-db", "nt", "-remote"])
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
            rows += sorted(hits[q].values(), key=lambda r: -float(r["bitscore"]))[:topn]

        # GenBank records of the hits, 100 per request; seq_start/seq_stop = 1: the annotation header only
        # (organism, taxonomy, title, references) - rDNA hits are often chromosome-scale records, whose sequences
        # would be downloaded otherwise (time-outs, then no metadata); the length comes from BLAST (slen)
        from Bio import Entrez, SeqIO
        Entrez.email = email or "wgs_skim@example.org"
        meta, accs = {}, sorted({r["sacc"] for r in rows})
        for i in range(0, len(accs), 100):
            for attempt in range(3):
                try:
                    with Entrez.efetch(db="nucleotide", id=",".join(accs[i:i + 100]), rettype="gb", retmode="text",
                                       seq_start=1, seq_stop=1) as h:
                        for rec in SeqIO.parse(h, "genbank"):
                            a = rec.annotations
                            refs = a.get("references") or []
                            meta[rec.id.split(".")[0]] = {
                                "organism": a.get("organism", ""), "taxonomy": "; ".join(a.get("taxonomy", [])),
                                "definition": rec.description, "seq_length": "",
                                "topology": a.get("topology", ""), "reference": refs[0].title if refs else ""}
                    break
                except Exception as e:                       # network: retry, then leave the metadata empty
                    print(f"efetch failed ({e}), attempt {attempt + 1}/3", file=sys.stderr, flush=True)
                    time.sleep(10 * (attempt + 1))
        for r in rows:
            r.update(meta.get(r["sacc"].split(".")[0], dict.fromkeys(GB_FIELDS, "")) | {"seq_length": r.get("slen", "")})
            if not r["organism"]:
                r["organism"] = "(efetch failed)"            # visible in the table, not only in the log

    write_tsv(out_tsv, rows, FIELDS + GB_FIELDS)
    print(f"{len(rows)} hits for {len({r['qseqid'] for r in rows})} sequences", file=sys.stderr)
    return rows


def main_snakemake(sm):
    """rule blast: the mitogenome / rDNA queries of every sample (see the module docstring)"""
    sys.stderr = open(sm.log[0], "w")
    out_dir = Path(sm.output.tsv).parent
    query = out_dir / "query.fasta"
    qs = queries("mt", sm.input.get("mt", []), by_sample(sm.input.get("gb", [])), sm.params.mt, mt_genes) + \
        queries("nr", sm.input.get("nr", []), by_sample(sm.input.get("nr_gb", [])), sm.params.nr, nr_regions)
    with open(query, "w") as fh:
        fh.writelines(f">{n}\n{q}\n" for n, q in qs if q)
    blast_topn(query, sm.output.tsv, sm.params.db, sm.params.topn, sm.threads, sm.params.email)


if __name__ == "__main__":
    if "snakemake" in globals():
        main_snakemake(snakemake)  # noqa: F821
    else:
        import argparse
        ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
        ap.add_argument("--query", required=True)
        ap.add_argument("-o", "--out", required=True)
        ap.add_argument("--db", default="")
        ap.add_argument("--topn", type=int, default=5)
        ap.add_argument("--threads", type=int, default=8)
        ap.add_argument("--email", default="")
        a = ap.parse_args()
        blast_topn(a.query, a.out, a.db, a.topn, a.threads, a.email)

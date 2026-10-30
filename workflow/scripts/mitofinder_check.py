"""Are the mitochondrial genes MitoFinder found really the sample's? BLAST every gene of every sample (the
mitofinder_summary.py genes/<gene>.fasta) against nt and compare the best hits: genes of one animal should all
hit the same family; a gene whose best hit is another family and close (a gut content, a parasite, a
cross-contaminated library - common in transcriptomes) is flagged (see check()).

    mitofinder_check.py --genes-dir mitofinder/rna/genes -o mitofinder/rna/blast [--genes all|COX1,CYTB,...]
                        [--db LOCAL_DB] [--topn 5] [--threads 8] [--email E] [--min-pident 85]
Writes, in the output directory:
  query.fasta         the queries, <sample>|<gene> (extra copies <sample>_copy2|<gene>)
  blast_topn.tsv      TopN hits per query with GenBank organism and taxonomy (blast_taxonomy.blast_topn)
  genes_check.tsv     per sample x gene: best hit, identity, coverage, family, verdict
  samples_check.tsv   per sample: majority family / organism, genes agreeing, flagged genes, verdict
"""
import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from blast_taxonomy import blast_topn  # noqa: E402
from seqtools import read_fasta, write_fasta  # noqa: E402
from tables import write_tsv  # noqa: E402


def family(taxonomy):
    """the family in a GenBank lineage (last rank ending in -idae / -aceae), else its last rank"""
    ranks = [x.strip() for x in taxonomy.split(";") if x.strip()]
    fam = [x for x in ranks if x.endswith(("idae", "aceae"))]
    return fam[-1] if fam else (ranks[-1] if ranks else "")


def queries(genes_dir, wanted):
    """[(id, seq)]: <sample>|<gene> for the wanted genes (all = every genes/<gene>.fasta)"""
    files = sorted(Path(genes_dir).glob("*.fasta"))
    keep = None if wanted in ("", "all") else {g.strip().lower() for g in wanted.split(",")}
    out = []
    for f in files:
        gene = f.stem
        if keep is not None and gene.lower() not in keep:
            continue
        out += [(f"{name}|{gene}", seq) for name, seq in read_fasta(f) if seq]
    return out


GENOME_RECORD = 10000          # a hit record this long is a (near-)complete mitogenome: it holds every gene
SAME_SPECIES = 97              # identity of a hit to the same species (the foreign-gene rule applies only then)
OTHERS_MATCH = 90              # another gene of that species matches its mitogenome at least this closely


def foreign_genes(gs, hits, min_pident, min_qcov):
    """queries of one sample that come from another organism, X = a database species the query matches
    near-identically (>= SAME_SPECIES %, over >= min_qcov %):
      two copies: another copy of the same gene in the sample does not match X (>= OTHERS_MATCH %) - two different
                  COX1s are two organisms (sim5 S5, Pisione RNA + Nephtys COX1/CYTB: COX1 copy4 100 % Nephtys,
                  copies 1-3 77-80 % to other families; the two Nephtys genes outvoted the sample's scattered own
                  genes, so by majority alone the sample became 'Nephtyidae')
      complete mitogenome: X's best record is a complete mitogenome (>= GENOME_RECORD bp) that none of the sample's
                  other genes match - with X's mitogenome in nt, X's other genes would match it too (sim5 S5 when
                  MitoFinder kept only the Nephtys COX1)
    and in both cases no other gene of the sample (outside these copies) supports X. Related genera are not
    flagged (S4 rrnS 88.7 % Micronephthys: below SAME_SPECIES). (min_pident: kept for the call signature.)"""
    def gene(q):
        return q.rsplit("|", 1)[1]

    def matches(q, b):
        return any((h["sacc"] == b["sacc"] or h["organism"] == b["organism"]) and float(h["pident"]) >= OTHERS_MATCH
                   for h in hits.get(q, []))
    cand = {}
    for g in gs:
        b = hits.get(g["query"], [None])[0]
        if not b or float(b["pident"]) < SAME_SPECIES or float(b.get("qcovs") or 0) < min_qcov:
            continue
        copies = [o["query"] for o in gs if gene(o["query"]) == gene(g["query"]) and o["query"] != g["query"]
                  and hits.get(o["query"])]
        others = [o["query"] for o in gs if gene(o["query"]) != gene(g["query"]) and hits.get(o["query"])]
        if any(not matches(c, b) for c in copies) or \
                (int(float(b.get("slen") or 0)) >= GENOME_RECORD and others and not any(matches(o, b) for o in others)):
            cand.setdefault(b["organism"], []).append((g["query"], b))
    # candidates do not vouch for each other (the Nephtys COX1 and CYTB of S5 both match the Nephtys mitogenome)
    every = {q for qs in cand.values() for q, _ in qs}
    out = set()
    for org, qs in cand.items():
        if not any(matches(o["query"], b) for o in gs if o["query"] not in every for _, b in qs):
            out |= {q for q, _ in qs}
    return out


def check(rows, qids, min_pident, min_qcov=50):
    """(gene rows, sample rows) from the BLAST hits.

    A gene supports a family only with a hit over >= min_qcov % of it. The sample's family is the one most genes
    support. A gene is flagged as possible contamination only when its best hit is another family AND close
    (>= min_pident): a contaminant matches its source closely, while a gene without a close relative in nt just
    hits whatever is least distant (MJ20260820041 CHY006-008: best COX1 82 % Travisia; ATP6 / COX2 / ND4 best hits
    Enchytraeidae / Drosophila / Bagridae at 72-83 %, often over 15-70 % of the gene - no contamination)."""
    best, hits = {}, {}
    for r in rows:                                    # rows are best first per query
        best.setdefault(r["qseqid"], r)
        hits.setdefault(r["qseqid"], []).append(r)
    genes, by_sample = [], {}
    for q in qids:
        sample, gene = q.rsplit("|", 1)
        base = sample.split("_copy")[0]
        b = best.get(q)
        g = {"sample": base, "query": q, "gene": gene, "best_hit": "", "organism": "", "family": "",
             "pident": "", "qcovs": "", "verdict": "no hit"}
        if b:
            g.update(best_hit=b["sacc"], organism=b["organism"], family=family(b["taxonomy"]),
                     pident=b["pident"], qcovs=b["qcovs"])
        genes.append(g)
        by_sample.setdefault(base, []).append(g)
    samples = []
    for s, gs in by_sample.items():
        foreign = foreign_genes(gs, hits, min_pident, min_qcov)
        hit = [g for g in gs if g["family"] and float(g["qcovs"] or 0) >= min_qcov and g["query"] not in foreign]
        votes = Counter(g["family"] for g in hit)
        major = max(votes, key=lambda f: (votes[f], sum(float(g["pident"]) for g in hit if g["family"] == f)),
                    default="")
        org = Counter(g["organism"] for g in hit if g["family"] == major).most_common(1)
        for g in gs:
            if not g["family"]:
                continue
            pid, cov = float(g["pident"]), float(g["qcovs"] or 0)
            if g["query"] in foreign:
                g["verdict"] = (f"possible contamination: {g['organism']} ({g['family']}, {pid:.1f}%; another copy of "
                                f"the gene or the sample's other genes are not this species)")
            elif cov < min_qcov:
                g["verdict"] = f"partial hit ({cov:.0f}% of the gene)"
            elif g["family"] == major:
                g["verdict"] = "ok" if pid >= min_pident else f"no close reference ({pid:.1f}%)"
            elif pid >= min_pident:
                g["verdict"] = f"possible contamination: {g['organism']} ({g['family']}, {pid:.1f}%)"
            else:
                g["verdict"] = f"no close reference (best {g['family']} {pid:.1f}%)"
        contam = sorted({g["gene"] for g in gs if g["verdict"].startswith("possible contamination")})
        close = [g for g in gs if g["verdict"] == "ok"]
        if not hit:
            verdict = "no hits"
        elif contam:
            verdict = f"possible contamination: {', '.join(contam)}"
        elif close:
            verdict = "consistent"
        else:
            verdict = f"consistent family, no close reference in the database (best {major})"
        samples.append({"sample": s, "family": major, "organism": org[0][0] if org else "",
                        "genes_agree": f"{votes.get(major, 0)}/{len(gs)}", "flagged": ", ".join(
                            f"{g['gene']}: {g['verdict']}" for g in gs if g["verdict"] != "ok"), "verdict": verdict})
    return genes, samples


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-dir", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--genes", default="all")
    ap.add_argument("--db", default="")
    ap.add_argument("--topn", type=int, default=5)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--email", default="")
    ap.add_argument("--min-pident", type=float, default=85)
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    qs = queries(a.genes_dir, a.genes)
    write_fasta(qs, out / "query.fasta", width=0)
    rows = blast_topn(out / "query.fasta", out / "blast_topn.tsv", a.db, a.topn, a.threads, a.email)
    genes, samples = check(rows, [n for n, _ in qs], a.min_pident)
    write_tsv(out / "genes_check.tsv", genes,
              ["sample", "gene", "query", "best_hit", "organism", "family", "pident", "qcovs", "verdict"])
    write_tsv(out / "samples_check.tsv", samples, ["sample", "family", "organism", "genes_agree", "flagged", "verdict"])
    for s in samples:
        print(f"{s['sample']}: {s['verdict']} ({s['genes_agree']} genes {s['family']})")


if __name__ == "__main__":
    main()

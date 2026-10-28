"""Build a small BLAST database for species checks from a local nt: one taxon, records up to a maximum length.

Why: BLAST reads the whole database for every search, and nt (3.4 Tbases, ~1.2 TB) is mostly chromosome-scale
genome records. On MJ20260820041 one search of COX1 + rDNA took 42-52 min whether limited to Metazoa (-taxids) or
not, with 8 or 32 threads: it is disk-bound. Barcode / marker queries (COI, 16S, 12S, 18S, ITS, 28S) are matched by
complete mitogenomes, rDNA records and gene sequences, which are all short; records longer than --max-length
(default 50 kb: genome contigs and chromosomes) are left out. Taxonomy ids are kept (-taxids still works).
nt v5 stores identical sequences once with all their accessions (NR_003027.2, an 18S, carries 920): the listing has
one line per accession, so records are written once per database entry (OID); otherwise makeblastdb stops at
"Duplicate seq_ids" (first build: 619 GB FASTA for 204 Gbases).

Usage: build_marker_db.py --nt /db/blast/nt --out /db/blast_marker/metazoa_marker [--taxid 33208] [--max-length 50000]
Needs blastdbcmd / makeblastdb (BLAST+ >= 2.10) and taxonomy4blast.sqlite3 next to nt.
"""
import argparse
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path


def descendants(sqlite, root):
    """Every taxid below root (inclusive), from taxonomy4blast.sqlite3 (table TaxidInfo: taxid, parent)."""
    con = sqlite3.connect(f"file:{sqlite}?mode=ro", uri=True)
    children = {}
    for t, p in con.execute("select taxid, parent from TaxidInfo"):
        if t != p:
            children.setdefault(p, []).append(t)
    out, todo = set(), [root]
    while todo:
        t = todo.pop()
        out.add(t)
        todo += children.get(t, [])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nt", required=True, help="nt database prefix, e.g. /db/blast/nt")
    ap.add_argument("--out", required=True, help="output database prefix")
    ap.add_argument("--taxid", type=int, default=33208, help="keep this taxon (default 33208 Metazoa)")
    ap.add_argument("--max-length", type=int, default=50000, help="keep records up to this length (default 50000)")
    ap.add_argument("--listing", help="existing listing of nt (blastdbcmd -outfmt '%%a %%l %%T', optionally '%%o' first)")
    ap.add_argument("--fasta", help="reuse an extracted FASTA (duplicate records are dropped)")
    a = ap.parse_args()
    nt, out = Path(a.nt), Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    keep = descendants(nt.parent / "taxonomy4blast.sqlite3", a.taxid)
    print(f"{len(keep):,} taxids below {a.taxid}", file=sys.stderr, flush=True)

    listing = a.listing or f"{out}.nt_listing.txt"
    if not a.listing:
        with open(listing, "w") as fh:
            subprocess.run(["blastdbcmd", "-db", str(nt), "-entry", "all", "-outfmt", "%o %a %l %T"], stdout=fh, check=True)
    acc_file, taxmap = f"{out}.accessions.txt", f"{out}.taxid_map.txt"
    n = bases = 0
    oids = set()
    with open(listing) as fh, open(acc_file, "w") as acc, open(taxmap, "w") as tm:
        for line in fh:
            x = line.split()
            oid, x = (x[0], x[1:]) if len(x) == 4 else (None, x)
            if len(x) == 3 and int(x[1]) <= a.max_length and int(x[2]) in keep:
                tm.write(f"{x[0]} {x[2]}\n")
                if oid is None or oid not in oids:
                    oids.add(oid)
                    acc.write(x[0] + "\n")
                    n += 1
                    bases += int(x[1])
    print(f"{n:,} accessions, {bases / 1e9:.1f} Gbases kept", file=sys.stderr, flush=True)

    fasta = a.fasta or f"{out}.fasta"
    if not a.fasta:
        with open(fasta, "w") as fh:
            subprocess.run(["blastdbcmd", "-db", str(nt), "-entry_batch", acc_file, "-outfmt", "%f"], stdout=fh, check=True)
    mk = subprocess.Popen(["makeblastdb", "-in", "-", "-dbtype", "nucl", "-parse_seqids", "-taxid_map", taxmap,
                           "-out", str(out), "-title", f"nt taxid {a.taxid} <= {a.max_length} bp"], stdin=subprocess.PIPE)
    seen, keep_rec, kept = set(), False, 0
    with open(fasta, "rb") as fh:
        for line in fh:                          # one copy per record (same header = same database entry)
            if line[:1] == b">":
                key = line.split(None, 1)[0]
                keep_rec = key not in seen
                seen.add(key)
                kept += keep_rec
            if keep_rec:
                mk.stdin.write(line)
    mk.stdin.close()
    if mk.wait():
        sys.exit(f"makeblastdb failed ({mk.returncode})")
    print(f"{kept:,} records in the database", file=sys.stderr, flush=True)
    if not a.fasta:
        Path(fasta).unlink()
    for f in ("taxdb.btd", "taxdb.bti", "taxonomy4blast.sqlite3"):     # organism names, -taxids on the new database
        if (nt.parent / f).exists() and not (out.parent / f).exists():
            shutil.copy2(nt.parent / f, out.parent / f)
    print(f"done: {out}", file=sys.stderr)


if __name__ == "__main__":
    main()

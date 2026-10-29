"""MitoFinder results of many samples gathered for analysis (the summary steps of the user's mitofinder-batch,
summary_basic_information / gather_genes / gather_mt_contigs / summarize_assembly_infos, in one pass).

    mitofinder_summary.py -o mitofinder_summary.tsv [--source spades] [--version V] mitofinder/A mitofinder/B ...
Each directory is a rule-mitofinder output: <s>.mitofinder.fasta / .gb and Final_Results/ (MitoFinder's
<s>_final_genes_NT.fasta, *.infos). Written next to the table:
  mitofinder_summary.tsv   sample, contigs, total length, circular, protein genes / rRNAs / tRNAs found, then the
                           length (bp) of every gene found in any sample (N = not found)
  mitofinder_contigs.tsv   one row per mitochondrial contig (.infos): initial contig, length, GC, circularization
  genes/<gene>.fasta       every sample's copy of the gene (<sample>_copy2 ... for extra copies, including the
                           transcripts of the gene MitoFinder did not annotate: <s>.mitofinder.extra_copies.fasta),
                           ready to align
  contigs/<s>_all_mt.fasta every sample's mitochondrial contigs (+ <s>_mt_contigsID_list.txt)
  README.txt               what the files are; MitoFinder version and source
"""
import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from seqtools import read_fasta, write_fasta  # noqa: E402
from tables import write_tsv  # noqa: E402

RRNA = {"rrnl", "rrns", "16s", "12s", "l-rrna", "s-rrna"}
INFOS = {"initial_contig": r"Initial contig name\s*:\s*(.+)", "length": r"Length\s*:\s*([\d.]+)",
         "gc": r"GC content\s*:\s*([\d.]+)%", "circularization": r"Circularization\s*:\s*(.+)"}


def gb_trnas(gb):
    """tRNA features in MitoFinder's GenBank output"""
    if not (gb.exists() and gb.stat().st_size):
        return 0
    return sum(1 for line in open(gb) if re.match(r"^\s{5}tRNA\s", line))


def sample_result(d):
    s = d.name
    res = d / "Final_Results"
    r = {"sample": s, "seqs": [], "genes": {}, "infos": []}
    fasta = d / f"{s}.mitofinder.fasta"
    if fasta.exists():
        r["seqs"] = list(read_fasta(fasta))
    for info in sorted(res.glob("*.infos")):
        text = info.read_text(errors="replace")
        r["infos"].append({k: (m.group(1).strip() if (m := re.search(p, text)) else "") for k, p in INFOS.items()})
    for f in sorted(res.glob("*_final_genes_NT.fasta")):
        for name, seq in read_fasta(f):
            r["genes"].setdefault(name.split("@", 1)[-1].split()[0], []).append(seq)
    r["trna"] = gb_trnas(d / f"{s}.mitofinder.gb")
    extra = d / f"{s}.mitofinder.extra_copies.fasta"        # joined transcripts: copies MitoFinder did not annotate
    r["extra"] = 0
    if extra.exists():
        for name, seq in read_fasta(extra):
            r["genes"].setdefault(name.split()[0], []).append(seq)
            r["extra"] += 1
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--source", default="", help="what MitoFinder searched (spades, rna, trinity, ...)")
    ap.add_argument("--version", default="", help="MitoFinder version, for README.txt")
    a = ap.parse_args()
    rows = [sample_result(Path(d)) for d in a.dirs]
    genes = sorted({g for r in rows for g in r["genes"]}, key=str.lower)
    out = Path(a.out)
    with open(out, "w") as fh:
        fh.write("\t".join(["sample", "contigs", "length", "circular", "protein_genes", "rRNAs", "tRNAs"] + genes)
                 + "\n")
        for r in rows:
            circ = ",".join("yes" if i["circularization"].lower().startswith("yes") else "no" for i in r["infos"])
            rrna = sum(1 for g in r["genes"] if g.lower() in RRNA)
            fh.write("\t".join([r["sample"], str(len(r["seqs"])), str(sum(len(x) for _, x in r["seqs"])), circ,
                                str(len(r["genes"]) - rrna), str(rrna), str(r["trna"])]
                               + [str(max(map(len, r["genes"][g]))) if g in r["genes"] else "N" for g in genes])
                     + "\n")
    write_tsv(out.parent / "mitofinder_contigs.tsv",
              ({"sample": r["sample"], "contig": i, "gc_pct": info["gc"], **info}
               for r in rows for i, info in enumerate(r["infos"], 1)),
              ["sample", "contig", "initial_contig", "length", "gc_pct", "circularization"])
    gdir, cdir = out.parent / "genes", out.parent / "contigs"
    gdir.mkdir(exist_ok=True)
    cdir.mkdir(exist_ok=True)
    for g in genes:
        write_fasta([(f"{r['sample']}{'' if k == 1 else f'_copy{k}'}", seq) for r in rows
                     for k, seq in enumerate(r["genes"].get(g, []), 1)], gdir / f"{g}.fasta", width=0)
    for r in rows:
        if r["seqs"]:
            write_fasta(r["seqs"], cdir / f"{r['sample']}_all_mt.fasta", width=0)
            (cdir / f"{r['sample']}_mt_contigsID_list.txt").write_text(
                "".join(re.sub(r"\(.*?\)", "", n).strip() + "\n" for n, _ in r["seqs"]))
    (out.parent / "README.txt").write_text(
        f"MitoFinder {a.version} on {a.source or 'the assemblies'}: {len(rows)} samples, {len(genes)} genes\n"
        "mitofinder_summary.tsv    per sample: contigs, total length, circular, genes found, length of every gene "
        "(N = not found)\n"
        "mitofinder_contigs.tsv    per contig: MitoFinder .infos (initial contig, length, GC, circularization)\n"
        "genes/<gene>.fasta        the gene of every sample (extra copies <sample>_copyN), ready to align\n"
        "contigs/<s>_all_mt.fasta  the mitochondrial contigs of each sample\n"
        "<s>/Final_Results/        MitoFinder's own output; <s>/<s>.mitofinder.fasta / .gb\n")
    print(f"{len(rows)} samples, {len(genes)} genes -> {out}, {gdir}/, {cdir}/")


if __name__ == "__main__":
    main()

"""One MitoFinder run: an assembly (-a) or reads assembled by MitoFinder itself (MEGAHIT / metaSPAdes), collected as
<out>/<sample>.mitofinder.fasta / .gb (every mitochondrial contig; empty when MitoFinder finds none) and
<out>/Final_Results/ (MitoFinder's own output). Used by rule mitofinder; importable (run_mitofinder) and a command:

    mitofinder_run.py -s SAMPLE -r ref.gb -o OUT (-a assembly.fasta | -1 R1 -2 R2 [--assembler megahit])
                      [--max-pairs 10000000] [--code 5] [--trna mitfi] [--min-contig 1000] [-p 8] [-m 16]
                      [--sif mitofinder_v1.4.2.sif] [--source spades]
Standard library only: started with Snakemake's Python 3 (the bioconda MitoFinder environment is Python 2.7).
MitoFinder 1.4.2 runs from the Singularity image (--sif; `singularity exec --cleanenv`, every input directory
bound); without it the `mitofinder` on PATH (bioconda 1.4.1: circularization not checked, docs/tools/MitoFinder.md).
Reads: MitoFinder's author recommends 5-10 M pairs for the mitogenome (MitoFinder#13), the first --max-pairs are used.
"""
import argparse
import json
import os
import shlex
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))   # the shared modules next to this script
from runtools import run, sh  # noqa: E402
from seqtools import rc, read_fasta, write_fasta  # noqa: E402
from tables import read_tsv, write_tsv  # noqa: E402

READ_ASSEMBLERS = ("megahit", "metaspades")
GENE_ALIAS = {"COI": "COX1", "CO1": "COX1", "COII": "COX2", "CO2": "COX2", "COIII": "COX3", "CO3": "COX3",
              "COB": "CYTB", "CYB": "CYTB", "ATPASE6": "ATP6", "ATPASE8": "ATP8"}


def gene_name(name):
    """reference / MitoFinder gene names on one scheme: COX1 ... ND4L CYTB ATP6 (nad1 -> ND1, coi -> COX1)"""
    n = name.upper().replace("NAD", "ND")
    return GENE_ALIAS.get(n, n)


def first_pairs(r1, r2, wd, max_pairs):
    """the first max_pairs read pairs as wd/R1.fq.gz, wd/R2.fq.gz (0 = all reads, linked)"""
    out = [Path(wd) / "R1.fq.gz", Path(wd) / "R2.fq.gz"]
    for src, dst in zip((r1, r2), out):
        if max_pairs and int(max_pairs) > 0:
            sh(f"(zcat -f {shlex.quote(str(src))} | head -n {int(max_pairs) * 4} || true) | gzip -1 > "
               f"{shlex.quote(str(dst))}")
        else:
            dst.symlink_to(Path(src).resolve())
    return out


def mitofinder_exe(sif, bind):
    """the MitoFinder command: the Singularity image with these directories bound, else `mitofinder`"""
    if not sif:
        return ["mitofinder"]
    binds = [x for d in sorted({str(Path(b).resolve()) for b in bind}) for x in ("-B", d)]
    return ["singularity", "exec", "--cleanenv", *binds, str(sif), "mitofinder"]


def collect(final, out, sample):
    """MitoFinder's *_Final_Results -> out/<sample>.mitofinder.fasta / .gb; returns the number of contigs (0: empty
    files - MitoFinder writes Final_Results even when it finds no mitochondrial contig)"""
    out = Path(out)
    contigs = sorted(f for f in Path(final).glob("*_mtDNA_contig*.fasta") if "_genes_" not in f.name) \
        if final and Path(final).is_dir() else []
    gbs = sorted(Path(final).glob("*_mtDNA_contig*.gb")) if contigs else []
    with open(out / f"{sample}.mitofinder.fasta", "w") as fh:
        fh.writelines(f.read_text() for f in contigs)
    with open(out / f"{sample}.mitofinder.gb", "w") as fh:
        fh.writelines(f.read_text() for f in gbs)
    return len(contigs)


def back_to_transcripts(out, sample, joined, places):
    """genes MitoFinder annotated on the joined transcripts (mitofinder_prefilter.py --concat) -> their transcripts:
    <s>.mitofinder.fasta becomes the transcripts holding a CDS / rRNA, <s>.mitofinder.transcripts.tsv lists every
    feature with its transcript and position there; returns the number of transcripts"""
    out = Path(out)
    gff = next(iter(sorted((out / "Final_Results").glob("*_mtDNA_contig*.gff"))), None)
    if not gff:
        return 0
    seq = read_fasta(joined)[0][1] if Path(joined).stat().st_size else ""
    rows, keep = [], {}
    for line in open(gff):
        t = line.rstrip("\n").split("\t")
        if len(t) < 9 or t[2] not in ("CDS", "rRNA", "tRNA"):
            continue
        a, b = int(t[3]), int(t[4])
        p = next((p for p in places if int(p["start"]) <= a and b <= int(p["end"])), None)
        name = t[8].split("Name=")[-1].split()[0] if "Name=" in t[8] else t[2]
        rows.append({"type": t[2], "name": name, "strand": t[6], "transcript": p["transcript"] if p else "(spans a gap)",
                     "start": a - int(p["start"]) + 1 if p else a, "end": b - int(p["start"]) + 1 if p else b})
        if p and t[2] != "tRNA":
            keep[p["transcript"]] = seq[int(p["start"]) - 1:int(p["end"])]
    write_tsv(out / f"{sample}.mitofinder.transcripts.tsv", rows, ["type", "name", "strand", "transcript", "start",
                                                                   "end"])
    write_fasta(list(keep.items()), out / f"{sample}.mitofinder.fasta", width=0)
    # MitoFinder annotates a gene once per sequence: a second transcript of the same gene (another isoform, or a
    # contaminant's - sim5 S5: Pisione COX1 annotated, the Nephtys COX1 not) would be lost. Candidates whose
    # blastx gene MitoFinder did not annotate on them are kept as extra copies (mitofinder_summary.py: _copyN).
    annotated = {(r["transcript"], gene_name(r["name"])) for r in rows if r["type"] != "tRNA"}
    extra = []
    for p in places:
        if p.get("gene") and (p["transcript"], gene_name(p["gene"])) not in annotated:
            a, b = int(p["gene_start"]), int(p["gene_end"])
            s = seq[int(p["start"]) - 1:int(p["end"])][min(a, b) - 1:max(a, b)]
            extra.append((f"{gene_name(p['gene'])} {p['transcript']} {a}-{b}", rc(s) if a > b else s))
    write_fasta(extra, out / f"{sample}.mitofinder.extra_copies.fasta", width=0)
    return len(keep)


def run_mitofinder(sample, ref, out, assembly=None, r1=None, r2=None, assembler="megahit", max_pairs=10_000_000,
                   code=5, trna="mitfi", min_contig=1000, threads=8, mem_gb=16, sif="", source="", log=None,
                   places=None):
    """MitoFinder on one sample (module doc); returns a status dict, also written to out/mitofinder.json"""
    out = Path(out).resolve()
    wd, res = out / "run", out / "Final_Results"
    log = Path(log).resolve() if log else out / "mitofinder.log"
    shutil.rmtree(wd, ignore_errors=True)
    shutil.rmtree(res, ignore_errors=True)
    wd.mkdir(parents=True)
    ref = Path(ref).resolve()
    status = {"sample": sample, "source": source or (assembler if r1 else "assembly"), "contigs": 0, "note": ""}
    if r1:
        if assembler not in READ_ASSEMBLERS:
            raise ValueError(f"assembler {assembler}: MitoFinder assembles reads with {' or '.join(READ_ASSEMBLERS)}")
        reads = first_pairs(r1, r2, wd, max_pairs)
        data, bind = ["-1", reads[0], "-2", reads[1], f"--{assembler}"], [r1, r2]
    else:
        if not (assembly and Path(assembly).exists() and Path(assembly).stat().st_size):
            log.write_text(f"[wgs_skim] no contig >= {min_contig} bp to search\n")
            collect(None, out, sample)
            shutil.rmtree(wd, ignore_errors=True)
            return save(out, dict(status, note=f"no contig >= {min_contig} bp"))
        shutil.copy(assembly, wd / "assembly.fasta")
        data, bind = ["-a", wd / "assembly.fasta"], [assembly]
    bind = [Path(b).resolve().parent for b in [*bind, ref]] + [wd, out]
    places = read_tsv(places) if places else []
    cmd = mitofinder_exe(sif, bind) + ["-j", sample, *data, "-r", ref, "-o", code, "-p", threads, "-m", int(mem_gb),
                                       "-t", trna, "--min-contig-size", min_contig]
    # joined transcripts: one long sequence (MitoFinder's default --max-contig-size 25000 would drop it), kept in
    # the given orientation so that the positions map back
    cmd += ["--max-contig-size", Path(assembly).stat().st_size + 1000] if places else ["--adjust-direction"]
    # MitoFinder exits 0 after "Aborting." (e.g. bioconda 1.4.1 without install.sh.ok: "MitoFinder is not
    # installed") - that is a failure, not "no mitochondrial contig"
    if run(cmd, log, cwd=wd) or "Aborting" in Path(log).read_text(errors="replace"):
        with open(log, "a") as fh:
            fh.write("[wgs_skim] MitoFinder failed\n")
        raise RuntimeError(f"MitoFinder failed on {sample}; see {log}")
    final = next(iter(sorted((wd / sample).glob("*_Final_Results"))), None)
    if final:
        shutil.move(str(final), str(res))
    status["contigs"] = collect(res if final else None, out, sample)
    if places and status["contigs"]:
        status["transcripts"] = back_to_transcripts(out, sample, wd / "assembly.fasta", places)
    if not status["contigs"]:
        status["note"] = "MitoFinder found no mitochondrial contig"
        with open(log, "a") as fh:
            fh.write(f"[wgs_skim] {status['note']}\n")
    own = next(iter(sorted((wd / sample).glob("*_link_*.scafSeq"))), None)
    if r1 and own:                                           # MitoFinder's own assembly is kept
        shutil.move(str(own.resolve()), str(out / f"{sample}.{status['source']}.contigs.fasta"))
    shutil.rmtree(wd, ignore_errors=True)
    return save(out, status)


def save(out, status):
    (Path(out) / "mitofinder.json").write_text(json.dumps(status, indent=1) + "\n")
    return status


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-s", "--sample", required=True)
    ap.add_argument("-r", "--ref", required=True, help="reference mitogenome(s), GenBank")
    ap.add_argument("-o", "--out", required=True, help="output directory of this sample")
    ap.add_argument("-a", "--assembly")
    ap.add_argument("-1", dest="r1")
    ap.add_argument("-2", dest="r2")
    ap.add_argument("--assembler", default="megahit", choices=READ_ASSEMBLERS)
    ap.add_argument("--max-pairs", type=int, default=10_000_000)
    ap.add_argument("--code", type=int, default=5, help="genetic code (5 = invertebrate mitochondrial)")
    ap.add_argument("--trna", default="mitfi", choices=["mitfi", "trnascan", "arwen"])
    ap.add_argument("--min-contig", type=int, default=1000)
    ap.add_argument("-p", "--threads", type=int, default=8)
    ap.add_argument("-m", "--mem-gb", type=int, default=16)
    ap.add_argument("--sif", default=os.environ.get("MITOFINDER_SIF", ""))
    ap.add_argument("--source", default="")
    ap.add_argument("--log")
    ap.add_argument("--places", help="map.tsv of mitofinder_prefilter.py --concat: -a holds joined transcripts")
    a = ap.parse_args(argv)
    if not (a.assembly or (a.r1 and a.r2)):
        ap.error("give an assembly (-a) or reads (-1 / -2)")
    st = run_mitofinder(a.sample, a.ref, a.out, a.assembly, a.r1, a.r2, a.assembler, a.max_pairs, a.code, a.trna,
                        a.min_contig, a.threads, a.mem_gb, a.sif, a.source, a.log, a.places)
    print(f"{st['sample']}: {st['contigs']} mitochondrial contig(s) {st['note']}".rstrip())


if __name__ == "__main__":
    main()

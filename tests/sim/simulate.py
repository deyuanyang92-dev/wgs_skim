"""Simulated genome-skimming (DNA) and RNA-seq samples with a known truth, for testing every wgs_skim function.

    simulate.py --sheet tests/sim/samples.tsv -o sim5/ [--nuclear-mb 20] [--depth 20] [--rna-nuclear FASTA]
                [--contaminant S5:nephtys_mt.fasta:COX1,CYTB]

Sheet (TSV): sample, mt (GenBank of the mitogenome), mt_topology (circular | linear), rdna (FASTA, one unit),
histone (FASTA of the cluster; the longest record is used, as for rdna), nuclear (assembly FASTA: the nuclear background is cut from it). Per sample:
  dna/<s>_R1.fq.gz, _R2.fq.gz   2 x 150 bp pairs, insert 350 +- 50, 0.2 % substitution errors:
                                nuclear 20 Mb at --depth x, diploid (0.2 % heterozygous SNPs), so SPAdes, BUSCO and
                                GenomeScope have a genome; mitogenome at --mt-depth x (circular: reads span the
                                junction), rDNA unit --rdna-copies x --depth, histone cluster --histone-copies x
                                --depth (tandem units: reads span unit junctions)
  rna/<s>_R1.fq.gz, _R2.fq.gz   transcripts: every mitochondrial gene of the sample (CDS, rRNA) at --mt-rna-depth x,
                                nuclear transcripts (--rna-nuclear, the same for every sample) at log-normal depth;
                                --contaminant adds genes of another animal's mitogenome (the regions matching the
                                sample's own genes, blastn) - the positive control of the MitoFinder contamination check
  truth/<s>/                    mt.fasta, rdna.fasta, histone.fasta, nuclear.fasta, mt_genes.fasta, contaminant.fasta
  truth.tsv                     lengths, topology, depths, read pairs, sources
Reads are drawn from all parts in random order (GetOrganelle --max-reads and MitoFinder use the first reads).
Standard library only (+ blastn on PATH for --contaminant). Deterministic (--seed)."""
import argparse
import gzip
import math
import random
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow" / "scripts"))
from seqtools import genbank_seqs, iter_fasta, rc, read_fasta, write_fasta  # noqa: E402
from tables import write_tsv  # noqa: E402

READ, INSERT, INSERT_SD, ERR = 150, 350, 50, 0.002
COV = re.compile(r"_cov_([\d.]+)")


# ── truth sequences ────────────────────────────────────────────────────────
def genbank_genes(path):
    """[(gene, SEQUENCE)] of the CDS and rRNA features of a GenBank file (join / complement locations)"""
    seq = genbank_seqs(path)[0][1]
    genes, feat, loc, key = [], None, "", None
    lines = open(path).read().split("ORIGIN")[0].splitlines()
    for line in lines + ["     end"]:
        if len(line) > 5 and line[5] != " " and line.startswith("     "):        # a new feature
            if feat in ("CDS", "rRNA") and key:
                genes.append((key, cut(seq, loc)))
            parts = line.split()
            feat, loc, key = parts[0], parts[1] if len(parts) > 1 else "", None
        elif line.startswith("                     /"):
            m = re.match(r'\s+/(gene|product)="([^"]+)"', line)
            if m and key is None:
                key = m.group(2).split()[0]
        elif line.startswith("                     ") and feat and not line.strip().startswith("/") and key is None:
            loc += line.strip()
    return [(g, s) for g, s in genes if s]


def cut(seq, loc):
    """the sequence of a GenBank location (join(...), complement(...), <, >)"""
    loc = loc.replace("<", "").replace(">", "")
    comp = loc.startswith("complement(")
    parts = [tuple(int(x) for x in p.split("..")) for p in re.findall(r"\d+\.\.\d+", loc)]
    s = "".join(seq[a - 1:b] for a, b in parts)
    return rc(s) if comp else s


def nuclear_background(asm, mb):
    """the longest scaffolds of an assembly up to `mb` Mb, at nuclear depth (0.5-2 x the length-weighted median
    k-mer coverage: no mitochondrial / rDNA / repeat scaffolds); [(name, SEQ)]"""
    stats = [(n, len(s), float(m.group(1)) if (m := COV.search(n)) else 1.0)
             for n, s in iter_fasta(asm, upper=False) if len(s) >= 5000]
    total, acc, med = sum(x[1] for x in stats), 0, 1.0
    for _, ln, c in sorted(stats, key=lambda x: x[2]):
        acc += ln
        if acc >= total / 2:
            med = c
            break
    keep, size = set(), 0
    for n, ln, c in sorted(stats, key=lambda x: -x[1]):
        if 0.5 * med <= c <= 2 * med and size < mb * 1e6:
            keep.add(n)
            size += ln
    return [(n.split()[0], s.replace("N", "")) for n, s in iter_fasta(asm) if n in keep]


def homologous(query_genes, subject_fasta, wanted):
    """the regions of subject_fasta matching the wanted genes of the sample (blastn dc-megablast, best HSP)"""
    out = []
    subj = read_fasta(subject_fasta)[0][1]
    for g, s in query_genes:
        if g.upper() not in wanted:
            continue
        with tempfile.NamedTemporaryFile("w", suffix=".fa") as q:
            q.write(f">{g}\n{s}\n")
            q.flush()
            r = subprocess.run(["blastn", "-task", "dc-megablast", "-query", q.name, "-subject", subject_fasta,
                                "-outfmt", "6 sstart send bitscore", "-max_hsps", "1"], capture_output=True,
                               text=True, check=True).stdout.split()
        if r:
            a, b = sorted((int(r[0]), int(r[1])))
            region = subj[a - 1:b]
            out.append((f"contaminant_{g}", region if int(r[0]) < int(r[1]) else rc(region)))
    return out


# ── reads ──────────────────────────────────────────────────────────────────
class Template:
    """a molecule reads are drawn from: circular (start anywhere, wraps), tandem (a circular unit) or linear;
    diploid: a second haplotype with SNPs, each fragment from one of the two"""

    def __init__(self, name, seq, circular, rng, het=0.0):
        self.name, self.circular = name, circular
        self.haps = [seq]
        if het:
            s = list(seq)
            for i in range(len(s)):
                if rng.random() < het:
                    s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
            self.haps.append("".join(s))
        self.length = len(seq)

    def fragment(self, rng, insert):
        h = self.haps[rng.randrange(len(self.haps))]
        n = min(insert, self.length)
        if self.circular:
            i = rng.randrange(self.length)
            f = (h + h[:n])[i:i + n] if i + n > self.length else h[i:i + n]
        else:
            i = rng.randrange(max(self.length - n, 0) + 1)
            f = h[i:i + n]
        return f if rng.random() < 0.5 else rc(f)


def poisson(lam, rng):
    k, p, limit = 0, rng.random(), math.exp(-lam)
    while p > limit:
        k += 1
        p *= rng.random()
    return k


def errors(s, rng):
    """substitution errors: Poisson(len x ERR) random positions"""
    k = poisson(len(s) * ERR, rng) if s else 0
    if not k:
        return s
    s = list(s)
    for _ in range(k):
        i = rng.randrange(len(s))
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def write_reads(parts, prefix, sample, rng):
    """parts: [(Template, pairs)]; pairs drawn in random order -> prefix_R1.fq.gz / _R2.fq.gz; returns pairs"""
    order = [j for j, (_, p) in enumerate(parts) for _ in range(p)]
    rng.shuffle(order)
    total = len(order)
    q = "F" * READ
    with gzip.open(f"{prefix}_R1.fq.gz", "wt", compresslevel=1) as f1, \
            gzip.open(f"{prefix}_R2.fq.gz", "wt", compresslevel=1) as f2:
        for k, j in enumerate(order):
            frag = parts[j][0].fragment(rng, max(READ + 20, int(rng.gauss(INSERT, INSERT_SD))))
            r1, r2 = errors(frag[:READ], rng), errors(rc(frag)[:READ], rng)
            f1.write(f"@{sample}:{k + 1} 1:N:0\n{r1}\n+\n{q[:len(r1)]}\n")
            f2.write(f"@{sample}:{k + 1} 2:N:0\n{r2}\n+\n{q[:len(r2)]}\n")
    return total


def pairs_for(length, depth):
    return max(1, int(round(length * depth / (2 * READ))))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--nuclear-mb", type=float, default=20)
    ap.add_argument("--depth", type=float, default=20, help="nuclear depth (x)")
    ap.add_argument("--het", type=float, default=0.002, help="heterozygous SNP rate of the nuclear genome")
    ap.add_argument("--mt-depth", type=float, default=300)
    ap.add_argument("--rdna-copies", type=float, default=100)
    ap.add_argument("--histone-copies", type=float, default=50)
    ap.add_argument("--rna-nuclear", help="FASTA of nuclear transcripts (the first --rna-transcripts >= 500 bp used)")
    ap.add_argument("--rna-transcripts", type=int, default=1500)
    ap.add_argument("--mt-rna-depth", type=float, default=500)
    ap.add_argument("--contaminant", action="append", default=[], metavar="SAMPLE:MT.fasta:GENE,GENE")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    out = Path(a.out)
    for d in ("dna", "rna", "truth"):
        (out / d).mkdir(parents=True, exist_ok=True)
    rows = [dict(zip(["sample", "mt", "mt_topology", "rdna", "histone", "nuclear"], line.rstrip("\n").split("\t")))
            for line in open(a.sheet) if line.strip() and not line.startswith(("#", "sample\t"))]
    contam = {c.split(":")[0]: (c.split(":")[1], {g.upper() for g in c.split(":")[2].split(",")})
              for c in a.contaminant}
    nuc_rna = []
    if a.rna_nuclear:
        nuc_rna = [(n.split()[0], s) for n, s in iter_fasta(a.rna_nuclear) if len(s) >= 500][:a.rna_transcripts]
    truth = []
    for i, r in enumerate(rows):
        s = r["sample"]
        rng = random.Random(a.seed * 1000 + i)
        t = out / "truth" / s
        t.mkdir(exist_ok=True)
        mt = genbank_seqs(r["mt"])[0][1]
        rdna = max((x for _, x in read_fasta(r["rdna"])), key=len)          # the longest record: one whole unit
        hist = max((x for _, x in read_fasta(r["histone"])), key=len)
        print(f"[sim] {s}: nuclear background from {r['nuclear']}", flush=True)
        nuc = nuclear_background(r["nuclear"], a.nuclear_mb)
        genes = genbank_genes(r["mt"])
        write_fasta([(f"{s}_mt {r['mt_topology']}", mt)], t / "mt.fasta")
        write_fasta([(f"{s}_rdna_unit", rdna)], t / "rdna.fasta")
        write_fasta([(f"{s}_histone_unit", hist)], t / "histone.fasta")
        write_fasta(nuc, t / "nuclear.fasta")
        write_fasta(genes, t / "mt_genes.fasta")
        nuc_len = sum(len(x) for _, x in nuc)
        dna = [(Template(n, x, False, rng, a.het), pairs_for(len(x), a.depth)) for n, x in nuc]
        dna += [(Template("mt", mt, r["mt_topology"] == "circular", rng), pairs_for(len(mt), a.mt_depth)),
                (Template("rdna", rdna, True, rng), pairs_for(len(rdna), a.rdna_copies * a.depth)),
                (Template("histone", hist, True, rng), pairs_for(len(hist), a.histone_copies * a.depth))]
        print(f"[sim] {s}: DNA {sum(p for _, p in dna):,} pairs", flush=True)
        n_dna = write_reads(dna, out / "dna" / s, s, rng)
        rna = [(Template(g, x, False, rng), pairs_for(len(x), a.mt_rna_depth)) for g, x in genes]
        rna += [(Template(n, x, False, rng), pairs_for(len(x), min(500, math.exp(rng.gauss(math.log(20), 1)))))
                for n, x in nuc_rna]
        cont = homologous(genes, contam[s][0], contam[s][1]) if s in contam else []
        if cont:
            write_fasta(cont, t / "contaminant.fasta")
            rna += [(Template(n, x, False, rng), pairs_for(len(x), a.mt_rna_depth * 2)) for n, x in cont]
        n_rna = write_reads(rna, out / "rna" / s, s, rng)
        truth.append(dict(sample=s, mt_length=len(mt), mt_topology=r["mt_topology"], mt_genes=len(genes),
                          rdna_unit=len(rdna), histone_unit=len(hist), nuclear_mb=f"{nuc_len / 1e6:.1f}",
                          nuclear_depth=a.depth, het=a.het, mt_depth=a.mt_depth,
                          rdna_depth=a.rdna_copies * a.depth, histone_depth=a.histone_copies * a.depth,
                          dna_pairs=n_dna, rna_pairs=n_rna, rna_mt_genes=len(genes),
                          rna_contaminant=",".join(n for n, _ in cont), mt_source=r["mt"], rdna_source=r["rdna"],
                          histone_source=r["histone"], nuclear_source=r["nuclear"]))
        write_tsv(out / "truth.tsv", truth, list(truth[0]))
        print(f"[sim] {s}: done ({n_dna:,} DNA, {n_rna:,} RNA pairs)", flush=True)


if __name__ == "__main__":
    main()

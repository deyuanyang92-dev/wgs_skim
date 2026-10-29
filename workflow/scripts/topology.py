"""Is a mitogenome / rDNA sequence really circular or linear? Read-based check of the junction (rule topology).

A. End overlap: k-mers of the first 2 kb that recur in the last 5 kb give the circle length P; the longest stretch with
   s[i] == s[i+P] (>= 25 bp, ragged ends tolerated) confirms it. No overlap: the sequence is joined end to end as is.
B. Junction: the circle is cut inside nad5 (else the longest CDS, else the middle) so that the original ends meet in
   the middle of the reference; reads spanning the junction (>= 20 bp each side, NM <= 5%) are compared with spanning
   reads 150-600 bp either side (Ahuja et al. 2024, GBE 16:evae048). Validated on MJ20260820041: true circles 0.56-1.02,
   SPAdes circles with a collapsed control-region repeat 0.02, truncated fragments 0.
C. If the ends do not join: read pairs leaving the end outwards with the mate leaving the start outwards = a circle with
   a gap (gap length estimated); none = linear_unjoined (true linear, or a gap longer than the insert size).

verdicts: circular_validated, circular_weak (0.2-0.5 x control), circular_repeat_junction (> 2 x: junction inside a
repeat, length uncertain), circular_rejected (overlap but no spanning reads), linear_gap_circular, linear_unjoined,
*_low_reads (control < 5 reads).
params.mode: auto = only linear mitogenomes are checked (circular ones are listed as not read-checked); yes = every
sequence, the rDNA too (its "circle" is the tandem repeat unit). Short reads can show that a junction is real, but
never prove a linear molecule: linear_unjoined may be a gap longer than the insert size (long reads / PCR needed).
"""
import json
import re
import shutil
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))   # the shared modules next to this script
if "snakemake" in globals():                                  # a Snakemake script: run from a temporary copy
    sys.path.insert(0, snakemake.scriptdir)  # noqa: F821
from runtools import sh  # noqa: E402
from seqtools import find_circle, rc, read_fasta, write_fasta  # noqa: E402
from tables import write_tsv  # noqa: E402

FLANK, PAIRS, THREADS = 20, 0, 4
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


# ── annotation: where to cut (nad5) ───────────────────────────────────────────
def genbank_cds(path):
    """[(sequence, [(gene, start0, end0)])] for the CDS of every record (minimal parser, no Biopython)."""
    recs, seq, feats, cur, in_seq = [], [], [], None, False
    for line in open(path):
        if line.startswith("LOCUS"):
            seq, feats, cur, in_seq = [], [], None, False
        elif line.startswith("ORIGIN"):
            in_seq = True
        elif line.startswith("//"):
            recs.append(("".join(seq).upper(), [f for f in feats if f[0]]))
        elif in_seq:
            seq.append(re.sub(r"[^A-Za-z]", "", line))
        elif line.startswith("     CDS "):
            iv = [(int(a) - 1, int(b)) for a, b in re.findall(r"<?(\d+)\.\.>?(\d+)", line)]
            cur = ["", iv[0][0], iv[-1][1]] if iv and iv[0][0] < iv[-1][1] else None
            if cur:
                feats.append(cur)
        elif cur is not None and "/gene=" in line and not cur[0]:
            cur[0] = line.split('"')[1].lower().replace("nd5", "nad5")
        elif not line.startswith(" " * 21):
            cur = None
    return recs


def split_point(s, ann):
    """Middle of nad5 (else the longest CDS) in the coordinates of s; the annotated record may be rotated / reversed."""
    L = len(s)
    for q, cds in ann:
        if len(q) != L or not cds:
            continue
        for strand, t in ((1, s), (-1, rc(s))):
            o = (t + t).find(q)
            if 0 <= o < L:
                g = next((c for c in cds if c[0] == "nad5"), None) or max(cds, key=lambda c: c[2] - c[1])
                m = ((g[1] + g[2]) // 2 + o) % L
                return g[0], (m if strand > 0 else L - 1 - m)
    return None


# ── reads ─────────────────────────────────────────────────────────────────────
def alignments(sam):
    for line in open(sam):
        if line[0] == "@":
            continue
        f = line.split("\t")
        if f[5] == "*":
            continue
        ops = [(int(n), op) for n, op in CIG.findall(f[5])]
        nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), 0)
        mc = next((t[5:].strip() for t in f[11:] if t.startswith("MC:Z:")), "")
        yield int(f[1]), int(f[3]) - 1, ops, nm, f[6], int(f[7]) - 1, mc, f[9]


def ref_span(ops):
    return sum(n for n, op in ops if op in "MDN=X")


def spanning(sam, points):
    """Reads whose alignment covers [p - FLANK, p + FLANK] with NM <= 5% (read level, so homopolymer indels count)."""
    cnt = dict.fromkeys(points, 0)
    for _, start, ops, nm, *_ in alignments(sam):
        end = start + ref_span(ops)
        if nm <= 0.05 * max(end - start, 1):
            for p in points:
                cnt[p] += start <= p - FLANK and end >= p + FLANK
    return cnt


def bridging_pairs(sam, L, window=600):
    """Pairs with the read leaving the end outwards and the mate leaving the start outwards, and the gap they imply
    (their assembled span vs the fragment lengths of pairs chosen the same way at internal points)."""
    pts = [int(L * f) for f in (0.25, 0.5, 0.75)] if L > 4 * window else []
    bridge, ctrl = [], []
    for flag, pos, ops, _, mref, mpos, mc, _ in alignments(sam):
        if not flag & 1 or flag & 0x91C or mref != "=" or not flag & 0x20 or not mc:
            continue                                # proper orientation: read forward, mate reverse
        mspan = ref_span([(int(n), op) for n, op in CIG.findall(mc)])
        if pos >= L - window and mpos < window:
            bridge.append(L - pos + mpos + mspan)
        ctrl += [mpos + mspan - pos for p in pts if p - window <= pos < p <= mpos < p + window]
    gap = None
    if len(bridge) >= 3 and len(ctrl) >= 20:        # the gap removes short fragments: fit, don't subtract medians
        obs, lo, best = statistics.median(bridge), min(bridge), None
        for g in range(-300, 3001, 5):
            pred = [x - g for x in ctrl if x - g >= lo]
            if len(pred) < 5:
                break
            d = abs(statistics.median(pred) - obs)
            if best is None or d < best[0]:
                best = (d, g)
        gap = best[1] if best else None
    return len(bridge), gap


def end_reads(sam, L, window=30):
    """Per end: reads reaching it and the parts of them soft-clipped outwards (>= 10 bp) = sequence beyond the end."""
    ends = {"L": [0, []], "R": [0, []]}
    for _, start, ops, *_, seq in alignments(sam):
        if start < window:
            ends["L"][0] += 1
            if ops[0][1] == "S" and ops[0][0] >= 10:
                ends["L"][1].append(seq[:ops[0][0]])
        if start + ref_span(ops) > L - window:
            ends["R"][0] += 1
            if ops[-1][1] == "S" and ops[-1][0] >= 10:
                ends["R"][1].append(seq[-ops[-1][0]:])
    return ends


def end_diagnosis(side, n, tails, terminal, med):
    """Why the assembly stops at this end, from the reads that run past it (clipped parts, oriented outwards).
    MJ20260820041 HN_Pisione_1: left end - 38/70 reads continue, (AT)n in many of them and the continuations
    disagree (the repeat occurs at several loci, so the assembler cannot choose); right end - poly-G, half the
    usual coverage."""
    notes = []
    tails = [t[::-1] for t in tails] if side == "L" else tails      # base 0 = next to the end
    if tails and len(tails) >= 0.3 * n:
        rep = sum(bool(re.search(r"(AT|TA|AC|CA|AG|GA|TG|GT|TC|CT|CG|GC)\1{4,}", t)) for t in tails) / len(tails)
        cols = [[t[i] for t in tails if len(t) > i] for i in range(20)]
        agree = statistics.mean(max(c.count(b) for b in "ACGT") / len(c) for c in cols if len(c) >= 3) \
            if any(len(c) >= 3 for c in cols) else 0
        what = f"{len(tails)}/{n} reads continue"
        if rep >= 0.3:
            what += f", {rep:.0%} into a dinucleotide repeat"
        what += (", consistently: extendable, unassembled" if agree >= 0.8 else
                 f", in different directions ({agree:.0%} agreement): repeat shared with other loci")
        notes.append(what)
    if med and n < 0.6 * med:
        notes.append(f"coverage drops ({n} reads vs {med:.0f})")
    run = max((len(m.group()) for m in re.finditer(r"A+|C+|G+|T+", terminal)), default=0)
    if run >= 7:
        notes.append(f"homopolymer at the end ({run} bp)")
    return f"{side}: " + ("; ".join(notes) if notes else "clean end (possible true terminus)")


def map_reads(ref_seq, name, wd, reads):
    ref = f"{wd}/{name}.fasta"
    write_fasta([(name, ref_seq)], ref, width=0)
    sh(f"bwa index {ref} 2>/dev/null && bwa mem -t {THREADS} {ref} {reads} 2>/dev/null "
       f"| samtools view -F 0x904 - > {wd}/{name}.sam")
    return f"{wd}/{name}.sam"


# ── main ──────────────────────────────────────────────────────────────────────
def main(sm):
    global PAIRS, THREADS
    sys.stderr = open(sm.log[0], "w")
    PAIRS, THREADS = sm.params.pairs, sm.threads
    # auto: GetOrganelle closed the circle itself - no reads needed; a SPAdes scaffold closed by circularize (graph
    # walk / end overlap / join) is checked, its junction can sit in a collapsed control-region repeat
    status = json.load(open(sm.input.status)) if sm.input.get("status") else {}
    assembler_circle = status.get("source", "getorganelle") == "getorganelle"
    out_dir = Path(sm.output.tsv).parent
    wd = out_dir / "work"
    wd.mkdir(parents=True, exist_ok=True)
    recs, skipped = [], []
    for f in sm.input.fasta:
        for n, q in read_fasta(f):
            # auto: reads only for linear mitogenomes; circular ones were closed by the assembler / circularize step
            if sm.params.mode == "auto" and ("|mt|" not in n or (" circular" in n and assembler_circle)):
                why = "circular: not read-checked" if " circular" in n else "rDNA: not checked"
                skipped.append((n.split()[0], len(q), why))
            else:
                recs.append((n.split()[0], q))
    ann = genbank_cds(sm.input.annotation) if sm.input.get("annotation") else []

    rows = []
    if recs:
        # 1. keep only the read pairs that map to any of the sequences (fast: one bwa run on the first PAIRS pairs)
        write_fasta([(f"s{i}", s) for i, (_, s) in enumerate(recs)], f"{wd}/all.fasta", width=0)
        head = f"| head -n {PAIRS * 4}" if PAIRS else ""
        sh(f"bwa index {wd}/all.fasta 2>/dev/null && bwa mem -t {sm.threads} {wd}/all.fasta "
           f"<(zcat -f {sm.input.r1} {head}) <(zcat -f {sm.input.r2} {head}) 2> {wd}/mem.log "
           f"| samtools view -u -G 12 - | samtools collate -O -u - {wd}/col "
           f"| samtools fastq -1 {wd}/s1.fq -2 {wd}/s2.fq -0 /dev/null -s /dev/null - 2>/dev/null")
        if not re.search(r"\[M::process\] read [1-9]", Path(f"{wd}/mem.log").read_text()):
            sys.exit(f"no reads could be read from {sm.input.r1} / {sm.input.r2}")
        reads = f"{wd}/s1.fq {wd}/s2.fq"

        for i, (name, s) in enumerate(recs):
            # 2. each sequence on its own index (similar sequences would share and dilute the reads)
            c = find_circle(s)
            circle = s[c["anchor"]:c["anchor"] + c["P"]] if c else s
            sp = split_point(s, ann)
            h = ((sp[1] - c["anchor"]) % c["P"] if c else sp[1]) if sp else len(circle) // 2
            h = h if 0 < h < len(circle) else len(circle) // 2
            rot, J = circle[h:] + circle[:h], len(circle) - h
            local = [J + d for d in (-600, -300, -150, 150, 300, 600) if FLANK < J + d < len(rot) - FLANK]
            cnt = spanning(map_reads(rot, f"rot{i}", wd, reads), [J] + local)
            junction, med = cnt[J], statistics.median([cnt[p] for p in local]) if local else 0
            lin = map_reads(s, f"lin{i}", wd, reads)
            bridge, gap = bridging_pairs(lin, len(s))
            ends = end_reads(lin, len(s))

            if med < 5:
                verdict = ("circular" if c else "linear") + "_low_reads"
            elif junction > 2 * med:
                verdict = "circular_repeat_junction"
            elif junction >= 0.5 * med:
                verdict = "circular_validated"
            elif junction >= 0.2 * med:
                verdict = "circular_weak"
            elif c:
                verdict = "circular_rejected"
            elif bridge >= max(5, 0.05 * med):
                verdict = "linear_gap_circular"
            else:
                verdict = "linear_unjoined"
            rows.append(dict(
                sample=sm.wildcards.sample, contig=name, length=len(s), circle_length=len(circle),
                overlap_bp=c["overlap"] if c else 0, split_at=f"{sp[0]}:{sp[1] + 1}" if sp else "midpoint",
                junction_reads=junction, control_reads_median=med, bridging_pairs=bridge,
                gap_estimate_bp="" if gap is None else gap,
                end_clip=f"L {len(ends['L'][1])}/{ends['L'][0]} R {len(ends['R'][1])}/{ends['R'][0]}", verdict=verdict,
                end_diagnosis="" if verdict.startswith("circular") else " | ".join(
                    end_diagnosis(e, ends[e][0], ends[e][1], s[:30] if e == "L" else s[-30:], med) for e in "LR")))
            print(rows[-1], file=sys.stderr, flush=True)
            if verdict in ("circular_validated", "circular_repeat_junction"):
                write_fasta([(name, circle)], out_dir / f"{name.replace('|', '_')}.circular.fasta", width=0)

    rows += [dict(sample=sm.wildcards.sample, contig=n, length=L, verdict=v) for n, L, v in skipped]
    fields = ["sample", "contig", "length", "circle_length", "overlap_bp", "split_at", "junction_reads",
              "control_reads_median", "bridging_pairs", "gap_estimate_bp", "end_clip", "verdict", "end_diagnosis"]
    write_tsv(sm.output.tsv, rows, fields)
    shutil.rmtree(wd, ignore_errors=True)


if __name__ == "__main__" and "snakemake" in globals():
    main(snakemake)  # noqa: F821

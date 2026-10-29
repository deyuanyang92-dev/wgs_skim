"""Reseed GetOrganelle until the mitogenome is circular (Snakemake script for rule mt_iterate).

Round 0 is the first assembly (rule getorganelle_mt). If it is not circular, each further round uses the previous
round's result as the seed (-s), a higher -R and all clean reads (a round with no sequence to reseed from is
tried once with all reads, not repeated); when the clean data exceed `max_gb` Gb only that
much is used (converted to GetOrganelle --max-reads with the mean read length from fastp). Stops at the first
circular round, and when a round returns the same sequence as the round before: more rounds cannot help then
(MJ20260820041: CHY007, HN_Pisione_1/2, WenSC001 changed by 0-74 bp after round 1 and never again; ~6 h of identical
rounds). The result is the first circular round, otherwise the last round that gave a result.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))   # the shared modules next to this script
if "snakemake" in globals():                                  # a Snakemake script: run from a temporary copy
    sys.path.insert(0, snakemake.scriptdir)  # noqa: F821
from collect_fasta import best, clean, getorganelle_cmd, renamed  # noqa: E402
from runtools import run  # noqa: E402
from seqtools import write_fasta  # noqa: E402


def main(sm):
    """the rule's work (importing this module runs nothing)"""
    sample = sm.wildcards.sample
    p = sm.params
    base = Path(sm.output.fasta).parent
    log = open(sm.log[0], "w")


    def max_reads():
        """Reads per file for the data cap, or None for all reads."""
        if not p.max_gb or float(p.max_gb) <= 0:
            return None
        s = json.load(open(sm.input.json))["summary"]["after_filtering"]
        if s["total_bases"] <= float(p.max_gb) * 1e9:
            return None
        mean_len = (s["read1_mean_length"] + s.get("read2_mean_length", s["read1_mean_length"])) / 2
        return int(float(p.max_gb) * 1e9 / 2 / mean_len)


    rounds = [("round0", *best(base / "round0"))]
    clean(base / "round0")
    log.write(f"round0: {len(rounds[0][1])} sequence(s), circular={rounds[0][2]}\n")
    cap = max_reads()
    stop = "circular" if rounds[0][2] else "rounds exhausted"
    for k in range(1, int(p.rounds) + 1):
        if any(c for _, _, c in rounds):
            stop = "circular"
            break
        if k >= 2 and rounds[-1][1] and sorted(s for _, s in rounds[-1][1]) == sorted(s for _, s in rounds[-2][1]):
            stop = f"no change after {rounds[-1][0]}"
            break
        if k >= 2 and not any(r for _, r, _ in rounds):
            # no sequence yet: no seed to reseed with, so this round would repeat the last one (same seeds, reads,
            # -R) - sim5 with GetOrganelle's own seeds: rounds 1-3 identical and empty, 75 min for nothing
            stop = f"no sequence in round0-{rounds[-1][0]} (seeds too distant: use MEANGS, --steps mito)"
            break
        prev = next((r for _, r, _ in reversed(rounds) if r), None)
        out = base / f"round{k}"
        seed = base / f"round{k}.seed.fasta" if prev else None
        if prev:
            write_fasta(prev, seed)
        cmd = getorganelle_cmd(sm.input.r1, sm.input.r2, out, "animal_mt", sm.threads, seed, p.R, p.k, cap,
                               sm.input.db)
        log.write(f"round{k}: {' '.join(cmd)}\n")
        rc = run(cmd, log)
        recs, circ = best(out)
        log.write(f"round{k}: rc={rc}, {len(recs)} sequence(s), circular={circ}, {clean(out) / 1e6:.0f} MB of "
                  "intermediates removed\n")
        rounds.append((f"round{k}", recs, circ))

    chosen = next((r for r in rounds if r[2]), None) or next((r for r in reversed(rounds) if r[1]), None) \
        or ("none", [], False)
    name, recs, circ = chosen
    write_fasta(renamed(recs, sample, "mt"), sm.output.fasta)
    json.dump({"sample": sample, "circular": circ, "round": name, "rounds_run": len(rounds) - 1,
               "n_sequences": len(recs), "length": sum(len(s) for _, s in recs),
               "max_reads_per_file": cap, "stop": stop}, open(sm.output.status, "w"), indent=1)
    log.write(f"final: {name} circular={circ} ({stop})\n")


if __name__ == "__main__" and "snakemake" in globals():
    main(snakemake)  # noqa: F821

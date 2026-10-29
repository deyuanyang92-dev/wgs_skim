"""Final mitogenome: GetOrganelle, or the SPAdes whole-genome scaffold when that one is better (rule mt_select).

GetOrganelle extends from seeds by recruiting reads and stops at regions it cannot extend through (AT-rich control
region repeats, coverage drops); a whole-genome assembly has no recruitment step and sometimes closes those
(MJ20260820041 CHY007: GetOrganelle 16,392 bp linear, SPAdes scaffold 17,040 bp with an end overlap, 100% identical
over the GetOrganelle sequence). Rules, in order:
  1. GetOrganelle circular                                   -> GetOrganelle
  2. a SPAdes scaffold contains the GetOrganelle sequence (>= 98% identity, >= 90% of it) and has an end overlap
     -> that circle (if not shorter than GetOrganelle)
  3. otherwise                                               -> GetOrganelle; a longer linear scaffold is only
     reported: its extra ends are not reliable (HNP1sub: +25 bp of (TG)n at 983x depth = collapsed repeat,
     +92 bp at 3.8x vs 50x = unsupported)
Without a GetOrganelle result, scaffolds matching the MEANGS seed (>= 95%) are used the same way, linear or not.
"""
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))   # the shared modules next to this script
if "snakemake" in globals():                                  # a Snakemake script: run from a temporary copy
    sys.path.insert(0, snakemake.scriptdir)  # noqa: F821
from circularize import circularize  # noqa: E402
from seqtools import read_fasta, write_fasta  # noqa: E402


def main(sm):
    """the rule's work (importing this module runs nothing)"""
    sys.stderr = open(sm.log[0], "w")
    sample = sm.wildcards.sample
    status = json.load(open(sm.input.status))
    go = [(n, s) for n, s in read_fasta(sm.input.fasta)]
    go_len = sum(len(s) for _, s in go)
    final, source, note = go, "getorganelle", ""


    def blast_hits(query, scaffolds, min_pid):
        """{scaffold: aligned query bp} for HSPs >= min_pid (a BLAST database: genomes are too big for -subject)."""
        db = Path(sm.output.status).parent / "scaffolds_db"
        subprocess.run(["makeblastdb", "-in", scaffolds, "-dbtype", "nucl", "-out", str(db)], check=True,
                       stdout=sys.stderr, stderr=sys.stderr)
        out = subprocess.run(["blastn", "-query", query, "-db", str(db), "-evalue", "1e-20", "-num_threads",
                              str(sm.threads), "-outfmt", "6 sseqid pident length"], capture_output=True, text=True,
                             check=True).stdout
        for f in db.parent.glob("scaffolds_db.*"):
            f.unlink()
        cov = {}
        for line in out.splitlines():
            sid, pid, length = line.split("\t")
            if float(pid) >= min_pid:
                cov[sid] = cov.get(sid, 0) + int(length)
        return cov


    def scaffold_seqs(path, ids):
        """Only the wanted scaffolds, streamed (the genome assembly does not need to fit in memory)."""
        seqs, cur = {}, None
        for line in open(path):
            if line.startswith(">"):
                cur = line[1:].split()[0]
                cur = cur if cur in ids else None
                if cur:
                    seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip().upper())
        return {k: "".join(v) for k, v in seqs.items()}


    scaffolds = sm.input.get("scaffolds")
    if not status["circular"] and scaffolds and Path(scaffolds).stat().st_size:
        query, min_pid = (sm.input.fasta, 98) if go else (sm.input.seed, 95)
        if Path(query).stat().st_size:
            qlen = go_len or sum(len(s) for _, s in read_fasta(query))
            hits = blast_hits(query, scaffolds, min_pid)
            need = 0.9 * qlen if go else 0.5 * qlen
            cand = {sid for sid, bp in hits.items() if bp >= need}
            seqs = scaffold_seqs(scaffolds, cand)
            # circularize (no reads): assembly-graph closure, end overlap, or several pieces joined end to end
            graph_dir = next((d for d in (Path(scaffolds).parent / "spades", Path(scaffolds).parent)
                              if (d / "scaffolds.paths").exists()), None)
            rows, outs = circularize(list(seqs.items()), graph_dir)
            how = {r["contig"]: r for r in rows}
            if any(r["method"] == "join" for r in rows):
                seqs = {"joined": outs[0][1]}
            best = None
            for sid, s in seqs.items():
                r = how.get(sid) or next(r for r in rows if r["method"] == "join")
                circle = dict(outs).get(sid, s) if r["method"] != "not_circular" else None
                if circle and len(circle) >= 0.98 * go_len:
                    score = (1, len(circle))
                elif len(s) >= go_len + 100:
                    score = (0, len(s))           # used only when GetOrganelle gave nothing
                else:
                    continue
                if best is None or score > best[0]:
                    best = (score, sid, circle or s, bool(circle))
            if best and not best[3] and go:
                note = f"SPAdes {best[1]}: linear {len(best[2])} bp, not used (longer ends not reliable)"
            elif best:
                _, sid, seq, circ = best
                final = [(f"{sample}|mt|1 {'circular' if circ else 'linear'} length={len(seq)} source=spades:{sid}", seq)]
                source = "spades"
                r = how.get(sid) or next(r for r in rows if r["method"] == "join")
                note = (f"SPAdes {sid} {'circular (' + r['method'] + ')' if circ else 'linear'} {len(seq)} bp vs "
                        f"GetOrganelle {go_len} bp" + (f"; repeat units {r['repeats']}" if r["repeats"] else ""))
            print(f"candidates {sorted(cand)}; chosen: {note or 'GetOrganelle'}", file=sys.stderr)

    with open(sm.output.fasta, "w") as fh:
        write_fasta(final, fh)
    status.update(source=source, getorganelle_length=go_len, length=sum(len(s) for _, s in final),
                  circular=status["circular"] or (source == "spades" and "circular" in final[0][0]), note=note)
    json.dump(status, open(sm.output.status, "w"), indent=1)


if __name__ == "__main__" and "snakemake" in globals():
    main(snakemake)  # noqa: F821

#!/bin/bash
# Every command and every step of wgs_skim, checked quickly (minutes): sub-commands on real small inputs, each
# --steps value and SPAdes mode as a dry run (the jobs Snakemake would run), and the batch Snakefiles.
#   tests/cli_matrix.sh READS_DIR RESULTS_DIR [OUT_TSV]
# READS_DIR: a small reads directory (tests/e2e.sh uses the same); RESULTS_DIR: a finished `wgs_skim run` of it
# (for nuclear / compare / circularize / status). Writes one row per check: command, expected, result, seconds.
set -uo pipefail
reads=${1:?usage: tests/cli_matrix.sh READS_DIR RESULTS_DIR [OUT_TSV]}
res=${2:?usage: tests/cli_matrix.sh READS_DIR RESULTS_DIR [OUT_TSV]}
tsv=${3:-cli_matrix.tsv}
tmp=$(mktemp -d)
echo -e "check\tcommand\texpected\tresult\tseconds" > "$tsv"
fail=0

check() {      # check NAME EXPECT(grep pattern in the output) COMMAND...
    local name=$1 expect=$2; shift 2
    local t0=$(date +%s) out
    out=$("$@" 2>&1)
    local ok=FAIL
    echo "$out" | grep -Eq -- "$expect" && ok=ok
    [ $ok = FAIL ] && fail=$((fail + 1))
    printf "%s\t%s\t%s\t%s\t%s\n" "$name" "$*" "$expect" "$ok" "$(( $(date +%s) - t0 ))" >> "$tsv"
    printf "  %-4s %s\n" "$ok" "$name"
}
dry() {        # dry NAME EXPECTED_RULES... -- ARGS : the rules Snakemake would run for these arguments
    local name=$1; shift
    local rules=()
    while [ "$1" != "--" ]; do rules+=("$1"); shift; done; shift
    local t0=$(date +%s) out ok=ok
    out=$(wgs_skim run -r "$reads" -o "$tmp/dry" -n "$@" 2>&1)
    for r in "${rules[@]}"; do echo "$out" | grep -q "^rule $r:" || ok=FAIL; done
    echo "$out" | grep -q "^total" || ok=FAIL
    [ $ok = FAIL ] && fail=$((fail + 1))
    printf "%s\twgs_skim run -n %s\t%s\t%s\t%s\n" "$name" "$*" "rules: ${rules[*]}" "$ok" "$(( $(date +%s) - t0 ))" >> "$tsv"
    printf "  %-4s %s\n" "$ok" "$name"
}

echo "[matrix] sub-commands"
check "version"                 "wgs_skim [0-9]"            wgs_skim --version
check "samples (reads dir)"     "[0-9]+ samples"           wgs_skim samples -r "$reads"
check "samples (no argument)"   "give the reads directory" wgs_skim samples
check "config"                  "reads_dir:"               wgs_skim config
check "run without reads"       "give the reads directory" wgs_skim run
check "every subcommand -h"     "usage: wgs_skim nuclear"  bash -c 'for c in run install samples circularize compare nuclear build-blastdb config; do wgs_skim $c -h || exit 1; done'

echo "[matrix] steps (dry runs: which jobs run)"
dry "everything"                fastp meangs getorganelle_mt mt_iterate mt_select spades busco jellyfish genomescope mitom_annotate getorganelle_nr nr_annotate getorganelle_histone histone_annotate nuc_markers nuc_copies nuc_compare blast summary --
dry "fastp"                     fastp -- --steps fastp
dry "meangs"                    fastp meangs -- --steps meangs
dry "mito (MEANGS-GetOrganelle)" meangs getorganelle_mt mt_iterate mt_select -- --steps mito
dry "getorganelle (no MEANGS)"  getorganelle_mt mt_iterate mt_select -- --steps getorganelle
dry "annotate"                  mitom_annotate -- --steps annotate
dry "topology"                  topology -- --steps topology
dry "nr (rDNA + histones)"      getorganelle_nr nr_annotate getorganelle_histone histone_annotate -- --steps nr
dry "spades"                    spades -- --steps spades
dry "spades-busco"              spades busco -- --steps spades,busco
dry "spades modes"              spades_mode busco_assembly -- --steps spades,busco --spades-modes genome,meta,plasmid,isolate
dry "rna (rnaSPAdes + Trinity + BUSCO)" spades_mode trinity busco_assembly -- --steps rna
dry "trinity"                   trinity -- --steps trinity
dry "genome_size"               jellyfish genomescope -- --steps genome_size
dry "nuclear"                   nuc_markers nuc_copies nuc_compare -- --steps nuclear
dry "blast"                     blast -- --steps blast
dry "summary alone"             summary -- --steps summary
dry "memory option"             fastp -- --steps fastp --mem-gb 50
check "unknown step"            "unknown step"             wgs_skim run -r "$reads" -o "$tmp/dry" -n --steps foo
check "unknown SPAdes mode"     "unknown SPAdes mode"      wgs_skim run -r "$reads" -o "$tmp/dry" -n --steps spades --spades-modes foo
check "--mem-gb reaches snakemake" "--resources mem_mb=51200" wgs_skim run -r "$reads" -o "$tmp/dry" -n --steps fastp --mem-gb 50

echo "[matrix] commands on the finished results ($res)"
s=$(ls "$res/qc" | head -1)
round=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['round'])" "$res/mito/$s/$s.getorganelle.json")
check "circularize (no reads)"  "contig|circular|method"    bash -c "wgs_skim circularize -i $res/mito/$s/$s.mt.fasta -o $tmp/circ --prefix $s && cat $tmp/circ/*.tsv"
check "compare (one sample)"    "getorganelle .* bp vs"     wgs_skim compare -a "$res/mito/$s/$round" -b "final=$res/mito/$s/$s.mt.fasta" -o "$tmp/cmp" --prefix "$s"
check "nuclear (Snakemake, reuse)" "Nothing to be done|steps \\(100%\\) done" wgs_skim nuclear --results "$res" -o "$res/nuclear_cli" -j 8
check "nuclear --status"        "rdna (compared|annotated)" wgs_skim nuclear --status -o "$res/nuclear_cli"
check "nuclear dry run"         "Nothing to be done|^total"  wgs_skim nuclear --results "$res" -o "$res/nuclear_cli" -n
check "unit tests"              "passed"                    python -m pytest -q "$(dirname "$(readlink -f "$0")")"

rm -rf "$tmp"
echo "[matrix] $fail failed -> $tsv"
exit $fail

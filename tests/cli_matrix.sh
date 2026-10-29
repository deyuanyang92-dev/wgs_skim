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

echo "[matrix] single tools and combinations (dry runs)"
asm=$tmp/asm; mkdir -p $asm $tmp/mt
for x in $(ls "$res/qc"); do
    ln -sf "$(realpath "$res/assembly/$x/$x.scaffolds.fasta")" "$asm/$x.fasta"
    ln -sf "$(realpath "$res/mito/$x/$x.mt.fasta")" "$tmp/mt/$x.fasta"
done
s0=$(ls "$res/qc" | head -1)
ref=$(realpath "$res/annotation/$s0/$s0.mt.gb")
tool() {       # tool NAME EXPECTED_RULES... -- wgs_skim ARGS : the rules a tool command would run
    local name=$1; shift
    local rules=()
    while [ "$1" != "--" ]; do rules+=("$1"); shift; done; shift
    local t0=$(date +%s) out ok=ok
    out=$(wgs_skim "$@" -n 2>&1)
    for r in "${rules[@]}"; do echo "$out" | grep -q "^rule $r:" || ok=FAIL; done
    echo "$out" | grep -q "^total" || ok=FAIL
    [ $ok = FAIL ] && fail=$((fail + 1))
    printf "%s\twgs_skim %s -n\t%s\t%s\t%s\n" "$name" "$*" "rules: ${rules[*]}" "$ok" "$(( $(date +%s) - t0 ))" >> "$tsv"
    printf "  %-4s %s\n" "$ok" "$name"
}
tool "tool fastp"                fastp -- fastp -r "$reads" -o $tmp/t1
tool "tool meangs"               meangs -- meangs -r "$reads" -o $tmp/t2
tool "tool getorganelle"         getorganelle_mt mt_iterate mt_select -- getorganelle -r "$reads" -o $tmp/t3
tool "tool getorganelle --meangs" meangs getorganelle_mt -- getorganelle -r "$reads" -o $tmp/t3b --meangs
tool "tool getorganelle --target nr" getorganelle_nr getorganelle_histone -- getorganelle -r "$reads" -o $tmp/t4 --target nr
tool "tool spades --mode genome,rna" spades spades_mode -- spades -r "$reads" -o $tmp/t5 --mode genome,rna
tool "tool trinity"              trinity -- trinity -r "$reads" -o $tmp/t6
tool "tool genome-size"          jellyfish genomescope -- genome-size -r "$reads" -o $tmp/t7
tool "tool busco -i"             assembly_given busco -- busco -i $asm -o $tmp/t8
tool "tool busco -m transcriptome" assembly_given busco -- busco -i $asm -o $tmp/t8b -m transcriptome
tool "tool mitom -i"             mt_given mitom_annotate -- mitom -i $tmp/mt -o $tmp/t9
tool "tool blast -i"             mt_given blast -- blast -i $tmp/mt -o $tmp/t10
tool "tool mitofinder -i"        assembly_given mitofinder mitofinder_summary -- mitofinder -i $asm -o $tmp/t11 --ref "$ref"
tool "tool mitofinder -r"        mitofinder -- mitofinder -r "$reads" -o $tmp/t12 --ref "$ref"
check "tool mitofinder without --ref" "needs a reference" wgs_skim mitofinder -i $asm -o $tmp/t13
check "tool busco without -i"   "give the FASTA"           wgs_skim busco -o $tmp/t14
dry "combo meangs-getorganelle" meangs getorganelle_mt mt_iterate mt_select -- --steps meangs-getorganelle
dry "combo spades-busco"        spades busco -- --steps spades-busco
dry "combo spades-mitofinder (own ref)" spades mitofinder mitom_annotate -- --steps meangs-getorganelle,spades-mitofinder
dry "combo rnaspades-busco"     spades_mode busco_assembly -- --steps rnaspades-busco
dry "combo trinity-busco"       trinity busco_assembly -- --steps trinity-busco
dry "combo rnaspades-mitofinder" spades_mode mitofinder_input mitofinder mitofinder_summary mitofinder_blast -- --steps rnaspades-mitofinder --mitofinder-ref "$ref"
dry "combo trinity-mitofinder"  trinity mitofinder_input mitofinder -- --steps trinity-mitofinder --mitofinder-ref "$ref"
tool "tool mitofinder --assembler rnaspades,trinity" spades_mode trinity mitofinder -- mitofinder -r "$reads" -o $tmp/t15 --ref "$ref" --assembler rnaspades,trinity
tool "tool mitofinder --no-check" mitofinder_summary -- mitofinder -i $asm -o $tmp/t17 --ref "$ref" --no-check
check "tool mitofinder bad --assembler" "unknown --assembler" wgs_skim mitofinder -r "$reads" -o $tmp/t16 --ref "$ref" --assembler foo
check "combo spades-mitofinder without ref" "needs a reference" wgs_skim run -r "$reads" -o "$tmp/dry" -n --steps spades-mitofinder

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

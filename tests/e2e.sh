#!/bin/bash
# End-to-end test: the whole workflow on a small dataset, then a check that every sample has every output.
#   tests/e2e.sh READS_DIR [OUTDIR]
# READS_DIR: two or more small samples (1 M read pairs each is enough: ~15 min with 32+ cores), e.g. made with
#   for s in A B; do for r in 1 2; do zcat $s.R$r.fq.gz | head -n 4000000 | gzip > small/${s}_R$r.fq.gz; done; done
# Test settings keep it short: one reseeded mitogenome round, 200 k read pairs for the topology check.
# Exit 0 = every check passed; the failed checks are listed otherwise.
set -uo pipefail
reads=${1:?usage: tests/e2e.sh READS_DIR [OUTDIR]}
out=${2:-$(dirname "$(readlink -f "$reads")")/e2e_$(date +%Y%m%d_%H%M)}
here=$(dirname "$(readlink -f "$0")")
start=$(date +%s)
wgs_skim run -r "$reads" -o "$out" --config mt_iter_rounds=1 topology_pairs=200000 > "$out.log" 2>&1
rc=$?
mins=$(( ($(date +%s) - start) / 60 ))
fail=0
check() {      # check DESCRIPTION FILE: the file exists and is not empty
    if [ -s "$2" ]; then echo "  ok    $1"; else echo "  FAIL  $1 ($2)"; fail=$((fail + 1)); fi
}
echo "[e2e] wgs_skim run exit $rc after $mins min; results $out; log $out.log"
for s in $(wgs_skim samples -r "$reads" 2>/dev/null | cut -f1); do
    echo "[e2e] $s"
    check "fastp"                "$out/qc/$s/$s.fastp.json"
    check "mitogenome"           "$out/mito/$s/$s.mt.fasta"
    check "MitoM annotation"     "$out/annotation/$s/$s.mt.gb"
    check "rDNA (GetOrganelle)"  "$out/rdna/$s/$s.nr.fasta"
    check "rDNA GenBank"         "$out/rdna/$s/$s.rdna.gbf"
    check "histone GenBank"      "$out/histone/$s/$s.histone.gbf"
    check "SPAdes"               "$out/assembly/$s/$s.scaffolds.fasta"
    check "BUSCO"                "$out/assembly/$s/$s.busco.txt"
    check "genome size"          "$out/genome_size/$s/$s.genomescope.txt"
    check "nuclear copies"       "$out/nuclear/$s/$s.copy_number.tsv"
    check "GetOrganelle vs SPAdes" "$out/nuclear/$s/$s.nuc_compare.tsv"
done
check "BLAST"   "$out/blast/blast_topn.tsv"
check "summary" "$out/summary/samples.tsv"
check "report"  "$out/summary/report.html"
n=$(($(wc -l < "$out/summary/samples.tsv" 2>/dev/null || echo 1) - 1))
echo "[e2e] summary rows: $n"
# the same results through `wgs_skim nuclear` (everything reused: seconds)
wgs_skim nuclear --results "$out" -o "$out/nuclear_cli" > "$out/nuclear_cli.log" 2>&1
check "wgs_skim nuclear summary" "$out/nuclear_cli/summary.tsv"
python -m pytest -q "$here" > "$out/unit_tests.log" 2>&1 && echo "  ok    unit tests" || { echo "  FAIL  unit tests"; fail=$((fail + 1)); }
[ $rc -eq 0 ] || fail=$((fail + 1))
echo "[e2e] $fail failed check(s) in $mins min"
exit $fail

#!/bin/bash
# Every tool on its own (wgs_skim <tool>) and the recommended combinations, run for real on small data; each
# check: the command succeeds and its per-sample outputs exist and are not empty.
#   tests/tools_e2e.sh RAW_READS_DIR RESULTS_DIR OUT_DIR [REF_GB] [RNA_READS_DIR]
# RAW_READS_DIR: 2+ small samples (as tests/e2e.sh); RESULTS_DIR: a finished `wgs_skim run` of them (its clean
# reads, assemblies and mitogenomes are the inputs of the single tools); REF_GB: MitoFinder reference
# (default: the first sample's MitoM annotation in RESULTS_DIR); RNA_READS_DIR: small RNA-seq samples for the
# transcriptome checks (rnaSPAdes / Trinity, BUSCO, MitoFinder on the transcripts). Writes OUT_DIR/tools_e2e.tsv;
# exit code = failures.
set -uo pipefail
raw=$(realpath "${1:?usage: tests/tools_e2e.sh RAW_READS_DIR RESULTS_DIR OUT_DIR [REF_GB]}")
res=$(realpath "${2:?}")
out=$(realpath -m "${3:?}")
mem=${MEM_GB:-120}
jobs=${JOBS:-32}
mkdir -p "$out/in/asm" "$out/in/mt"
samples=$(ls "$res/qc")
first=$(echo "$samples" | head -1)
ref=$(realpath "${4:-$res/annotation/$first/$first.mt.gb}")
rna=${5:+$(realpath "$5")}
for s in $samples; do                     # inputs of the FASTA tools: files named by sample, as a user has them
    ln -sf "$res/assembly/$s/$s.scaffolds.fasta" "$out/in/asm/$s.fasta"
    ln -sf "$res/mito/$s/$s.mt.fasta" "$out/in/mt/$s.fasta"
done
clean="$res/qc"                           # clean reads: <s>/<s>.R1.fq.gz
tsv="$out/tools_e2e.tsv"
echo -e "check\tcommand\tresult\tseconds\tmissing" > "$tsv"
fail=0

check() {   # check NAME "OUTPUT PATTERNS (with {s})" COMMAND...
    local name=$1 pats=$2; shift 2
    local t0=$(date +%s) log="$out/logs/${name// /_}.log" ok=ok missing=""
    mkdir -p "$out/logs"
    "$@" -j "$jobs" --mem-gb "$mem" > "$log" 2>&1 || { ok=FAIL; missing="exit $?"; }
    for s in $samples; do for p in $pats; do
        f=${p//\{s\}/$s}
        [ -s "$f" ] || { ok=FAIL; missing+=" ${f#$out/}"; }
    done; done
    [ $ok = FAIL ] && fail=$((fail + 1))
    printf "%s\t%s\t%s\t%s\t%s\n" "$name" "$*" "$ok" "$(( $(date +%s) - t0 ))" "$missing" >> "$tsv"
    printf "  %-4s %-34s %5ss %s\n" "$ok" "$name" "$(( $(date +%s) - t0 ))" "$missing"
}

o=$out/single
echo "[tools] single tools ($samples)"
check "fastp" "$o/fastp/qc/{s}/{s}.R1.fq.gz $o/fastp/qc/{s}/{s}.fastp.html" \
    wgs_skim fastp -r "$raw" -o "$o/fastp"
check "meangs (clean reads)" "$o/meangs/mito/{s}/meangs/{s}.seed.fasta" \
    wgs_skim meangs -r "$clean" -o "$o/meangs"
# GetOrganelle's own animal_mt seeds: a sample far from them can give no sequence (data2 pisione: 0 sequences,
# 15,776 bp with the MEANGS seed) - only the run record is required here
check "getorganelle mt (own seeds, no rounds)" "$o/go/mito/{s}/{s}.mt.json" \
    wgs_skim getorganelle -r "$clean" -o "$o/go" --rounds 0
check "getorganelle mt --meangs" "$o/gom/mito/{s}/{s}.mt.fasta" \
    wgs_skim getorganelle -r "$clean" -o "$o/gom" --meangs
check "getorganelle nr (rDNA+histones)" "$o/gonr/rdna/{s}/{s}.rdna.gbf" \
    wgs_skim getorganelle -r "$clean" -o "$o/gonr" --target nr
check "spades genome" "$o/spades/assembly/{s}/{s}.scaffolds.fasta" \
    wgs_skim spades -r "$clean" -o "$o/spades"
check "spades meta" "$o/spadesm/assembly/{s}/{s}.meta.fasta" \
    wgs_skim spades -r "$clean" -o "$o/spadesm" --mode meta
check "busco genome (-i)" "$o/busco/assembly/{s}/{s}.busco.txt" \
    wgs_skim busco -i "$out/in/asm" -o "$o/busco"
check "mitom (-i)" "$o/mitom/annotation/{s}/{s}.mt.gb" \
    wgs_skim mitom -i "$out/in/mt" -o "$o/mitom"
check "blast full (-i)" "$o/blast/blast/blast_topn.tsv" \
    wgs_skim blast -i "$out/in/mt" -o "$o/blast"
check "blast cox1 (-i, MitoM first)" "$o/blastc/blast/blast_topn.tsv $o/blastc/annotation/{s}/{s}.mt.gb" \
    wgs_skim blast -i "$out/in/mt" -o "$o/blastc" --genes cox1
check "mitofinder (-i assemblies)" "$o/mf/mitofinder/given/{s}/{s}.mitofinder.fasta $o/mf/mitofinder/given/mitofinder_summary.tsv" \
    wgs_skim mitofinder -i "$out/in/asm" -o "$o/mf" --ref "$ref"
check "mitofinder (-r reads, MEGAHIT)" "$o/mfr/mitofinder/megahit/{s}/{s}.mitofinder.fasta" \
    wgs_skim mitofinder -r "$clean" -o "$o/mfr" --ref "$ref"
check "mitofinder bioconda 1.4.1" "$o/mfc/mitofinder/given/{s}/{s}.mitofinder.fasta" \
    wgs_skim mitofinder -i "$out/in/asm" -o "$o/mfc" --ref "$ref" --config mitofinder_sif=
check "genome-size" "$o/gs/genome_size/{s}/{s}.genomescope.txt" \
    wgs_skim genome-size -r "$clean" -o "$o/gs"

o=$out/combo
echo "[tools] combinations (raw reads)"
check "spades-busco" "$o/sb/assembly/{s}/{s}.busco.txt" \
    wgs_skim run -r "$raw" -o "$o/sb" --steps spades-busco
check "spades-mitofinder (own MitoM ref)" "$o/smf/mitofinder/spades/{s}/{s}.mitofinder.fasta $o/smf/mito/{s}/{s}.mt.fasta" \
    wgs_skim run -r "$raw" -o "$o/smf" --steps meangs-getorganelle,spades-mitofinder

if [ -n "$rna" ]; then
    samples=$(wgs_skim samples -r "$rna" 2>/dev/null | cut -f1)
    o=$out/rna
    echo "[tools] transcriptome ($samples)"
    check "rnaspades-busco" "$o/rb/assembly/{s}/{s}.rna.fasta $o/rb/assembly/{s}/{s}.rna.busco.txt" \
        wgs_skim run -r "$rna" -o "$o/rb" --steps rnaspades-busco
    check "trinity-busco" "$o/tb/assembly/{s}/{s}.trinity.fasta $o/tb/assembly/{s}/{s}.trinity.busco.txt" \
        wgs_skim run -r "$rna" -o "$o/tb" --steps trinity-busco
    check "rnaspades-mitofinder" "$o/rm/mitofinder/rna/{s}/{s}.mitofinder.fasta $o/rm/mitofinder/rna/mitofinder_summary.tsv" \
        wgs_skim run -r "$rna" -o "$o/rm" --steps rnaspades-mitofinder --mitofinder-ref "$ref"
    check "mitofinder --assembler trinity (tool)" "$o/tm/mitofinder/trinity/{s}/{s}.mitofinder.fasta" \
        wgs_skim mitofinder -r "$o/rb/qc" -o "$o/tm" --ref "$ref" --assembler trinity
    check "busco -m transcriptome (-i)" "$o/bt/assembly/{s}.rna/{s}.rna.busco.txt" \
        wgs_skim busco -i "$o/rb/assembly"/*/*.rna.fasta -o "$o/bt" -m transcriptome
fi

echo "[tools] $fail failed -> $tsv"
exit $fail

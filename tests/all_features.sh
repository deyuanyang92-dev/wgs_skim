#!/bin/bash
# Every wgs_skim function on the simulated samples (tests/sim/simulate.py), one output directory per function:
#   tests/all_features.sh SIM_DIR OUT_DIR [dna|rna|all]
# SIM_DIR: simulate.py output (dna/, rna/, truth/, truth.tsv). OUT_DIR/NN_<function>/ holds that function's results
# exactly as wgs_skim writes them (<tool>/<sample>/ + <tool>/<tool>_summary.tsv, pipeline_info/);
# OUT_DIR/features.<lane>.tsv lists the runs. Then tests/sim/check_truth.py and tests/sim/report.py (README.md).
# The DNA and RNA functions are two lanes (run them at the same time: one call with dna, one with rna).
# ONLY="13 36" reruns just these functions (their directories are replaced; the other directories are inputs).
# Then: python tests/sim/check_truth.py SIM_DIR OUT_DIR  (assemblies against the truth -> OUT_DIR/truth_check.tsv)
set -uo pipefail
sim=$(realpath "${1:?usage: tests/all_features.sh SIM_DIR OUT_DIR [dna|rna|all]}")
out=$(realpath -m "${2:?}")
lane=${3:-all}
jobs=${JOBS:-64}
mem=${MEM_GB:-300}
dna=$sim/dna rna=$sim/rna
mkdir -p "$out/inputs"
samples=$(wgs_skim samples -r "$dna" 2>/dev/null | cut -f1)
first=$(echo "$samples" | head -1)
py=$(head -1 "$(command -v wgs_skim)" | sed 's/^#!//')      # the Python wgs_skim is installed in (API check)
tsv=$out/features.$lane${ONLY:+.rerun_$(date +%Y%m%d_%H%M%S)}.tsv   # each rerun its own record
echo -e "id\tfunction\tresult\tseconds\toutput\tcommand\tmissing\twhat" > "$tsv"
fail=0

# feat ID NAME "OUTPUT PATTERNS ({s} = each sample, relative to the output dir)" "WHAT TO CHECK" COMMAND...
# (@ in the command = this function's output directory)
feat() {
    local id=$1 name=$2 pats=$3 what=$4; shift 4
    local d="$out/${id}_${name}" t0=$(date +%s) ok=ok missing="" cmd=()
    if [ -n "${ONLY:-}" ]; then
        [[ " $ONLY " == *" $id "* ]] || return 0
        rm -rf "$d"
    fi
    for x in "$@"; do cmd+=("${x//@/$d}"); done
    mkdir -p "$d/pipeline_info"
    echo "[$id] $name ..." >&2
    "${cmd[@]}" -j "$jobs" --mem-gb "$mem" > "$d/pipeline_info/wgs_skim.log" 2>&1 || { ok=FAIL; missing="exit $?"; }
    for s in $samples; do for p in $pats; do
        f="$d/${p//\{s\}/$s}"
        compgen -G "$f" > /dev/null && [ -s "$(compgen -G "$f" | head -1)" ] || { ok=FAIL; missing+=" ${p//\{s\}/$s}"; }
    done; done
    local sec=$(( $(date +%s) - t0 ))
    [ $ok = FAIL ] && fail=$((fail + 1))
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$id" "$name" "$ok" "$sec" "$d" "${cmd[*]}" "$missing" "$what" >> "$tsv"
    printf "  %-4s %s_%s %ss %s\n" "$ok" "$id" "$name" "$sec" "$missing" >&2
}
# extra arguments some commands do not take (-j / --mem-gb): wrap them
noj() { local a=(); while [ $# -gt 0 ]; do case $1 in -j|--mem-gb) shift 2;; *) a+=("$1"); shift;; esac; done; "${a[@]}"; }

links() {   # links DIR PATTERN: OUT/inputs/DIR/<s>.fasta -> the per-sample file (input of the -i functions)
    mkdir -p "$out/inputs/$1"
    for s in $samples; do f=${2//\{s\}/$s}; [ -s "$f" ] && ln -sf "$f" "$out/inputs/$1/$s.fasta"; done
}

if [ "$lane" = dna ] || [ "$lane" = all ]; then
echo "[features] DNA: $(echo $samples)" >&2
feat 01 samples "samples.tsv" "sample names detected from a reads directory (5 samples, R1/R2 paired): samples.tsv" \
    bash -c "wgs_skim samples -r $dna > @/samples.tsv; cat @/samples.tsv" --
feat 02 fastp "fastp/{s}/{s}.R1.fq.gz fastp/{s}/{s}.fastp.html fastp/fastp_summary.tsv" \
    "fastp, once: clean reads and reports; fastp/fastp_summary.tsv = data and Q30 per sample. Every later DNA function reads these clean reads where they are" \
    wgs_skim fastp -r "$dna" -o @
clean=$out/02_fastp/fastp            # clean reads of every later DNA function (used where they are, no QC again)
feat 03 meangs "meangs/{s}/{s}.seed.fasta meangs/meangs_summary.tsv" "MEANGS mitochondrial seed per sample; meangs/meangs_summary.tsv" \
    wgs_skim meangs -r "$clean" -o @
feat 04 getorganelle_mt "getorganelle/{s}/{s}.getorganelle.json getorganelle/getorganelle_summary.tsv" \
    "GetOrganelle with its own animal_mt seeds; getorganelle/getorganelle_summary.tsv = length, circular, rounds. Seeds far from the sample can give nothing - the reason for 05" \
    wgs_skim getorganelle -r "$clean" -o @
feat 05 meangs-getorganelle "mitogenome/{s}/{s}.mt.fasta mitogenome/mitogenome_summary.tsv" \
    "recommended: MEANGS seed + GetOrganelle, reseeded while not circular; mitogenome/mitogenome_summary.tsv = final length and topology (truth: S1-S4 circular, S5 linear)" \
    wgs_skim run -r "$clean" --no-qc -o @ --steps meangs-getorganelle
feat 06 getorganelle_nr "rdna/{s}/{s}.rdna.gbf histone/{s}/{s}.histone.gbf rdna/rdna_summary.tsv" \
    "rDNA unit (18S-ITS1-5.8S-ITS2-28S) and histone cluster per sample, GenBank annotation (.gbf)" \
    wgs_skim getorganelle -r "$clean" -o @ --target nr
feat 07 spades_genome "spades/{s}/{s}.scaffolds.fasta spades/spades_summary.tsv" "SPAdes genome assembly; spades/spades_summary.tsv = total length (truth ~20 Mb), N50" \
    wgs_skim spades -r "$clean" -o @
feat 08 spades_meta "spades/{s}/{s}.meta.fasta spades/spades_summary.tsv" "metaSPAdes; compare spades_summary.tsv with 07" \
    wgs_skim spades -r "$clean" -o @ --mode meta
feat 09 genome-size "genomescope/{s}/{s}.genomescope.txt genomescope/genomescope_summary.tsv" \
    "GenomeScope; genomescope/genomescope_summary.tsv = genome size (truth ~20 Mb), heterozygosity (truth 0.2 %)" \
    wgs_skim genome-size -r "$clean" -o @
links asm "$out/07_spades_genome/spades/{s}/{s}.scaffolds.fasta"
links mt "$out/05_meangs-getorganelle/mitogenome/{s}/{s}.mt.fasta"
feat 10 busco_genome "busco/{s}/{s}.busco.txt busco/busco_summary.tsv" "BUSCO of given assemblies (-i); busco/busco_summary.tsv (a 20 Mb piece holds few genes: low completeness is expected)" \
    wgs_skim busco -i "$out/inputs/asm" -o @
feat 11 mitom "mitom/{s}/{s}.mt.gb mitom/mitom_summary.tsv" "MitoM annotation of given mitogenomes (-i); mitom/mitom_summary.tsv = 13 PCG + 2 rRNA + 22 tRNA" \
    wgs_skim mitom -i "$out/inputs/mt" -o @
ref=$out/11_mitom/mitom/$first/$first.mt.gb
feat 12 blast "blast/blast_topn.tsv" "BLAST of cox1 and 16S against the local metazoa_marker database: best hit species and identity" \
    wgs_skim blast -i "$out/inputs/mt" -o @ --genes cox1,16s
feat 13 mitofinder_assembly "mitofinder_assembly/{s}/{s}.mitofinder.fasta mitofinder_assembly/mitofinder_summary.tsv" \
    "MitoFinder on assemblies, one distant reference (S1). mitofinder_assembly/: <sample>/ + mitofinder_summary.tsv (gene table) + genes/ + contamination/" \
    wgs_skim mitofinder -i "$out/inputs/asm" -o @ --ref "$ref"
feat 14 mitofinder_reads "mitofinder_megahit/{s}/{s}.mitofinder.fasta mitofinder_megahit/mitofinder_summary.tsv" \
    "MitoFinder assembling the reads itself (MEGAHIT); compare with 13" \
    wgs_skim mitofinder -r "$clean" -o @ --ref "$ref"
feat 15 mitofinder_bioconda "mitofinder_assembly/{s}/{s}.mitofinder.fasta mitofinder_assembly/mitofinder_summary.tsv" \
    "bioconda MitoFinder 1.4.1 (no Singularity image): tRNAs from tRNAscan-SE; compare with 13" \
    wgs_skim mitofinder -i "$out/inputs/asm" -o @ --ref "$ref" --config mitofinder_sif=
if [ -d "$sim/refs" ]; then          # related species' mitogenomes (not the simulated ones): SIM_DIR/refs/*.gb
feat 16 mitofinder_multiref "mitofinder_assembly/{s}/{s}.mitofinder.fasta mitofinder_assembly/mitofinder_summary.tsv" \
    "as 13 with 5 related references (--ref directory): gene counts of the distant samples against 13 (recommended use)" \
    wgs_skim mitofinder -i "$out/inputs/asm" -o @ --ref "$sim/refs"
fi
feat 17 spades-busco "spades/{s}/{s}.scaffolds.fasta busco/{s}/{s}.busco.txt busco/busco_summary.tsv" "combination SPAdes -> BUSCO" \
    wgs_skim run -r "$clean" --no-qc -o @ --steps spades-busco
feat 18 spades-mitofinder "mitofinder_spades/{s}/{s}.mitofinder.fasta mitogenome/{s}/{s}.mt.fasta mitofinder_spades/mitofinder_summary.tsv" \
    "combination: GetOrganelle mitogenome + SPAdes-MitoFinder (reference = each sample's own MitoM annotation)" \
    wgs_skim run -r "$clean" --no-qc -o @ --steps meangs-getorganelle,spades-mitofinder
feat 19 circularize "{s}/{s}.circularize.tsv" "circularization without reads: method (walk / overlap / join) and circle length per sample" \
    noj bash -c "for s in $(echo $samples); do wgs_skim circularize -i $out/18_spades-mitofinder/mitofinder_spades/\$s/\$s.mitofinder.fasta --spades $out/18_spades-mitofinder/spades/\$s/spades -o @/\$s --prefix \$s || exit 1; done" --
feat 20 compare "{s}/{s}.compare.tsv all.compare.html" "GetOrganelle vs SPAdes-MitoFinder: verdict in all.compare.tsv; differences and figure per sample" \
    wgs_skim compare --results "$out/18_spades-mitofinder" -o @
feat 21 run_all "summary/samples.tsv summary/report.html fastp/fastp_summary.tsv" \
    "the full workflow from raw reads (every step): summary/samples.tsv, plus each tool's *_summary.tsv" \
    wgs_skim run -r "$dna" -o @
feat 22 nuclear "summary.tsv" "nuclear markers: GetOrganelle vs SPAdes-MitoM per gene (verdict)" \
    wgs_skim nuclear --results "$out/21_run_all" -o @
printf "sample\tfastq_1\tfastq_2\tscaffolds\tmitogenome\n" > "$out/inputs/sheet_given.tsv"
for s in $samples; do printf "%s\t%s\t%s\t%s\t%s\n" "$s" "$clean/$s/$s.R1.fq.gz" "$clean/$s/$s.R2.fq.gz" \
    "$out/inputs/asm/$s.fasta" "$out/inputs/mt/$s.fasta" >> "$out/inputs/sheet_given.tsv"; done
feat 23 sheet_given "summary/samples.tsv mitom/{s}/{s}.mt.gb busco/busco_summary.tsv" \
    "sample sheet with existing assemblies and mitogenomes: no assembly, BUSCO / MitoM / BLAST on them" \
    wgs_skim run -s "$out/inputs/sheet_given.tsv" --no-qc -o @ --steps busco,annotate,blast,summary
cat "$out"/inputs/mt/*.fasta > "$out/inputs/merged_mt.fasta"
feat 24 split_by_id "mitom/{s}*/*.mt.gb" "one merged FASTA of several samples, split by sequence name, annotated per sample" \
    wgs_skim mitom -i "$out/inputs/merged_mt.fasta" --split-by-id -o @
feat 25 python_api "busco/{s}/{s}.busco.txt busco/busco_summary.tsv" "Python interface wgs_skim.tool(...): same as 10" \
    noj "$py" -c "import sys, wgs_skim; sys.exit(wgs_skim.tool('busco', input=['$out/inputs/asm'], outdir='@', cores=$jobs, mem_gb=$mem))" --
if command -v nextflow > /dev/null; then
feat 26 nextflow "mitogenome/{s}/{s}.mt.fasta summary/samples.tsv" "the Nextflow engine (raw reads, default steps): compare summary/samples.tsv with 21" \
    wgs_skim run -r "$dna" -o @ --engine nextflow
fi
fi

if [ "$lane" = rna ] || [ "$lane" = all ]; then
echo "[features] RNA: $(echo $samples)" >&2
# MitoFinder reference of the RNA lane (runs beside the DNA lane, so not a DNA-lane output): the first sample's
# source GenBank (truth.tsv mt_source) - as a related species' mitogenome would be in a real study
ref=$(awk -F'\t' 'NR == 2 {for (i = 1; i <= NF; i++) if (h[i] == "mt_source") print $i} NR == 1 {for (i = 1; i <= NF; i++) h[i] = $i}' "$sim/truth.tsv")
feat 30 fastp_rna "fastp/{s}/{s}.R1.fq.gz fastp/fastp_summary.tsv" "fastp of the RNA-seq reads, once: every later RNA function reads these clean reads" \
    wgs_skim fastp -r "$rna" -o @
rclean=$out/30_fastp_rna/fastp
feat 31 rnaspades "spades/{s}/{s}.rna.fasta spades/spades_summary.tsv" "rnaSPAdes: transcripts per sample in spades/spades_summary.tsv" \
    wgs_skim spades -r "$rclean" -o @ --mode rna
feat 32 trinity "trinity/{s}/{s}.trinity.fasta trinity/trinity_summary.tsv" "Trinity: trinity/trinity_summary.tsv, compare with 31" \
    wgs_skim trinity -r "$rclean" -o @
links rna "$out/31_rnaspades/spades/{s}/{s}.rna.fasta"
feat 33 busco_transcriptome "busco/{s}/{s}.busco.txt busco/busco_summary.tsv" "BUSCO -m transcriptome (nuclear transcripts taken from a real RNA-seq assembly)" \
    wgs_skim busco -i "$out/inputs/rna" -o @ -m transcriptome
feat 34 rnaspades-busco "spades/{s}/{s}.rna.fasta busco/{s}/{s}.rna.busco.txt busco/busco_summary.tsv" "combination (recommended) rnaSPAdes -> BUSCO" \
    wgs_skim run -r "$rclean" --no-qc -o @ --steps rnaspades-busco
feat 35 trinity-busco "trinity/{s}/{s}.trinity.fasta busco/{s}/{s}.trinity.busco.txt busco/busco_summary.tsv" "combination Trinity -> BUSCO" \
    wgs_skim run -r "$rclean" --no-qc -o @ --steps trinity-busco
feat 36 rnaspades-mitofinder "mitofinder_rnaspades/{s}/{s}.mitofinder.fasta mitofinder_rnaspades/mitofinder_summary.tsv mitofinder_rnaspades/contamination/samples_check.tsv" \
    "mitochondrial genes from RNA-seq: mitofinder_summary.tsv (gene table, N = not found); contamination/samples_check.tsv: S5 must report COX1, CYTB, the others nothing" \
    wgs_skim run -r "$rclean" --no-qc -o @ --steps rnaspades-mitofinder --mitofinder-ref "$ref"
feat 37 trinity-mitofinder "mitofinder_trinity/{s}/{s}.mitofinder.fasta mitofinder_trinity/mitofinder_summary.tsv mitofinder_trinity/contamination/samples_check.tsv" \
    "as 36 with Trinity transcripts" \
    wgs_skim run -r "$rclean" --no-qc -o @ --steps trinity-mitofinder --mitofinder-ref "$ref"
feat 38 mitofinder_rna_tool "mitofinder_rnaspades/mitofinder_summary.tsv mitofinder_trinity/mitofinder_summary.tsv" \
    "the tool command wgs_skim mitofinder --assembler rnaspades,trinity (one directory per source)" \
    wgs_skim mitofinder -r "$rclean" -o @ --ref "$ref" --assembler rnaspades,trinity
fi

echo "[features] $lane: $fail failed -> $tsv" >&2
exit $fail

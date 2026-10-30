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
tsv=$out/features.$lane${ONLY:+.rerun}.tsv
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
feat 01 samples "samples.tsv" "样品名识别是否正确（5 个，R1/R2 配对）" \
    bash -c "wgs_skim samples -r $dna > @/samples.tsv; cat @/samples.tsv" --
feat 02 fastp "fastp/{s}/{s}.R1.fq.gz fastp/{s}/{s}.fastp.html fastp/fastp_summary.tsv" \
    "fastp/fastp_summary.tsv：每个样品 clean 数据量、Q30；fastp/<样品>/：clean reads 与报告（后面的功能都直接用这里的 clean reads）" \
    wgs_skim fastp -r "$dna" -o @
clean=$out/02_fastp/fastp            # clean reads of every later DNA function (used where they are, no QC again)
feat 03 meangs "meangs/{s}/{s}.seed.fasta meangs/meangs_summary.tsv" "meangs/meangs_summary.tsv：每个样品的线粒体种子长度" \
    wgs_skim meangs -r "$clean" -o @
feat 04 getorganelle_mt "getorganelle/{s}/{s}.getorganelle.json getorganelle/getorganelle_summary.tsv" \
    "GetOrganelle 自带种子（animal_mt）：getorganelle_summary.tsv 的长度、环化、轮次；种子太远时可能没有结果——这正是推荐 05 的原因" \
    wgs_skim getorganelle -r "$clean" -o @
feat 05 meangs-getorganelle "mitogenome/{s}/{s}.mt.fasta mitogenome/mitogenome_summary.tsv" \
    "推荐组合：MEANGS 种子 + GetOrganelle，未成环时迭代。mitogenome/mitogenome_summary.tsv：最终线粒体长度、是否环状（真值 S1-S4 环状，S5 线状）" \
    wgs_skim run -r "$clean" --no-qc -o @ --steps meangs-getorganelle
feat 06 getorganelle_nr "rdna/{s}/{s}.rdna.gbf histone/{s}/{s}.histone.gbf rdna/rdna_summary.tsv" \
    "rdna/、histone/：每个样品的 rDNA 单元（18S-ITS1-5.8S-ITS2-28S）与组蛋白簇，GenBank 注释（.gbf）" \
    wgs_skim getorganelle -r "$clean" -o @ --target nr
feat 07 spades_genome "spades/{s}/{s}.scaffolds.fasta spades/spades_summary.tsv" "spades/spades_summary.tsv：总长（真值 ~20 Mb）、N50" \
    wgs_skim spades -r "$clean" -o @
feat 08 spades_meta "spades/{s}/{s}.meta.fasta spades/spades_summary.tsv" "metaSPAdes 模式：spades_summary.tsv 与 07 比较" \
    wgs_skim spades -r "$clean" -o @ --mode meta
feat 09 genome-size "genomescope/{s}/{s}.genomescope.txt genomescope/genomescope_summary.tsv" \
    "genomescope/genomescope_summary.tsv：基因组大小（真值 ~20 Mb）、杂合度（真值 0.2 %）" \
    wgs_skim genome-size -r "$clean" -o @
links asm "$out/07_spades_genome/spades/{s}/{s}.scaffolds.fasta"
links mt "$out/05_meangs-getorganelle/mitogenome/{s}/{s}.mt.fasta"
feat 10 busco_genome "busco/{s}/{s}.busco.txt busco/busco_summary.tsv" "busco/busco_summary.tsv：完整度（20 Mb 片段只含少量基因，值低是正常的）" \
    wgs_skim busco -i "$out/inputs/asm" -o @
feat 11 mitom "mitom/{s}/{s}.mt.gb mitom/mitom_summary.tsv" "mitom/mitom_summary.tsv：13 PCG + 2 rRNA + 22 tRNA" \
    wgs_skim mitom -i "$out/inputs/mt" -o @
ref=$out/11_mitom/mitom/$first/$first.mt.gb
feat 12 blast "blast/blast_topn.tsv" "blast/blast_topn.tsv：cox1 / 16S 的最佳命中物种与相似度（本地 metazoa_marker 库）" \
    wgs_skim blast -i "$out/inputs/mt" -o @ --genes cox1,16s
feat 13 mitofinder_assembly "mitofinder_assembly/{s}/{s}.mitofinder.fasta mitofinder_assembly/mitofinder_summary.tsv" \
    "MitoFinder 在组装上，参考 = 1 个远缘种（$first）。mitofinder_assembly/：<样品>/ + mitofinder_summary.tsv（基因表）+ genes/ + contamination/" \
    wgs_skim mitofinder -i "$out/inputs/asm" -o @ --ref "$ref"
feat 14 mitofinder_reads "mitofinder_megahit/{s}/{s}.mitofinder.fasta mitofinder_megahit/mitofinder_summary.tsv" \
    "MitoFinder 自己组装 reads（MEGAHIT）：与 13 比较" \
    wgs_skim mitofinder -r "$clean" -o @ --ref "$ref"
feat 15 mitofinder_bioconda "mitofinder_assembly/{s}/{s}.mitofinder.fasta mitofinder_assembly/mitofinder_summary.tsv" \
    "bioconda MitoFinder 1.4.1（没有 Singularity 镜像时）：tRNA 用 tRNAscan-SE；与 13 比较" \
    wgs_skim mitofinder -i "$out/inputs/asm" -o @ --ref "$ref" --config mitofinder_sif=
if [ -d "$sim/refs" ]; then          # related species' mitogenomes (not the simulated ones): SIM_DIR/refs/*.gb
feat 16 mitofinder_multiref "mitofinder_assembly/{s}/{s}.mitofinder.fasta mitofinder_assembly/mitofinder_summary.tsv" \
    "与 13 相同，但参考 = 5 个近缘种（--ref 目录）：远缘样品的基因数与 13 比较（推荐用法）" \
    wgs_skim mitofinder -i "$out/inputs/asm" -o @ --ref "$sim/refs"
fi
feat 17 spades-busco "spades/{s}/{s}.scaffolds.fasta busco/{s}/{s}.busco.txt busco/busco_summary.tsv" "组合：SPAdes → BUSCO" \
    wgs_skim run -r "$clean" --no-qc -o @ --steps spades-busco
feat 18 spades-mitofinder "mitofinder_spades/{s}/{s}.mitofinder.fasta mitogenome/{s}/{s}.mt.fasta mitofinder_spades/mitofinder_summary.tsv" \
    "组合：GetOrganelle 线粒体 + SPAdes-MitoFinder（参考 = 每个样品自己的 MitoM 注释）" \
    wgs_skim run -r "$clean" --no-qc -o @ --steps meangs-getorganelle,spades-mitofinder
feat 19 circularize "{s}/{s}.circularize.tsv" "无 reads 环化：每个样品的 method（walk / overlap / join）与环长" \
    noj bash -c "for s in $(echo $samples); do wgs_skim circularize -i $out/18_spades-mitofinder/mitofinder_spades/\$s/\$s.mitofinder.fasta --spades $out/18_spades-mitofinder/spades/\$s/spades -o @/\$s --prefix \$s || exit 1; done" --
feat 20 compare "{s}/{s}.compare.tsv all.compare.html" "GetOrganelle 与 SPAdes-MitoFinder 是否一致：all.compare.tsv 的 verdict，<样品>/ 中差异与图" \
    wgs_skim compare --results "$out/18_spades-mitofinder" -o @
feat 21 run_all "summary/samples.tsv summary/report.html fastp/fastp_summary.tsv" \
    "完整流程（原始 reads，默认全部步骤）：summary/samples.tsv 每样品一行；各工具目录下有自己的 *_summary.tsv" \
    wgs_skim run -r "$dna" -o @
feat 22 nuclear "summary.tsv" "核多拷贝基因：GetOrganelle 与 SPAdes-MitoM 逐基因比较（verdict）" \
    wgs_skim nuclear --results "$out/21_run_all" -o @
printf "sample\tfastq_1\tfastq_2\tscaffolds\tmitogenome\n" > "$out/inputs/sheet_given.tsv"
for s in $samples; do printf "%s\t%s\t%s\t%s\t%s\n" "$s" "$clean/$s/$s.R1.fq.gz" "$clean/$s/$s.R2.fq.gz" \
    "$out/inputs/asm/$s.fasta" "$out/inputs/mt/$s.fasta" >> "$out/inputs/sheet_given.tsv"; done
feat 23 sheet_given "summary/samples.tsv mitom/{s}/{s}.mt.gb busco/busco_summary.tsv" \
    "样品表给出已有组装与线粒体：不重新组装，直接 BUSCO / MitoM / BLAST" \
    wgs_skim run -s "$out/inputs/sheet_given.tsv" --no-qc -o @ --steps busco,annotate,blast,summary
cat "$out"/inputs/mt/*.fasta > "$out/inputs/merged_mt.fasta"
feat 24 split_by_id "mitom/{s}*/*.mt.gb" "一个合并的 FASTA（多个样品）按序列名拆分后批量注释：mitom/<样品>/" \
    wgs_skim mitom -i "$out/inputs/merged_mt.fasta" --split-by-id -o @
feat 25 python_api "busco/{s}/{s}.busco.txt busco/busco_summary.tsv" "Python 接口 wgs_skim.tool(...)：与 10 相同" \
    noj "$py" -c "import sys, wgs_skim; sys.exit(wgs_skim.tool('busco', input=['$out/inputs/asm'], outdir='@', cores=$jobs, mem_gb=$mem))" --
if command -v nextflow > /dev/null; then
feat 26 nextflow "mitogenome/{s}/{s}.mt.fasta summary/samples.tsv" "Nextflow 引擎（原始 reads，默认步骤）：与 21 的汇总表比较" \
    wgs_skim run -r "$dna" -o @ --engine nextflow
fi
fi

if [ "$lane" = rna ] || [ "$lane" = all ]; then
echo "[features] RNA: $(echo $samples)" >&2
# MitoFinder reference of the RNA lane (runs beside the DNA lane, so not a DNA-lane output): the first sample's
# source GenBank (truth.tsv mt_source) - as a related species' mitogenome would be in a real study
ref=$(awk -F'\t' 'NR == 2 {for (i = 1; i <= NF; i++) if (h[i] == "mt_source") print $i} NR == 1 {for (i = 1; i <= NF; i++) h[i] = $i}' "$sim/truth.tsv")
feat 30 fastp_rna "fastp/{s}/{s}.R1.fq.gz fastp/fastp_summary.tsv" "RNA-seq 质控（只做一次，后面的 RNA 功能都用这里的 clean reads）" \
    wgs_skim fastp -r "$rna" -o @
rclean=$out/30_fastp_rna/fastp
feat 31 rnaspades "spades/{s}/{s}.rna.fasta spades/spades_summary.tsv" "rnaSPAdes：spades_summary.tsv 的转录本数、总长" \
    wgs_skim spades -r "$rclean" -o @ --mode rna
feat 32 trinity "trinity/{s}/{s}.trinity.fasta trinity/trinity_summary.tsv" "Trinity：trinity_summary.tsv 与 31 比较" \
    wgs_skim trinity -r "$rclean" -o @
links rna "$out/31_rnaspades/spades/{s}/{s}.rna.fasta"
feat 33 busco_transcriptome "busco/{s}/{s}.busco.txt busco/busco_summary.tsv" "BUSCO -m transcriptome（核转录本取自真实 RNA-seq 组装）" \
    wgs_skim busco -i "$out/inputs/rna" -o @ -m transcriptome
feat 34 rnaspades-busco "spades/{s}/{s}.rna.fasta busco/{s}/{s}.rna.busco.txt busco/busco_summary.tsv" "组合（推荐）：rnaSPAdes → BUSCO" \
    wgs_skim run -r "$rclean" --no-qc -o @ --steps rnaspades-busco
feat 35 trinity-busco "trinity/{s}/{s}.trinity.fasta busco/{s}/{s}.trinity.busco.txt busco/busco_summary.tsv" "组合：Trinity → BUSCO" \
    wgs_skim run -r "$rclean" --no-qc -o @ --steps trinity-busco
feat 36 rnaspades-mitofinder "mitofinder_rnaspades/{s}/{s}.mitofinder.fasta mitofinder_rnaspades/mitofinder_summary.tsv mitofinder_rnaspades/contamination/samples_check.tsv" \
    "转录组线粒体基因：mitofinder_summary.tsv（基因表，N = 未找到）；contamination/samples_check.tsv：S5 应报 COX1、CYTB 可能污染，其余不报" \
    wgs_skim run -r "$rclean" --no-qc -o @ --steps rnaspades-mitofinder --mitofinder-ref "$ref"
feat 37 trinity-mitofinder "mitofinder_trinity/{s}/{s}.mitofinder.fasta mitofinder_trinity/mitofinder_summary.tsv mitofinder_trinity/contamination/samples_check.tsv" \
    "同 36，Trinity 转录本" \
    wgs_skim run -r "$rclean" --no-qc -o @ --steps trinity-mitofinder --mitofinder-ref "$ref"
feat 38 mitofinder_rna_tool "mitofinder_rnaspades/mitofinder_summary.tsv mitofinder_trinity/mitofinder_summary.tsv" \
    "单独命令：wgs_skim mitofinder --assembler rnaspades,trinity（两种来源各一个目录）" \
    wgs_skim mitofinder -r "$rclean" -o @ --ref "$ref" --assembler rnaspades,trinity
fi

echo "[features] $lane: $fail failed -> $tsv" >&2
exit $fail

#!/bin/bash
# Every wgs_skim function on the simulated samples (tests/sim/simulate.py), one output directory per function:
#   tests/all_features.sh SIM_DIR OUT_DIR [dna|rna|all]
# SIM_DIR: simulate.py output (dna/, rna/, truth/, truth.tsv). OUT_DIR/NN_<function>/ holds that function's results
# and CHECK.md (command, result, time, the files to look at and what to check); OUT_DIR/features.tsv lists them all.
# The DNA and RNA functions are two lanes (run them at the same time: one call with dna, one with rna).
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
tsv=$out/features.$lane.tsv
echo -e "id\tfunction\tresult\tseconds\toutput\tcommand\tmissing" > "$tsv"
fail=0

# feat ID NAME "OUTPUT PATTERNS ({s} = each sample, relative to the output dir)" "WHAT TO CHECK" COMMAND...
# (@ in the command = this function's output directory)
feat() {
    local id=$1 name=$2 pats=$3 what=$4; shift 4
    local d="$out/${id}_${name}" t0=$(date +%s) ok=ok missing="" cmd=()
    for x in "$@"; do cmd+=("${x//@/$d}"); done
    mkdir -p "$d"
    echo "[$id] $name ..." >&2
    "${cmd[@]}" -j "$jobs" --mem-gb "$mem" > "$d/run.log" 2>&1 || { ok=FAIL; missing="exit $?"; }
    for s in $samples; do for p in $pats; do
        f="$d/${p//\{s\}/$s}"
        compgen -G "$f" > /dev/null && [ -s "$(compgen -G "$f" | head -1)" ] || { ok=FAIL; missing+=" ${p//\{s\}/$s}"; }
    done; done
    local sec=$(( $(date +%s) - t0 ))
    [ $ok = FAIL ] && fail=$((fail + 1))
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$id" "$name" "$ok" "$sec" "$d" "${cmd[*]}" "$missing" >> "$tsv"
    {
        echo "# $id $name: $ok ($sec s)"
        echo
        echo '```bash'
        echo "${cmd[*]}"
        echo '```'
        echo
        echo "**看什么**：$what"
        echo
        echo "**每个样品的输出**（应存在且非空）：\`${pats}\`"
        [ -n "$missing" ] && echo && echo "**缺少**：$missing"
        echo
        echo "运行日志：run.log；每个任务的日志：logs/；用时：benchmarks/"
    } > "$d/CHECK.md"
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
feat 01 samples "" "样品名识别是否正确（5 个，R1/R2 配对）" \
    bash -c "wgs_skim samples -r $dna > @/samples.tsv; cat @/samples.tsv" --
feat 02 fastp "qc/{s}/{s}.R1.fq.gz qc/{s}/{s}.fastp.html" "clean 数据量、Q30（模拟数据 ~0.45 Gb / 样品）" \
    wgs_skim fastp -r "$dna" -o @
clean=$out/02_fastp/qc
feat 03 meangs "mito/{s}/meangs/{s}.seed.fasta" "MEANGS 种子：长度、是否线粒体（与 truth/<s>/mt.fasta 比）" \
    wgs_skim meangs -r "$clean" -o @
feat 04 getorganelle_mt "mito/{s}/{s}.mt.json" \
    "GetOrganelle 自带种子（animal_mt）：长度与环化（truth: S1-S4 环状，S5 线状）；种子太远时可能没有结果（mt.json 的 stop），这正是推荐 05 的原因" \
    wgs_skim getorganelle -r "$clean" -o @
feat 05 meangs-getorganelle "mito/{s}/{s}.mt.fasta mito/{s}/{s}.mt.json mito/{s}/meangs/{s}.seed.fasta" \
    "推荐组合：MEANGS 种子 + GetOrganelle，未成环时迭代；与 04 比较" \
    wgs_skim run -r "$dna" -o @ --steps meangs-getorganelle
feat 06 getorganelle_nr "rdna/{s}/{s}.rdna.gbf histone/{s}/{s}.histone.gbf" \
    "rDNA 单元（18S-ITS1-5.8S-ITS2-28S）与组蛋白簇（H1-H4）：长度、基因是否完整（rounds.tsv）" \
    wgs_skim getorganelle -r "$clean" -o @ --target nr
feat 07 spades_genome "assembly/{s}/{s}.scaffolds.fasta" "SPAdes 基因组组装：总长 ~ truth nuclear_mb，N50" \
    wgs_skim spades -r "$clean" -o @
feat 08 spades_meta "assembly/{s}/{s}.meta.fasta" "metaSPAdes 模式：与 07 比较" \
    wgs_skim spades -r "$clean" -o @ --mode meta
feat 09 genome-size "genome_size/{s}/{s}.genomescope.txt" \
    "GenomeScope：基因组大小（truth ~20 Mb，另含线粒体/rDNA 高拷贝）、杂合度（truth 0.2 %）" \
    wgs_skim genome-size -r "$clean" -o @
links asm "$out/07_spades_genome/assembly/{s}/{s}.scaffolds.fasta"
links mt "$out/05_meangs-getorganelle/mito/{s}/{s}.mt.fasta"
feat 10 busco_genome "assembly/{s}/{s}.busco.txt" "BUSCO（metazoa_odb12）完整度：20 Mb 片段只含少量基因，值低是正常的" \
    wgs_skim busco -i "$out/inputs/asm" -o @
feat 11 mitom "annotation/{s}/{s}.mt.gb" "MitoM 注释：13 PCG + 2 rRNA + 22 tRNA（truth/<s>/mt_genes.fasta）" \
    wgs_skim mitom -i "$out/inputs/mt" -o @
ref=$out/11_mitom/annotation/$first/$first.mt.gb
feat 12 blast "blast/blast_topn.tsv" "cox1 / 16S 的最佳命中：物种与相似度（本地 metazoa_marker 库）" \
    wgs_skim blast -i "$out/inputs/mt" -o @ --genes cox1,16s
feat 13 mitofinder_assembly "mitofinder/given/{s}/{s}.mitofinder.fasta mitofinder/given/mitofinder_summary.tsv" \
    "MitoFinder 在组装上（参考 = $first 的 MitoM 注释）：基因表、genes/<基因>.fasta、污染检查 blast/samples_check.tsv" \
    wgs_skim mitofinder -i "$out/inputs/asm" -o @ --ref "$ref"
feat 14 mitofinder_reads "mitofinder/megahit/{s}/{s}.mitofinder.fasta mitofinder/megahit/mitofinder_summary.tsv" \
    "MitoFinder 自己组装 reads（MEGAHIT）：与 13 比较" \
    wgs_skim mitofinder -r "$clean" -o @ --ref "$ref"
feat 15 mitofinder_bioconda "mitofinder/given/{s}/{s}.mitofinder.fasta" \
    "bioconda MitoFinder 1.4.1（无镜像时）：tRNA 用 tRNAscan-SE；与 13 比较" \
    wgs_skim mitofinder -i "$out/inputs/asm" -o @ --ref "$ref" --config mitofinder_sif=
feat 16 spades-busco "assembly/{s}/{s}.scaffolds.fasta assembly/{s}/{s}.busco.txt" "组合：原始 reads → fastp → SPAdes → BUSCO" \
    wgs_skim run -r "$dna" -o @ --steps spades-busco
feat 17 spades-mitofinder "mitofinder/spades/{s}/{s}.mitofinder.fasta mito/{s}/{s}.mt.fasta" \
    "组合：GetOrganelle + SPAdes-MitoFinder（参考 = 每个样品自己的 MitoM 注释）" \
    wgs_skim run -r "$dna" -o @ --steps meangs-getorganelle,spades-mitofinder
feat 18 circularize "{s}.circularize.tsv" "无 reads 环化：MitoFinder contigs + SPAdes 组装图（method 列）" \
    noj bash -c "for s in $(echo $samples); do wgs_skim circularize -i $out/17_spades-mitofinder/mitofinder/spades/\$s/\$s.mitofinder.fasta --spades $out/17_spades-mitofinder/assembly/\$s/spades -o @ --prefix \$s || exit 1; done" --
feat 19 compare "{s}.compare.tsv all.compare.html" "GetOrganelle 与 SPAdes-MitoFinder 是否一致（verdict）、差异的位置与基因" \
    wgs_skim compare --results "$out/17_spades-mitofinder" -o @
feat 20 run_all "summary/samples.tsv summary/report.html" "完整流程（默认全部步骤）：汇总表每一列" \
    wgs_skim run -r "$dna" -o @
feat 21 nuclear "summary.tsv" "核多拷贝基因批量：GetOrganelle vs SPAdes-MitoM（verdict）" \
    wgs_skim nuclear --results "$out/20_run_all" -o @
printf "sample\tfastq_1\tfastq_2\tscaffolds\tmitogenome\n" > "$out/inputs/sheet_given.tsv"
for s in $samples; do printf "%s\t%s\t%s\t%s\t%s\n" "$s" "$dna/${s}_R1.fq.gz" "$dna/${s}_R2.fq.gz" \
    "$out/inputs/asm/$s.fasta" "$out/inputs/mt/$s.fasta" >> "$out/inputs/sheet_given.tsv"; done
feat 22 sheet_given "summary/samples.tsv annotation/{s}/{s}.mt.gb" \
    "样品表给出已有组装与线粒体：不重新组装，直接 BUSCO / MitoM / BLAST" \
    wgs_skim run -s "$out/inputs/sheet_given.tsv" -o @ --steps busco,annotate,blast,summary
cat "$out"/inputs/mt/*.fasta > "$out/inputs/merged_mt.fasta"
feat 23 split_by_id "annotation/{s}*/*.mt.gb" "一个合并的 FASTA（多个样品）按序列名拆分后批量注释" \
    wgs_skim mitom -i "$out/inputs/merged_mt.fasta" --split-by-id -o @
feat 24 python_api "assembly/{s}/{s}.busco.txt" "Python 接口 wgs_skim.tool(...) 与命令行结果相同（对比 10）" \
    noj "$py" -c "import sys, wgs_skim; sys.exit(wgs_skim.tool('busco', input=['$out/inputs/asm'], outdir='@', cores=$jobs, mem_gb=$mem))" --
if command -v nextflow > /dev/null; then
feat 25 nextflow "mito/{s}/{s}.mt.fasta summary/samples.tsv" "Nextflow 引擎（默认步骤）：与 20 的汇总表比较" \
    wgs_skim run -r "$dna" -o @ --engine nextflow
fi
fi

if [ "$lane" = rna ] || [ "$lane" = all ]; then
echo "[features] RNA: $(echo $samples)" >&2
# MitoFinder reference of the RNA lane (runs beside the DNA lane, so not a DNA-lane output): the first sample's
# source GenBank (truth.tsv mt_source) - as a related species' mitogenome would be in a real study
ref=$(awk -F'\t' 'NR == 2 {for (i = 1; i <= NF; i++) if (h[i] == "mt_source") print $i} NR == 1 {for (i = 1; i <= NF; i++) h[i] = $i}' "$sim/truth.tsv")
feat 31 rnaspades "assembly/{s}/{s}.rna.fasta" "rnaSPAdes 转录本数" \
    wgs_skim spades -r "$rna" -o @ --mode rna --qc
feat 32 trinity "assembly/{s}/{s}.trinity.fasta" "Trinity 转录本数；与 31 比较" \
    wgs_skim trinity -r "$rna" -o @ --qc
links rna "$out/31_rnaspades/assembly/{s}/{s}.rna.fasta"
feat 33 busco_transcriptome "assembly/{s}/{s}.busco.txt" "BUSCO -m transcriptome（核转录本取自真实 RNA-seq 组装）" \
    wgs_skim busco -i "$out/inputs/rna" -o @ -m transcriptome
feat 34 rnaspades-busco "assembly/{s}/{s}.rna.fasta assembly/{s}/{s}.rna.busco.txt" "组合（推荐）：rnaSPAdes → BUSCO" \
    wgs_skim run -r "$rna" -o @ --steps rnaspades-busco
feat 35 trinity-busco "assembly/{s}/{s}.trinity.fasta assembly/{s}/{s}.trinity.busco.txt" "组合：Trinity → BUSCO" \
    wgs_skim run -r "$rna" -o @ --steps trinity-busco
feat 36 rnaspades-mitofinder "mitofinder/rna/{s}/{s}.mitofinder.fasta mitofinder/rna/mitofinder_summary.tsv mitofinder/rna/blast/samples_check.tsv" \
    "转录组线粒体基因：基因表（15 个，N = 未找到）；污染检查：S5 应报 COX1、CYTB 可能污染（Nephtyidae），其余不报" \
    wgs_skim run -r "$rna" -o @ --steps rnaspades-mitofinder --mitofinder-ref "$ref"
feat 37 trinity-mitofinder "mitofinder/trinity/{s}/{s}.mitofinder.fasta mitofinder/trinity/mitofinder_summary.tsv mitofinder/trinity/blast/samples_check.tsv" \
    "同 36，Trinity 转录本" \
    wgs_skim run -r "$rna" -o @ --steps trinity-mitofinder --mitofinder-ref "$ref"
feat 38 mitofinder_rna_tool "mitofinder/rna/mitofinder_summary.tsv mitofinder/trinity/mitofinder_summary.tsv" \
    "单独命令：wgs_skim mitofinder --assembler rnaspades,trinity（两种来源各一套结果）" \
    wgs_skim mitofinder -r "$rna" -o @ --ref "$ref" --assembler rnaspades,trinity --qc
fi

echo "[features] $lane: $fail failed -> $tsv" >&2
exit $fail

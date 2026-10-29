# 输入数据

## 原始数据目录 `-r`

双端 FASTQ（可 gzip），放在一个目录里，或每个样品一个子目录。样品名取文件名中读段编号之前的部分。可识别的命名：

| 例子 | 来源 | 样品名 |
|---|---|---|
| `S_R1_001.fastq.gz` / `S_R2_001.fastq.gz` | Illumina | S |
| `S_S1_L001_R1_001.fastq.gz` | Illumina bcl2fastq | S |
| `S.R1.raw.fastq.gz` | 诺禾 / 美吉 | S |
| `S_1.fq.gz`、`S_1.clean.fq.gz` | 华大 | S |
| `S_R1_trimmed.fq.gz` | Macrogen | S |
| `S/S_1.fq.gz` | 每样品一个子目录 | S |

先检查：`wgs_skim samples -r DIR`。只有 R1 没有 R2 的文件会被跳过并提示；同一样品出现两次（如两个 lane）会报错。

## 样品表 `-s samples.tsv`

命名不规则、要合并多个目录、或只跑部分样品时用。TSV 或 CSV，相对路径相对于表所在目录：

```text
sample	fastq_1	fastq_2	scaffolds
CHY006	raw/CHY006_1.fq.gz	raw/CHY006_2.fq.gz
CHY007	/data/x/CHY007_R1.fq.gz	/data/x/CHY007_R2.fq.gz	/data/old/CHY007/scaffolds.fasta
```

- 可选列 `scaffolds`：已有的组装，填了就不再组装（直接用于 BUSCO、MitoFinder、核基因、线粒体兜底）。
- 可选列 `mitogenome`：已有的线粒体序列，填了就不再组装（直接用于 MitoM 注释、BLAST）。
- 只有组装 / 序列、没有 reads 的样品：`fastq_1`、`fastq_2` 留空（只能运行 FASTA 类工具）。
- Excel 导出的 "CSV UTF-8" 可以直接用。

## FASTA 输入 `-i`（BUSCO、MitoFinder、MitoM、BLAST）

```bash
wgs_skim busco -i assemblies/                 # 目录：每个 .fasta / .fa / .fna / .fas（可 .gz）= 一个样品，样品名 = 文件名
wgs_skim busco -i A.fasta B.fasta             # 列出文件
wgs_skim busco -i assemblies/ --only A,B      # 只要部分样品
wgs_skim mitofinder -i all.fasta --split-by-id --ref ref.gb   # 一个文件含多个样品：按序列名第一列拆分
```

`--split-by-id`：`>B2 topology=circular`、`>B17 …` 这样的合并文件按第一列（B2、B17）拆成每个样品一个文件，
同一样品的多条序列全部保留（改名为 `B2_contig1`、`B2_contig2` …），`<输出>/.wgs_skim/split/split_manifest.tsv`
记录每条序列的原始名称。

## reads 是否已清洗

`wgs_skim run`（组合流程）总是先运行 fastp。单个 reads 类工具（`spades`、`trinity`、`getorganelle`、`meangs`、`mitofinder -r`、
`genome-size`）默认直接用给出的 reads（视为已清洗，fastp 只统计数据量），加 `--qc` 先清洗。

## 已有结果目录（`wgs_skim nuclear` / `wgs_skim compare`）

两种结果目录可以直接作为输入：本软件 `wgs_skim run` 的结果，或旧流程 Temp_scripts/wgs_skim_pipeline 的结果
（`01.fastp`、`02.organelle/03_nr`、`03.spades`、`10.mitofinder` …）。软件自动识别目录结构。

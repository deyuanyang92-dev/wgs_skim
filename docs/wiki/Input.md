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

- 可选列 `scaffolds`：已有的 SPAdes 组装，填了就不再重新组装（直接用于 BUSCO、核基因、线粒体兜底）。
- Excel 导出的 "CSV UTF-8" 可以直接用。

## 已有结果目录（`wgs_skim nuclear` / `wgs_skim compare`）

两种结果目录可以直接作为输入：本软件 `wgs_skim run` 的结果，或旧流程 Temp_scripts/wgs_skim_pipeline 的结果
（`01.fastp`、`02.organelle/03_nr`、`03.spades`、`10.mitofinder` …）。软件自动识别目录结构。

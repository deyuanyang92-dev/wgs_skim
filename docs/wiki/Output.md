# 输出

结果目录按工具分（与 nf-core、Snakemake 社区流程的习惯相同）：**每个工具一个目录，其中每个样品一个子目录，
该工具的汇总表在工具目录这一层**；所有工具合在一起的汇总在 `summary/`；日志、配置等运行记录在 `pipeline_info/`。
只运行一个工具（如 `wgs_skim busco`）时，结果目录里只有该工具的目录和 `pipeline_info/`。

```
res/
├── fastp/                    <s>/  clean reads（<s>.R1.fq.gz / R2）、<s>.fastp.html、<s>.fastp.json
│                             fastp_summary.tsv      每个样品：原始 / clean 数据量、Q30、GC、去接头的 reads
├── meangs/                   <s>/<s>.seed.fasta     MEANGS 线粒体种子；meangs_summary.tsv
├── getorganelle/             <s>/  round0、round1…（每轮 GetOrganelle）、<s>.getorganelle.fasta / .json（选中的一轮）
│                             getorganelle_summary.tsv   长度、是否成环、用了几轮、停止原因
├── mitogenome/               <s>/<s>.mt.fasta       最终线粒体（GetOrganelle 或 SPAdes 中更好的一个）；<s>.mt.json
│                             mitogenome_summary.tsv
├── mitom/                    <s>/<s>.mt.gb          MitoM 注释（GenBank）、<s>.mitom.tsv；mitom_summary.tsv
├── topology/                 <s>/                   reads 检查环 / 线性、末端诊断
├── rdna/                     <s>/  <s>.nr.fasta、<s>.rdna.gbf（18S ITS1 5.8S ITS2 28S）、nr_features.tsv、rounds.tsv
│                             rdna_summary.tsv
├── histone/                  <s>/  <s>.histone.fasta、<s>.histone.gbf（H1 H2A H2B H3 H4）；histone_summary.tsv
├── nuclear_markers/          <s>/  SPAdes 中的 rDNA / 组蛋白拷贝、拷贝数、与 GetOrganelle 的比较
├── spades/                   <s>/  <s>.scaffolds.fasta（基因组）、<s>.<模式>.fasta（rna、meta…）、spades/（SPAdes 目录）
│                             spades_summary.tsv     每个样品每个组装：条数、总长、N50、最长
├── trinity/                  <s>/<s>.trinity.fasta；trinity_summary.tsv
├── busco/                    <s>/<s>.busco.txt（基因组）、<s>.<组装>.busco.txt（rna、trinity…）
│                             busco_summary.tsv      C / S / D / F / M
├── genomescope/              <s>/  k-mer 直方图、GenomeScope2 图与结果；genomescope_summary.tsv
├── mitofinder_<来源>/        <s>/  MitoFinder 结果（来源：spades、rnaspades、trinity、megahit、assembly = 用户给的组装）
│                             mitofinder_summary.tsv（基因表）、genes/<基因>.fasta、contigs/、contamination/（污染检查）
├── blast/blast_topn.tsv      线粒体 / rDNA 的前 N 个命中，带 GenBank 物种与分类
├── summary/                  samples.tsv（所有工具，每个样品一行，列见下）、report.html
└── pipeline_info/            运行记录（不是结果）
    ├── logs/<任务>/<s>.log         每个任务的工具输出
    ├── benchmarks/<任务>/<s>.tsv   每个任务的时间与峰值内存
    ├── config_<时间>.yaml          这次运行实际使用的全部参数
    ├── run_<时间>.txt              命令、wgs_skim 与 Snakemake 版本、样品
    ├── read_stats/<s>/             reads 不清洗时（--no-qc / 单个工具）的 reads 数，用于估计内存
    └── inputs/                     命令行准备的输入：-i FASTA 的样品表、拆分的 FASTA、合并的 MitoFinder 参考
```

**质控只做一次**：`wgs_skim run` 默认先 fastp，后面所有步骤直接读取 `fastp/<s>/` 中的 clean reads；reads 已经是 clean 的
（例如上一次运行的 `fastp/` 目录），用 `wgs_skim run --no-qc` 或单个工具（默认不质控），工具直接读取原来的文件，
不会再生成 `fastp/`。

## summary/samples.tsv 的列

| 列 | 含义 | 例（cwsc010，测试数据 0.3 Gb） |
|---|---|---|
| `sample` | 样品名 | cwsc010 |
| `clean_gb`、`q30` | fastp 后的数据量（Gb）、Q30 % | 0.30、96.5 |
| `mt_length`、`mt_circular` | 线粒体长度、是否成环 | 8954、no |
| `mt_source` | 来源：GetOrganelle 第几轮（停止原因）或 SPAdes | GetOrganelle round1 (rounds exhausted) |
| `mt_genes` | MitoM 找到的基因（完整的动物线粒体：13 PCG / 22 tRNA / 2 rRNA） | 9 PCG / 13 tRNA / 2 rRNA |
| `mt_topology`、`mt_end_diagnosis` | reads 检查的结论（每条序列一个），末端为什么没接上（见 [线粒体](Combo-mitogenome.md)） | linear_low_reads… |
| `nr_length`、`nr_topology` | rDNA 单元长度、拓扑 | 8004、circular |
| `spades_n50`、`spades_total_mb`、`busco` | SPAdes N50、总长、BUSCO | 274、44.3、C:0.0% |
| `genome_size_mb`、`heterozygosity_pct` | GenomeScope2 | 0.3、3.85（数据太少，不可信） |
| `rdna_copies` | rDNA 拷贝数（需要 BUSCO 单拷贝基因作标准） | 空 |
| `histones` | 组蛋白基因与拷贝数（`x?`：无法估计） | H1 x? (5'-partial), H2A x? … |
| `blast_mt`、`blast_nr` | 每个查询的最佳命中：物种（一致度、覆盖度） | cox1: Nephtys sp. (93.6%, qcov 100%) |

没有运行的步骤，对应的列为空。

## 旧流程目录里的零散文件

旧流程（Temp_scripts/wgs_skim_pipeline）在每个结果目录中写了 `returncode_*.txt`、`stderr_*.log`、`stdout_*.log`、`cmd_*.txt` 等。
本软件不删除、不改动旧结果：这些文件在旧流程运行完之后统一移到各步骤目录下的 `_run_logs/`（`batch_report.py` 运行时还会读取它们，
所以旧流程运行中不移动），移动记录写进项目的 `CLEANUP.md`。本软件自己的日志只在 `pipeline_info/logs/` 和 `pipeline_info/benchmarks/` 中。

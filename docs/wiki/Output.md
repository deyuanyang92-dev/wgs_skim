# 输出

```
res/
├── summary/samples.tsv        每个样品一行（列见下）
│   summary/report.html        同一内容的网页报告
├── qc/<s>/                    clean reads（<s>.R1.fq.gz / R2）、fastp.html、fastp.json
├── mito/<s>/<s>.mt.fasta      最终线粒体；<s>.mt.json：来源、是否成环、迭代停止原因；round0..N 每轮 GetOrganelle
│   mito/<s>/<s>.getorganelle.json   GetOrganelle 选中的那一轮
├── annotation/<s>/<s>.mt.gb   MitoM 注释（GenBank）；<s>.mitom.tsv 基因数
├── topology/<s>/              reads 检查环 / 线性、末端诊断
├── rdna/<s>/                  <s>.nr.fasta、<s>.rdna.gbf（18S ITS1 5.8S ITS2 28S）、nr_features.tsv、rounds.tsv
├── histone/<s>/               <s>.histone.fasta、<s>.histone.gbf（H1 H2A H2B H3 H4，组蛋白类别）、rounds.tsv
├── nuclear/<s>/               SPAdes 中的 rDNA / 组蛋白拷贝、consensus/、copy_number.tsv、nuc_compare.tsv
├── assembly/<s>/              <s>.scaffolds.fasta、<s>.busco.txt；其他模式 <s>.<mode>.fasta + .<mode>.busco.txt；
│                              Trinity <s>.trinity.fasta + .trinity.busco.txt
├── genome_size/<s>/           k-mer 直方图、GenomeScope2 图和 summary.txt
├── blast/blast_topn.tsv       COX1 / rDNA 前 N 个命中，带 GenBank 物种与分类
├── logs/<任务>/<s>.log        每个任务的工具输出
└── benchmarks/<任务>/<s>.tsv  每个任务的时间与峰值内存
```

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
所以旧流程运行中不移动），移动记录写进项目的 `CLEANUP.md`。本软件自己的日志只在 `logs/` 和 `benchmarks/` 中。

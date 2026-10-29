# 步骤详解

`--steps` 选择要跑的步骤（逗号分隔），它们依赖的上游步骤会自动一起运行，已完成的结果直接复用。不写 `--steps` = 全部。
时间与内存为实测值：双样品测试数据（每个样品 1 M 对 reads = 0.3 Gb，见 [[测试|Testing]]），128 核机器；真实数据
（10–50 Gb / 样品）的时间随数据量增加，内存按数据量自动估计（见 [[资源调度与续跑|Resources-and-resume]]）。

| 步骤 | 做什么 | 主要输出 | 测试数据：时间 / 内存峰值（每样品） |
|---|---|---|---|
| `fastp` | 质控、去接头 | `qc/<s>/<s>.R1.fq.gz`、`.fastp.html` | 30 s / 0.4 GB |
| `meangs` | MEANGS 找线粒体种子 | `mito/<s>/meangs/<s>.seed.fasta` | 4–5 min / 5.9 GB |
| `mito` | MEANGS → GetOrganelle → 未成环则迭代 → SPAdes 兜底 | `mito/<s>/<s>.mt.fasta` + `.mt.json` | GetOrganelle 2–3 min / 3.6 GB；迭代 2–10 min / 5.8 GB |
| `getorganelle` | 同上，但不用 MEANGS 种子 | 同上 | 同上 |
| `annotate` | MitoM 注释线粒体 | `annotation/<s>/<s>.mt.gb` | 15 s – 2 min / 0.3 GB |
| `topology` | 用 reads 验证环 / 线性、末端诊断 | `topology/<s>/<s>.topology.tsv` | 5 s / 0.1 GB |
| `nr` | GetOrganelle 组装 rDNA 与组蛋白基因簇，MitoM 注释 | `rdna/<s>/<s>.rdna.gbf`、`histone/<s>/<s>.histone.gbf` | rDNA 3–9 min / 0.9 GB；组蛋白 10–11 min / 2.4 GB |
| `spades` | SPAdes 组装（默认 genome；`--spades-modes` 选其他模式） | `assembly/<s>/<s>.scaffolds.fasta`、`<s>.<模式>.fasta` | 6–8 min / 2.8 GB |
| `busco` | BUSCO 评估每个组装 | `assembly/<s>/<s>.busco.txt`、`<s>.<模式>.busco.txt` | 40 s / 1.7 GB |
| `trinity` | Trinity 转录组组装 | `assembly/<s>/<s>.trinity.fasta` | 见 [[SPAdes 各模式与转录组|Assembly-modes]] |
| `rna` | 转录组：rnaSPAdes + Trinity + 两者的 BUSCO（transcriptome 模式） | `<s>.rna.fasta`、`<s>.trinity.fasta`、`*.busco.txt` | 同上 |
| `genome_size` | jellyfish + GenomeScope2 | `genome_size/<s>/<s>.genomescope.txt` + 图 | 30 s – 2 min / 6.8 GB |
| `nuclear` | SPAdes 组装中的 rDNA / 组蛋白、拷贝数、与 GetOrganelle 结果逐基因比较 | `nuclear/<s>/<s>.copy_number.tsv`、`.nuc_compare.tsv` | 15–50 s / 1.4 GB |
| `blast` | COX1 与 rDNA 比对 nt，Top 命中及 GenBank 分类 | `blast/blast_topn.tsv` | 本地小库 2.5 min（全批次） |
| `summary` | 汇总表与 HTML 报告（单独运行时汇总已有的全部结果） | `summary/samples.tsv`、`report.html` | 1 s |

别名：`qc`=fastp，`mitogenome` / `meangs-getorganelle`=mito，`rdna` / `histone` / `histones`=nr，`assembly`=spades，
`jellyfish` / `genomescope`=genome_size，`mitom` / `mitoz`=annotate，`nuc`=nuclear，`rnaspades`=rna，`report`=summary。

## 每个步骤的批量示例

```bash
# 质控（fastp）
wgs_skim run -r raw/ -o res/ --steps fastp

# SPAdes 全基因组组装，再 BUSCO
wgs_skim run -r raw/ -o res/ --steps spades
wgs_skim run -r raw/ -o res/ --steps busco                     # 会先跑 spades（已完成则复用）
wgs_skim run -r raw/ -o res/ --steps spades,busco -- --set-threads spades=48 --set-resources spades:mem_mb=250000

# 线粒体：MEANGS-GetOrganelle（默认）或只用 GetOrganelle
wgs_skim run -r raw/ -o res/ --steps mito
wgs_skim run -r raw/ -o res/ --steps getorganelle
wgs_skim run -r raw/ -o res/ --steps annotate,topology         # + MitoM 注释 + 环 / 线性验证

# 核多拷贝基因
wgs_skim run -r raw/ -o res/ --steps nr                        # rDNA + 组蛋白：GetOrganelle + MitoM .gbf
wgs_skim run -r raw/ -o res/ --steps nuclear                   # + SPAdes 中的同一批基因、拷贝数、两种方法比较

# BLAST 物种鉴定
wgs_skim run -r raw/ -o res/ --steps blast --blast-db /db/metazoa_marker/metazoa_marker

# 基因组大小、杂合度
wgs_skim run -r raw/ -o res/ --steps genome_size

# 只重新生成报告
wgs_skim run -r raw/ -o res/ --steps summary
```

多个步骤：`--steps fastp,spades,busco`。只跑部分样品：用样品表。修改参数后重跑某一步：`-- --forcerun busco`。
干运行看会跑哪些任务：加 `-n`。

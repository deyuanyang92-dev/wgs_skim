# MitoFinder：从组装中找出线粒体

MitoFinder 用参考线粒体（GenBank）在一个组装的全部 contig 中找出线粒体 contig，并注释基因。可以搜索的组装：

| 来源（`--assembler` / `mitofinder_assembly`） | 是什么 | 适合 |
|---|---|---|
| `-i` 给出的 FASTA（`given`） | 已有的任何组装 | 以前跑过的 SPAdes / MEGAHIT / Trinity 结果 |
| `spades` | SPAdes 全基因组组装 | 基因组浅测序（组合 spades-mitofinder） |
| `rnaspades` | rnaSPAdes 转录本 | **转录组中获取线粒体基因**（组合 rnaspades-mitofinder） |
| `trinity` | Trinity 转录本 | 同上（组合 trinity-mitofinder） |
| `megahit`（默认）/ `metaspades` | MitoFinder 自己组装 reads（前 1000 万对） | 只有 reads、不需要其他组装 |

```bash
# 已有组装（目录中每个 FASTA 一个样品）
wgs_skim mitofinder -i assemblies/ -o mf/ --ref related.gb
# 合并在一个 FASTA 中的多个样品（>B2 ...、>B17 ...）：按序列名第一列拆分
wgs_skim mitofinder -i all_contigs.fasta --split-by-id -o mf/ --ref related.gb
# 从 RNA-seq 开始：rnaSPAdes 和 / 或 Trinity 组装，再 MitoFinder（每种组装一次）
wgs_skim mitofinder -r rnaseq_clean/ -o mf/ --ref related.gb --assembler rnaspades,trinity
# 从 DNA reads 开始：MitoFinder 用 MEGAHIT 组装
wgs_skim mitofinder -r clean_reads/ -o mf/ --ref related.gb
```

## 参考 `--ref`

近缘物种的线粒体 GenBank，必须带基因注释。`--ref` 可给多个文件或一个目录，自动合并成一个多记录 GenBank
（MitoFinder 每个基因选最相近的参考）。**建议给几个近缘种**：模拟测试（sim5）中，参考是另一个科的一条线粒体时，
远缘样品只找到 12–13 / 15 个基因（缺 ND2、ND4L、rrnL）；换成 5 个近缘种（不含被测物种本身）→ 15 / 15。
在组合 `meangs-getorganelle,spades-mitofinder` 中默认用该样品自己的 MitoM 注释作参考。

## 送给 MitoFinder 之前的过滤

- 全基因组组装（spades / given）：只保留 ≥ 1000 bp、k-mer 覆盖度 ≥ 5 × 长度加权中位覆盖度的 scaffold
  （线粒体比核基因组深几十到几百倍；来自旧脚本 spades_mitofinder.py：CHY006 的 2 Gb 组装 → 3 871 条，MitoFinder 5 分钟、
  1 条线粒体 contig，环化检查才会运行）。没有 SPAdes 覆盖度信息的组装只按长度过滤。`mitofinder_cov_factor`（0 = 不按覆盖度）。
- 低覆盖度：参考序列 BLAST 命中的 scaffold 不论覆盖度都保留（小数据量时线粒体 scaffold 可能低于覆盖度阈值：
  cwsc010 子样本的线粒体只有 1.5–2.9×，阈值 6.3×，原来全部被滤掉，MitoFinder 什么也找不到）。
- 转录本（rnaspades / trinity）：长度 ≥ `mitofinder_min_contig_rna`（默认 300 bp），并且**把候选线粒体转录本拼成一条序列**
  （`mitofinder_concat_transcripts`，默认开）：MitoFinder 用核苷酸 BLAST 判断哪些 contig 是线粒体，一条完整线粒体靠几个保守基因
  就能通过，然后整条注释；但每个转录本只有一个基因，必须单独通过，参考较远时大多通不过（sim5 S3：15 个基因只找到 2 个）。
  wgs_skim 先用 blastx（参考的蛋白）和 dc-megablast 挑出候选转录本，以 200 个 N 隔开连成一条，交给 MitoFinder
  （S3：2 → 12 个基因；其他样品 12–14 个），再把每个基因映射回它所在的转录本。MitoFinder 每条序列每个基因只注释一次，
  同一基因的其他转录本（另一个异构体，或污染物的拷贝）作为额外拷贝保留（`<s>.mitofinder.extra_copies.fasta`，
  汇总中为 `<样品>_copyN`），污染检查因此仍能看到它们。`--config mitofinder_concat_transcripts=false` = 旧方式（逐条）。
  ATP8（约 160 bp）短于 300 bp，转录本中一般找不到。

## 实际命令与软件

`mitofinder -j S -a input.fasta -r ref.gb -o <遗传密码> -p 8 -m <内存> -t mitfi --min-contig-size <N> --adjust-direction`；
`--trna mitfi | trnascan | arwen`；遗传密码来自 `--taxon`（环节动物 5）。

- 配置 `mitofinder_sif` 时用 Singularity 镜像 MitoFinder **1.4.2**（`singularity exec --cleanenv`；本机 `~/mybiosofts/mitofinder_v1.4.2.sif`）。
- 否则自动建 bioconda 环境 MitoFinder **1.4.1**：它的环化检查实际没有运行（启动检查脚本时没有传参数），结果中的
  "Circularization" 不可信，汇总 README.txt 会注明；环化请以 `wgs_skim circularize` / topology 为准。

## 输出（按来源分目录，`mitofinder_<来源>/`）

| 文件 | 内容 |
|---|---|
| `<s>/<s>.mitofinder.fasta`、`.gb` | 该样品找到的全部线粒体 contig 及注释（转录本：带基因的转录本） |
| `<s>/<s>.mitofinder.transcripts.tsv` | （转录本）每个注释的基因在哪条转录本、位置 |
| `<s>/<s>.mitofinder.extra_copies.fasta` | （转录本）MitoFinder 未注释的同一基因的其他转录本 |
| `<s>/Final_Results/` | MitoFinder 原始结果（基因 NT / AA、.infos、.png、tRNA） |
| `mitofinder_summary.tsv` | 每个样品：contig 数、总长、是否环化、蛋白基因 / rRNA / tRNA 数、每个基因的长度（N = 未找到） |
| `mitofinder_contigs.tsv` | 每条线粒体 contig：原始 contig 名、长度、GC、环化 |
| `genes/<基因>.fasta` | 同一基因的全部样品（多拷贝 `<样品>_copy2`…，名称唯一），可直接比对、建树 |
| `contigs/<s>_all_mt.fasta` | 每个样品的线粒体 contig（+ ID 列表） |
| `README.txt` | 版本、来源、各文件说明 |
| `<s>/<s>.megahit.contigs.fasta` | （从 reads 组装时）MitoFinder 的组装 |

同一批样品用不同来源（例如 rnaSPAdes 与 Trinity）的结果在不同目录中，可以直接对比。

## 污染检查（BLAST，默认开启）

转录组中提取的线粒体基因可能来自污染（肠道内容物、寄生虫、文库交叉污染）。汇总之后，每个样品的每个基因都比对 nt
（与 BLAST 步骤相同：blast2metadata 的 Top N 命中 + GenBank 物种、分类），比较同一样品各基因的最佳命中：

| 判断 | 条件 |
|---|---|
| `ok` | 最佳命中与该样品多数基因同科，一致度 ≥ 85 %（`mitofinder_min_pident`） |
| `no close reference` | 同科（或别的科）但一致度低：数据库中没有近缘序列，不是污染 |
| `partial hit` | 命中只覆盖基因 < 50 %，不计入判断 |
| `possible contamination` | 最佳命中是**别的科**且一致度 ≥ 85 %：污染物与其来源高度相似 |
| `possible contamination`（同种规则） | 某拷贝与数据库中某物种 X 几乎相同（≥ 97 %），并且：同一基因的另一个拷贝不是 X（< 90 %）——同一样品有两种不同的 COX1 就是两种生物；或 X 的最佳记录是完整线粒体，而样品其他基因都不与它相似。样品中这些拷贝之外没有任何基因支持 X 时才报（候选拷贝之间不互相作证）。不受多数投票影响：sim5 S5 中污染的 COX1、CYTB 两个基因曾以 2 票把样品判为 Nephtyidae。近缘属（< 97 %）不报 |

输出 `mitofinder_<来源>/blast/`：`samples_check.tsv`（每个样品：多数科与物种、一致的基因数、被标记的基因、结论）、
`genes_check.tsv`（每个样品每个基因：最佳命中、物种、科、一致度、覆盖度、结论）、`blast_topn.tsv`（全部命中）。

实测（MJ20260820041 CHY006–008 的 MitoFinder 结果，本地后生动物库，45 个基因 105 s）：nt 中没有近缘种（COX1 最佳 82 %
*Travisia*，Opheliidae），部分基因以 72–83 % 命中果蝇、鱼、芋螺等 → 判为 "no close reference"，不误报；把一个样品的 COX1
换成 *Nephtys*（另一科）的 COX1 → 判为 `possible contamination: COX1 (Nephtyidae, 93.6%)`。

模拟测试（sim5，`tests/sim/`）：S5 的 RNA 中加入 Nephtys（EU293739）的 COX1、CYTB → 报 `possible contamination: COX1, CYTB`（rnaSPAdes 与 Trinity 都是）；同一样品的 DNA 结果不报；
其余 4 个样品不报（S4 的 rrnS 88.7 % 命中 *Micronephthys*，同科近缘属，不是同一物种，不报）。

关闭：`--no-check`（工具）或 `--config mitofinder_blast=false`；只查部分基因：`--config mitofinder_blast_genes=COX1,CYTB`。

MitoFinder 的线粒体常常未闭合：用 `wgs_skim circularize`（SPAdes 组装图，无需 reads）尝试闭合，用 `wgs_skim compare`
与 GetOrganelle 结果比较，见 [线粒体组合](Combo-mitogenome.md)。

深入资料（版本差异、参数、已知 issues）：[docs/tools/MitoFinder.md](../tools/MitoFinder.md)。

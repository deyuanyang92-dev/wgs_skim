# wgs_skim 测试记录

最后更新：2026-10-29。代码版本：本文件所在提交（GitHub `main`）。MitoM 版本：`6f75f73`（`config/config.yaml` 的 `mitom_nuc_commit`）。

本文记录：测了哪些功能、怎么测、用什么数据、结果如何、测试中发现了哪些错误以及如何修复，最后列出**尚未测试**的部分。
自己重复测试的方法见第 6 节。

---

## 1. 测试数据

| 名称 | 内容 | 大小 | 用途 |
|---|---|---|---|
| HNP1sub | 一个多毛类（*Pisione*）样品的 2 M 对 reads（MJ20260820041 `smk_test/00.rawdata`） | 0.6 Gb | 单样品完整流程、各步骤 |
| data2（双样品） | `pisione`（HNP1sub 前 1 M 对）+ `cwsc010`（CWSC010_0729 前 1 M 对），两个不同物种 | 各 0.3 Gb | 批量流程（多样品）、快速回归测试 |
| CHY006 合成数据 | CHY006 线粒体参考上人为制造的错误（见 3.3） | — | `compare` 结构检查 |
| MJ20260820041 | 22 个环节动物样品的真实结果（旧流程 Temp_scripts 的结果目录） | 每样品 10–50 Gb | 真实批量：核基因、比较 |
| 单元测试 | `tests/test_scripts.py` 中的人工序列 | — | 各脚本函数 |

---

## 2. 单元测试（`python -m pytest -q tests/`，16 项，0.3 秒，全部通过）

不需要任何外部软件，用人工序列检查各函数的输入输出。

| 测试 | 检查的功能 | 方法 / 期望 |
|---|---|---|
| `test_find_pairs_recognises_provider_names` | 样品名自动识别（`samples.py`） | 6 种测序公司命名（`_R1_001`、`_1`、`.R1.raw`、Illumina `_S3_L001_`、`_trimmed`、子目录 `_1.clean`）全部配对正确；无 R2 的文件和无关文件被跳过 |
| `test_find_pairs_duplicate_sample_raises` | 同一样品两种命名 | 报错 "found twice"，不会静默覆盖 |
| `test_sheet_relative_paths_and_optional_columns` | 样品表（sample sheet） | 相对路径按表所在目录解析；可选列 `scaffolds`（已有组装）只对填写的样品生效 |
| `test_fasta_roundtrip_and_circle_detection` | FASTA 读写、首尾重叠识别环（`seqtools.find_circle`） | 3000 bp 环 + 400 bp 首尾重叠 → 识别 P=3000、overlap=400；随机线性序列 → 不判为环 |
| `test_collect_fasta_prefers_complete_graph` | GetOrganelle 结果选择（`collect_fasta.py`） | 同时有 complete 与 scaffolds 结果时选 complete；重命名格式 `S1\|mt\|1 circular length=8` |
| `test_summary_helpers` | 报告汇总（`summary.py`） | N50、FASTA 长度、BUSCO 行、GenomeScope 基因组大小 / 杂合度的解析 |
| `test_split_point_follows_rotation_and_strand` | 环化验证的切点（`topology.py`） | 注释记录被旋转 700 bp、或在反链时，nad5 中点仍映射到正确位置 |
| `test_end_diagnosis_classes` | 线性末端诊断（`topology.py`） | 干净末端 / 重复序列（AT)n 延伸 / 同聚物 / 覆盖度下降 四类判断 |
| `test_circularize_overlap_and_join` | 无 reads 环化（`circularize.py`） | 首尾重叠 300 bp 的 contig → 环 4000 bp；三段（其中一段反向，如 MitoFinder 输出）→ 拼成 4000 bp 环；随机序列 → 不成环 |
| `test_diff_blocks_and_description` | 两组装差异（`compare_mt.py`） | 1 bp 插入与 3 bp 缺失被正确定位；落在 cox1 内 / 基因间区（cox1..nad2）/ GetOrganelle 重复区的标注 |
| `test_orient_turns_and_rotates` | 方向与起点统一 | 反向且旋转 1000 bp 的序列 → 转回参考方向和起点 |
| `test_homopolymer_difference_is_explained` | 同聚物长度差异 | T×12 中多一个 T → 判为 "homopolymer"（二代数据无法确定，归为一致） |
| `test_gff_cleaning_as_modify_gff2` | GFF 注释读取（按用户 `modify_gff2.py` 规则） | 去掉 region / gene 行；`Name=COX1 gene` → COX1；`trnN%28guu%29` → trnN |
| `test_svg_figure_is_well_formed` | 比较图（`plot_compare.py`） | 生成的 SVG 是合法 XML，包含基因与反向比对图例 |
| `test_reads_cap_and_completeness` | 核基因数据量上限与完整性判断（`nuclear_markers.py`） | 3 Gb → 每文件 reads 数按读长换算；环 = 完整；线性但 18S/5.8S/28S 全覆盖 = 完整；28S 只覆盖 80 % = 不完整；组蛋白不要求 H1；标签序列名到基因名的识别（含基因簇记录不计） |
| `test_compare_and_status_from_consensus_dirs` | GetOrganelle vs SPAdes-MitoM 比较与进度 | 相同 → identical；一条包含另一条 → contained；1 个错配 → differ；只有一方 → only SPAdes-MitoM；`--status` 显示 "compared +spades" |

MitoM 自身的测试（`~/mybiosofts/Mitom/tests/`，mitotools 环境）：**485 项通过，3 项跳过**（2026-10-28，提交 6f75f73）。

---

## 3. 功能测试（端到端）

### 3.1 完整流程（测试 E：HNP1sub，`wgs_skim run`，全部步骤）

- 命令：`wgs_skim run -r smk_test/00.rawdata -o run_release -j 64`
- 结果：**21/21 个任务全部成功，0 个错误**，28 分钟（1697 s），内存峰值 49.7 GB（来自 BLAST 读取本地库，其他步骤 < 8 GB）。
- 生成：`summary/samples.tsv`、`summary/report.html`；线粒体 15 662 bp（线性，GetOrganelle round 2，"no change after round2"），13 PCG / 22 tRNA / 2 rRNA；topology `linear_unjoined`（末端诊断：37/68 条 reads 继续延伸进入重复区）。
- 覆盖步骤：fastp → MEANGS → GetOrganelle 线粒体 → 迭代 → SPAdes → BUSCO → mt_select → MitoM 注释 → topology → rDNA（GetOrganelle + MitoM）→ 组蛋白（GetOrganelle + MitoM）→ SPAdes 中的核基因（mitom nuc）→ 拷贝数 → GetOrganelle vs SPAdes 比较 → jellyfish + GenomeScope → BLAST → 汇总报告。

### 3.2 核基因：rDNA（18S-ITS1-5.8S-ITS2-28S）与组蛋白（H1–H4）

| 测试 | 内容 | 结果 |
|---|---|---|
| A（0.05 / 0.1 Gb） | 故意给极小数据量，看第 2 轮是否触发 | 第 2 轮按设计触发；数据太少时组蛋白无结果，比较表正确显示 "only SPAdes" |
| A2（默认 5 / 10 Gb） | rDNA 迭代 | 第 1 轮未闭合 → 第 2 轮用第 1 轮结果作种子 → 被选用；rDNA 5 个基因与 SPAdes-MitoM **完全相同**（18S 1882、ITS1 377、5.8S 158、ITS2 226、28S 3773 nt） |
| A3–A4 | 组蛋白迭代 | 第 2 轮得到**闭合的组蛋白基因簇单元**；H2B、H3、H4 与 SPAdes **完全相同**，H2A 1 个碱基差异（簇内不同拷贝），H1 只在 SPAdes 中（单拷贝，不在基因簇里） |
| B | `wgs_skim nuclear` 用于 wgs_skim 结果目录 | 与流程内结果一致；复用流程已有的 SPAdes-MitoM 结果 |
| C | `wgs_skim nuclear` 用于旧流程（Temp_scripts）结果目录 | 与 B 结果相同；已有的 rDNA 组装被当作第 1 轮，未闭合时自动进行第 2 轮 |
| 断点续跑 | 对 C 再运行一次 | 0.97 秒完成（全部复用） |
| D | 样品表 `scaffolds` 列（使用已有 SPAdes 组装）+ `--steps nuclear` + `--mem-gb 60` | SPAdes 未重新运行（链接已有组装）；核基因比较结果与 A4 相同；447 秒 |

### 3.3 `wgs_skim compare`（GetOrganelle vs SPAdes-MitoFinder 线粒体）

真实数据（CHY006 / CHY007 / CHY008，同一物种）：

| 样品 | 结论 | 说明 |
|---|---|---|
| CHY006 | identical | 完全相同 |
| CHY007 | consistent | GetOrganelle 线性 16 318 bp，MitoFinder 环化后 17 228 bp：多出的 910 bp 在控制区 |
| CHY008 | consistent | 两处 1 bp 差异，均在 poly-T（14 与 11）中 |

合成错误（CHY006 参考上人为制造，测 QUAST 式结构检查）：

| 制造的错误 | 检测结果 |
|---|---|
| 环旋转 / 反向互补 | 无事件（正确） |
| COX1 中 1 bp 缺失 | COX1 frameshift |
| ND4 区段倒位 | 2 个 inversion 断点，ND4 标为反向 |
| ATP8 删除 182 bp | ATP8 missing，182 bp deletion |
| 2 kb 区段移位 | 3 个 relocation（+2000、−6000、+4000） |
| 切成 3 段线性片段 | COX1 split、rrnL partial、两段参考缺口 |

### 3.4 资源与调度

| 测试 | 结果 |
|---|---|
| 干运行（dry run）`--mem-gb 100` | 命令中出现 `--resources mem_mb=102400`；每个任务显示自己的 mem_mb |
| `--steps summary` 单独运行 | 只汇总结果目录中已有的内容（12 个任务，不重跑组装） |
| `--steps getorganelle` | 运行到最终线粒体（getorganelle_mt → mt_iterate → mt_select） |
| 测试 F：中间文件清理 | SPAdes 目录 2.1 GB → 817 MB（删除 K* / tmp）；每轮 GetOrganelle 23 MB → 112 KB（删除 extended reads 等） |
| 实测各步骤内存峰值（0.6 Gb 样品） | fastp 0.4 GB、MEANGS 7.6、GetOrganelle 线粒体 3.2、rDNA 0.8、组蛋白 1.6、SPAdes 6.9、BUSCO 2.0、jellyfish 6.8、MitoM 0.2、topology 0.9 GB |

### 3.5 真实批量（MJ20260820041，22 个样品）

- rDNA：22 个样品全部注释，22 个 `.rdna.gbf`；**NCBI table2asn 校验 22/22 个 0 错误 0 警告**。18S 在 CHY / CWSC 样品中均为 1817 nt。
- 组蛋白与 SPAdes 比较：批量运行中（`results/11.nuclear_markers/`，`wgs_skim nuclear --status -o ...` 查看进度）。
- 第 2 轮（10 Gb）的实际效果：CHY006 组蛋白第 2 轮改善（标签覆盖 1292 → 2794 bp，并闭合）；CHY008 rDNA 第 2 轮没有任何变化（5721 = 5721 bp）→ 据此改为：线性但基因完整的单元不再做第 2 轮。

---

## 4. 测试中发现并已修复的错误

| # | 发现于 | 错误 | 修复 |
|---|---|---|---|
| 1 | 测试 A2 | 组蛋白在另一科样品（*Pisione*）上 "Too few seed reads"，无结果 | 组蛋白用 bowtie2 `--very-sensitive-local` 招募种子 reads：H2A–H4 完整组装在 9.7 kb 上 |
| 2 | 测试 A3 | 第 2 轮"闭合"成一个 155 bp 的简单重复环，取代了真正的基因簇 | 每轮按目标基因覆盖碱基数打分；闭合的结果必须保留 ≥ 90 % 的基因覆盖才采用 |
| 3 | 测试 A3 | GetOrganelle 失败时解压的 reads 不删除（每个失败样品数 GB） | 每轮结束统一清理（`collect_fasta.clean`） |
| 4 | 真实数据 CHY007 | MitoM：contig 末端 20 nt 的 ITS1 被误注为"5′ 不完整 18S"，ITS1 被截短 | MitoM：无比对支持的片段需与参考基因末端相似 ≥ 70 % 才保留（MitoM PR #2） |
| 5 | 真实数据 6 个 SPAdes 组装 | MitoM：3 个样品的 H3 共识实为单拷贝 H3.3 变体，H2A 实为 H2A.X | MitoM：组蛋白拷贝分为基因簇 / 替换变体 / 单拷贝，共识优先取基因簇拷贝（MitoM PR #1） |
| 6 | 代码审查 | 报告在拷贝数为 "<1 (fragments)" 时崩溃，整个汇总失败 | 安全解析 |
| 7 | 代码审查 | BLAST 元数据下载整条序列（染色体级记录超时）→ 物种名静默为空 | 只下载注释头；失败时表中写 "(efetch failed)" |
| 8 | 代码审查 | 两个样品名互为前缀（X 与 X.1）时 BLAST 用错样品的序列 | 按样品目录对应 |
| 9 | 代码审查 | auto 模式下 SPAdes 闭合的"环"不经 reads 验证 | SPAdes 来源的环也做 reads 验证 |
| 10 | 代码审查 | 弱相关的两条序列被判为 "identical"（端部空位不计分） | 较短序列须被覆盖 ≥ 90 % 才判相同，否则 "no overlap" / "differ" |
| 11 | 代码审查 | GetOrganelle 中途被杀仍写下 `done`，以后永不重跑 | 每轮都失败时不写 `done`，下次自动重试 |
| 12 | 代码审查 | 组装更新后旧注释被复用 | 注释标记中加入输入文件大小与时间 |
| 13 | 代码审查 | MitoFinder scaffold 含缺口时环化崩溃 | 跳过图遍历并注明 |
| 14 | 代码审查 | 个人配置 `~/.wgs_skim/config.yaml` 中的 `db_dir` 被命令行默认值覆盖 | 用户写的值优先 |
| 15 | 代码审查 | Excel 导出的 UTF-8 样品表（BOM）报 KeyError | 按 utf-8-sig 读取 |
| 16 | 运行中 | 整个机器内存紧张时后台任务被停止，一夜无结果 | 每个规则声明内存，Snakemake 按核数和内存调度；核基因批量每完成一个样品就写 `summary.tsv`；`--status` 查看进度 |

---

## 5. 尚未测试（如实列出）

| 部分 | 状态 |
|---|---|
| 其他类群预设（mollusc / fish / crustacean） | 未用真实数据测试，只有环节动物 |
| 全新机器上从零安装（`wgs_skim install` 含数据库下载） | 未在干净机器上测试 |
| 在线 NCBI BLAST（`-remote`） | 测试均用本地小库（`build-blastdb`） |
| 集群（Slurm 等 executor） | 未测试 |
| 读深度图（旧脚本 using_bwa_draw_coverage） | 未移植（MitoM `mito_plot.py` 可画） |
| `wgs_skim compare` 在最新代码上的真实数据重跑 | 单元测试覆盖；真实数据最后一次在 2026-10-28 |
| Nextflow 版本 | 开发中（见 README） |

---

## 6. 如何重复测试

```bash
# 1. 单元测试（秒级，不需要外部软件；比较测试需要 Biopython）
python -m pytest -q tests/

# 2. 做一个双样品小数据集（每个 1 M 对 reads）
mkdir small
for s in A B; do for r in 1 2; do zcat $s.R$r.fq.gz | head -n 4000000 | gzip > small/${s}_R$r.fq.gz; done; done

# 3. 端到端：完整流程 + 逐样品检查每个输出 + wgs_skim nuclear + 单元测试；退出码 = 失败项数
tests/e2e.sh small/ e2e_out/
```

---

## 7. 双样品快速回归（data2，`tests/e2e.sh data2 e2e2`）

- 两个样品（`pisione`、`cwsc010`，不同物种，各 1 M 对 reads = 0.3 Gb），完整流程 **39/39 个任务，0 个错误，26 分钟**。
- 逐样品检查 11 项输出（fastp、线粒体、MitoM 注释、rDNA、rDNA GenBank、组蛋白 GenBank、SPAdes、BUSCO、基因组大小、
  核基因拷贝数、GetOrganelle vs SPAdes）+ BLAST、汇总表（2 行）、HTML 报告、`wgs_skim nuclear`、单元测试：**0 项失败**。
- 线粒体在 0.3 Gb 下未闭合（cwsc010 8 954 bp、pisione 15 776 bp，topology `linear_low_reads`）：数据量太小，属预期；
  这个数据集用于检查流程是否完整跑通，不用于评价组装质量。

## 8. 批量管理：Snakemake

- `wgs_skim run`、`wgs_skim nuclear --results`、`wgs_skim compare --results` 全部由 Snakemake 调度（核数 + 内存、每个任务
  一个日志、断点续跑）。`nuclear`：小样品 8 个任务成功，再次运行 "Nothing to be done"（2.2 秒）；`compare`：CHY006/007/008
  identical / consistent / consistent，CHY002 无 MitoFinder 结果被列为 skipped，110 秒。
- SPAdes 各模式与 Trinity（`--spades-modes`、`--steps rna`）：干运行检查依赖图（rna = rnaSPAdes + Trinity + 4 个 BUSCO；
  genome + meta + plasmid 各自组装 + BUSCO；未知模式报错并列出可用模式）；真实 RNA 双样品测试见第 9 节。

## 9. 命令矩阵与转录组

### 9.1 命令矩阵（`tests/cli_matrix.sh data2 e2e2`，33 项，全部通过，约 1 分钟）

- 子命令：`--version`、`samples`（有 / 无参数）、`config`、`run` 无输入时的提示、每个子命令的 `-h`。
- 每个 `--steps` 值的干运行，检查 Snakemake 将运行的规则：everything、fastp、meangs、mito、getorganelle、annotate、topology、
  nr、spades、spades,busco、`--spades-modes genome,meta,plasmid,isolate`、rna、trinity、genome_size、nuclear、blast、summary 单独。
- 错误输入：未知步骤、未知 SPAdes 模式（都给出可用值）；`--mem-gb 50` → `--resources mem_mb=51200`。
- 对已完成结果：`circularize`（无 reads）、`compare` 单样品、`nuclear` 复用已有结果（Nothing to be done）、`nuclear --status`、
  `nuclear -n`、单元测试。

### 9.2 转录组（`wgs_skim run -r rna2 -o rna2_out --steps rna`）

两个 RNA-seq 样品（A1、R1，各约 50 MB × 2，约 0.1 Gb）。rnaSPAdes（`--rna -k 21,33,55`）+ Trinity 2.15.2 + 各自 BUSCO（metazoa_odb12，`-m transcriptome`）。

| 样品 | rnaSPAdes 转录本 | 时间 / 峰值内存 | BUSCO |
|---|---|---|---|
| A1 | 20 573 | 87 s / 1.3 GB | C:12.1% [S:11.9%, D:0.1%], F:22.8%, M:65.2%, n:672 |
| R1 | 27 620 | 92 s / 1.2 GB | C:13.2% [S:12.9%, D:0.3%], F:29.8%, M:57.0%, n:672 |

BUSCO 低是数据量太小（完整转录组通常需要 ≥ 2 Gb）所致；本测试检查流程本身。

Trinity：第二阶段为每个 read 簇单独运行一次（A1 有 13 961 个簇），在机器被其他 SPAdes 任务占满时很慢（20 分钟约 7 %）。
结果（转录本数、BUSCO）运行结束后补在这里。

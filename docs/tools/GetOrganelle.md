# GetOrganelle 排错参考

## 1. 简介与用途

GetOrganelle 从基因组浅测序 reads 中按种子招募目标 reads（bowtie2 种子 → 按 word 迭代延伸 `-R` 轮），用 SPAdes 组装
招募到的 reads，再用标签库（`--genes`，BLAST）在组装图中筛出目标并“解缠”（disentangle）成环或线性路径。种子只用来
开始招募，组装是 de novo 的（官方 FAQ）。

**在 wgs_skim 中的作用**（三处，同一个 conda 环境 `envs/getorganelle.yaml`）：

| 目标 | 规则 / 脚本 | 模式 |
|---|---|---|
| 线粒体第 0 轮 | `mito.smk` rule `getorganelle_mt` | `-F animal_mt`，MEANGS 种子（无则内置） |
| 线粒体迭代（不成环时） | rule `mt_iterate` → `scripts/mt_iterate.py` | `-F animal_mt`，上一轮结果作种子，全部 reads（≤ 20 Gb） |
| 核 rDNA（18S-ITS-28S）、组蛋白基因簇 | `nuclear.smk` → `scripts/nuclear_markers.py assemble` | `-F anonym` + 自带 seed/label（`resources/rdna`、`resources/histone`） |

- GitHub：https://github.com/Kinggerm/GetOrganelle；Wiki（Usage / FAQ / Examples）：https://github.com/Kinggerm/GetOrganelle/wiki
- 论文：Jin J-J, Yu W-B, Yang J-B, et al. (2020) GetOrganelle: a fast and versatile toolkit for accurate de novo
  assembly of organelle genomes. *Genome Biology* 21: 241. doi:10.1186/s13059-020-02154-5
- 数据库：`get_organelle_config.py -a animal_mt`（rule `getorganelle_db`，`~/.wgs_skim/db/getorganelle_db`，animal_mt 0.0.1：
  种子 1853 条、标签 7373 条）

## 2. 版本

| | 版本 | 说明 |
|---|---|---|
| wgs_skim 固定 | **1.7.7.1**（bioconda，2024-04-03）+ SPAdes 3.15.5、bowtie2 2.4.5、**blast 2.15.0**、python 3.10 | 最新正式版 |
| GitHub master | **1.8.0.1**（未发正式 release，未上 bioconda；1.8.0.0 开发于 2023-08，2024-10-29 合并 intermediate_graph 分支 PR #281，2025-02-13 兼容性修复） | 只能从源码安装 |

1.8.0.x（`GetOrganelleLib/versions.py`）与我们相关的改动：读入文件改为**软链接而不是复制**、不同轮次的临时文件分开、
去掉 scipy/sympy（改用 gekko）、覆盖度过滤去掉上限、`--continue` 时索引检查加速；intermediate_graph 分支据作者说
修复了 anonym 模式线性结果偏短（标签只分给一条 contig，issue #239）。

**是否升级：暂不升级，值得单独测试 rDNA。** 理由：未正式发布、无 bioconda 包；已有用户报告 1.8.0.1 的 disentangle
断言错误（#390、#388）。但 #239 的修复正对应我们 18/22 个 rDNA 单元是线性的情况（TESTING 3.5），可在一个独立环境里用
3–5 个 MJ 样品对比 rDNA 长度与基因完整度后再决定。

## 3. wgs_skim 中的调用方式

线粒体第 0 轮（16 线程，mem 16 GB）：

```bash
get_organelle_from_reads.py -1 R1 -2 R2 -o mito/<s>/round0 -F animal_mt -R 10 -k 21,45,65,85,105 -t 16 \
    --config-dir <db>/getorganelle_db [-s mito/<s>/meangs/<s>.seed.fasta] --overwrite
```

线粒体迭代第 k 轮（`mt_iterate.py`；第 0 轮不成环才运行，最多 `mt_iter_rounds` 轮）：

```bash
get_organelle_from_reads.py ... -F animal_mt -R 20 -k 21,45,65,85,105 --overwrite \
    --max-reads <cap 或 inf> --reduce-reads-for-coverage inf -s round<k>.seed.fasta   # 种子 = 上一轮结果
```

停止条件：某轮成环；或连续两轮序列完全相同（"no change after roundN"）。cap = `mt_iter_max_data_gb` Gb ÷ 2 ÷ 平均读长
（来自 fastp JSON）。

rDNA / 组蛋白（`nuclear_markers.py assemble`，每轮）：

```bash
get_organelle_from_reads.py -1 R1 -2 R2 -o <dir> -F anonym -s <seed 或上一轮结果> --genes <label> -P 0 -t 16 \
    -R 10 -k 21,45,65,85,105 --max-reads <5 Gb 的 reads 数> --reduce-reads-for-coverage inf --overwrite \
    --bowtie2-options "--very-fast -t"            # 组蛋白: "--very-sensitive-local -t"
```

第 1 轮 5 Gb；单元不完整（不成环，且核心基因未被标签 ≥ 95 % 覆盖）才做第 2 轮 10 Gb（`nuclear_max_gb: "5,10"`）。

可调配置键：`getorganelle_kmers`、`getorganelle_rounds`（-R 第 0 轮）、`mt_iter_rounds`、`mt_iter_R`、`mt_iter_max_data_gb`
（0 = 不限）、`nr_seed` / `nr_label` / `histone_seed` / `histone_label`、`nuclear_max_gb`、`meangs`、`db_dir`。
每轮结束 `collect_fasta.clean()` 删除 `extended_*.fq`、`extended_spades/`、`seed/`、`[0-9]-*.fastq`、`*temp.indices*`。

## 4. 重要参数说明

| 参数 | 1.7.7.1 默认 | 说明 |
|---|---|---|
| `-F` | — | animal_mt / anonym（anonym 必须给 `-s` 和 `--genes`） |
| `-s` | `SeedDatabase/<F>.fasta` | 只用于招募初始 reads；动物线粒体需要近缘种子，没有时用上一次的结果作种子（FAQ） |
| `--genes` | `LabelDatabase/<F>.fasta` | 图筛选 / 解缠用；自建标签库可能反而不成环（#371），线粒体保持默认标签 |
| `-R` | 10（animal_mt），-P 0 时 inf | 延伸轮数；深度低、种子远时多给 |
| `-k` | 21,55,85,115 | SPAdes k；结果不完整时用更宽更密的梯度（FAQ）。150 bp reads 可加 115/127（ISSUES #19 的待试项） |
| `-w` | 自动估计 | 估计过小会耗内存（#214）；出环就不用调 |
| `--max-reads` | 3E8（animal_mt）、7.5E7（anonym）；> 536870911 视为 inf | 每个文件最多用多少 reads（硬上限） |
| `--reduce-reads-for-coverage` | 500 | 估计目标覆盖度超过此值就减少 reads；`inf` 关闭。日志里的覆盖度是**子集**估计（ISSUES #24） |
| `-P` | 2E5 | 预分组高覆盖 reads；目标拷贝少或位点短（< 3 kb）时设 0（FAQ 自定义位点） |
| `--max-extending-len` | 自动 | FAQ 建议自定义位点按需设（如 100）以免延伸进旁侧区 |
| `--expected-max-size` | 25000（animal_mt） | 解缠时的最大目标长度 |
| `--memory-save` | — | = `--out-per-round -P 0 --remove-duplicates 0`，MemoryError 时用 |
| `--keep-temp` / `--continue` | 关 | 不加 keep-temp 时临时文件只在**正常结束**时删除 |
| `--bowtie2-options` | `--very-fast -t` | 种子招募的 bowtie2 参数 |

## 5. 已知问题与 issues

| 现象 / 报错 | 原因 | 解决 | 链接 |
|---|---|---|---|
| `Too few seed reads found!` | 线粒体 reads 太少，或种子太远 | 近缘种子 `-s`；加数据；`--max-reads` 无帮助（并没有在减 reads） | [#51](https://github.com/Kinggerm/GetOrganelle/issues/51) |
| 自建种子 / 标签不成环，默认 animal_mt 能成环 | 标签库影响解缠 | 线粒体用默认标签；种子可换 | [#371](https://github.com/Kinggerm/GetOrganelle/issues/371) |
| 用类群专用种子 / 标签反而更短 | 同上，作者未答复 | 对比默认库 | [#377](https://github.com/Kinggerm/GetOrganelle/issues/377) |
| 动物线粒体重复区不成环 | 短 reads 跨不过重复 | 加大 k 只能解决部分；以 N 代替或长读长 | [#55](https://github.com/Kinggerm/GetOrganelle/issues/55) |
| anonym（ITS / rDNA）线性结果偏短 | 1.7.x 标签只分给一条 contig | 手工检查 `filtered*.fastg`；intermediate_graph 分支（1.8.0 master）已改 | [#239](https://github.com/Kinggerm/GetOrganelle/issues/239) |
| `-F anonym` 不给 `--target-genome-size` 时 IndexError | 1.7.7.1 参数解析 | 加 `--target-genome-size`（本项目未遇到） | [#337](https://github.com/Kinggerm/GetOrganelle/issues/337) |
| `zsh: killed` / 进程被杀 / MemoryError | 内存；-w 估得太小 | 加大 `-w`、`--memory-save`、减 reads | [#214](https://github.com/Kinggerm/GetOrganelle/issues/214)、[#212](https://github.com/Kinggerm/GetOrganelle/issues/212)、[#347](https://github.com/Kinggerm/GetOrganelle/issues/347) |
| 只有 `extended_*.fq`，没有组装 | SPAdes 前后失败；看 get_org.log.txt | 附完整日志排查 | [#341](https://github.com/Kinggerm/GetOrganelle/issues/341) |
| `Disentangling unsuccessful` / `Incomplete/Complicated graph` | 图不完整；1.7.7.1 起由 "failed" 改为 "unsuccessful" | 输出 scaffolds 路径；加 reads / 轮数 / k | [#308](https://github.com/Kinggerm/GetOrganelle/issues/308)、[#322](https://github.com/Kinggerm/GetOrganelle/issues/322) |
| `ERROR: Blast is not available!` | blastn 动态库缺失（libgomp 等） | 终端直接运行 `blastn` 看报错 | [#208](https://github.com/Kinggerm/GetOrganelle/issues/208) |
| conda 安装失败 / glibc 不兼容 / scipy 导入错 | 求解器 / 旧依赖 | 新建独立环境；1.7.7.1 已修 scipy 导入 | [#216](https://github.com/Kinggerm/GetOrganelle/issues/216)、[#132](https://github.com/Kinggerm/GetOrganelle/issues/132) |
| 解压 .gz 很慢 | 开始时把 reads 解压 / 复制进输出目录 | 1.8.0 改为软链接（未发布） | [#53](https://github.com/Kinggerm/GetOrganelle/issues/53) |
| 没有 animal_nr 库 | 作者只做了植物 / 真菌 nr | 用 `-F anonym` + 自建 seed/label | [#136](https://github.com/Kinggerm/GetOrganelle/issues/136) |

## 6. 本项目中遇到的问题

1. **blast 未固定 → 一条线粒体都没有（ISSUES #45，已修）**：yaml 只钉 python/getorganelle/spades 时求解器选了 blast 2.5.0
   （2018），各轮都是 `No sequence hit our LabelDatabase!`。现钉 blast 2.15.0、bowtie2 2.4.5（与可用的旧环境相同）。
2. **迭代只换种子是白跑（ISSUES #14、#19、#47）**：第一版迭代沿用 -R 10 和相同 reads，三轮结果完全一样；改为全部 reads +
   -R 20 后 4 个不成环样品（CHY007、HN_Pisione_1/2、WenSC001）也只变长 0–74 bp，第 1 轮后即收敛。HNP1sub round1→3
   md5 相同 → 加“与上一轮相同即停止”。原因是生物学的：末端停在 (AT)n 微卫星和 poly-G（软剪切 reads 诊断），短 reads
   跨不过；改 `-w 45/55/65/75` 招募量不变（817–821 KB）。此时改用 SPAdes 兜底（`mt_select`）和 topology 诊断。
3. **磁盘：解压 reads 与索引（ISSUES #25、#51；TROUBLESHOOTING）**：GetOrganelle 把 .gz 解压成 `<n>-<文件名>.fastq` 并建
   `temp.indices.*`，只在正常结束且无 `--keep-temp` 时删。中途被杀的 CHY007 残留 233.7 GB；全数据（1.6 亿对）一轮占
   221 GB，把系统盘写满（`OSError: [Errno 28] No space left on device`，在 make_read_index）。对策：输出放大盘，预留约
   1.5× 未压缩 reads；`mt_iter_max_data_gb` / `nuclear_max_gb` 限量；每轮 `clean()`；被杀后确认进程已退出再删这三类文件。
4. **失败后仍写 done（TESTING #11，已修）**：被杀的 GetOrganelle 也 `touch .done`，以后永不重跑 → 每轮都失败时不写。
5. **组蛋白 "Too few seed reads"（TESTING #1，已修）**：另一科（*Pisione*，HNP1sub）的种子 reads 只有 391 bytes；改
   `--very-sensitive-local` 后 153 kb，H2A–H4 在 9.7 kb 上完整。只含 H3 的标签只留下 H3 片段（2.2 kb）→ 标签用 H1–H4 簇。
6. **第 2 轮“成环”成 155 bp 简单重复（TESTING #2，已修）**：按标签覆盖碱基数打分，成环结果需保留 ≥ 90 % 基因覆盖才采用。
7. **线性 rDNA 不必做第 2 轮**：22 个 rDNA 中 18 个线性但基因完整；CHY008 第 2 轮 5721 = 5721 bp 无变化 → 核心基因被标签
   覆盖 ≥ 95 % 即视为完整。
8. **多路径结果**：一个样品有 graph1.1、1.2…（只差 flip-flop 方向，FAQ），`collect_fasta.best()` 只取第一个 complete 文件，
   没有 complete 才用 scaffolds（ISSUES #5）。
9. **从 SPAdes 全图提取不可行（ISSUES #49）**：`get_organelle_from_assembly.py -g`（4.95 GB fastg）“Slimming” 5 h 50 min 无进展。
10. **线程利用率低（ISSUES #10）**：延伸阶段基本单线程，6 样品 × 16 线程整机负载约 62/128 → 多样品少线程更划算（未改）。
11. 资源实测（HNP1sub 0.6 Gb）：第 0 轮 254–276 s / 3.2–6.1 GB；迭代 245–266 s / 4.7–5.2 GB；rDNA 381–474 s / 0.5–0.8 GB；
    组蛋白 395–483 s / 1.5–1.6 GB。规则内存 = max(8, 1.2 × 实际用到的 Gb) GB（`common.smk`）。

## 7. 排错清单

- 日志：`logs/getorganelle_mt/<s>.log`、`logs/mt_iterate/<s>.log`（每轮命令、rc、序列数、circular、删除的 MB、stop 原因）、
  `logs/getorganelle_nr|histone/<s>.log`；每个输出目录里的 `get_org.log.txt` 保留。状态：`mito/<s>/<s>.getorganelle.json`
  （`circular`、`round`、`stop`、`max_reads_per_file`）。
- 关键字（get_org.log.txt）：
  - `Seed reads made: ... (N bytes)` —— 种子招募量；几十 KB 以下基本组装不全。
  - `Too few seed reads found` —— 种子太远 / 数据太少。
  - `No sequence hit our LabelDatabase!` —— 标签不匹配，或 **blast 版本错**（先查 `blastn -version`）。
  - `Hit the round limit N and terminated` —— -R 用完（不一定是错）。
  - `Disentangling unsuccessful: 'Incomplete/Complicated graph: please check around EDGE_x'` → 输出 `*scaffolds*.path_sequence.fasta`
    = 不完整；`complete.graph1.1.path_sequence.fasta` 且标题含 `(circular)` = 成环。
  - `Reads used` / `Estimated ... base-coverage` —— 实际用量（子集估计，不等于测序深度）。
  - `MemoryError` / `Killed` / `No space left on device` / `Errno 28`。
  - `Blast is not available` / `bowtie2-align died`。
- 先看磁盘：`du -sh mito/*/round* rdna/*/getorganelle`，残留 `[0-9]-*.fastq`、`*temp.indices*` 说明上次异常退出。
- 复现单轮：从 `logs/mt_iterate/<s>.log` 复制 `roundk:` 后的命令，在 getorganelle 环境里换一个 `-o` 运行。

## 8. 资料来源与查阅日期（2026-10-29）

- GitHub Releases 1.7.7.0 / 1.7.7.1；master `GetOrganelleLib/versions.py`（1.8.0.0、1.8.0.1）；commits 列表
- Wiki FAQ：https://github.com/Kinggerm/GetOrganelle/wiki/FAQ （de novo 说明、自定义位点、MemoryError、不完整结果、nr 是否成环）
- Issues #51 #53 #55 #132 #136 #208 #212 #214 #216 #239 #308 #322 #337 #341 #347 #371 #377 #388 #390（gh api）
- 本地 `get_organelle_from_reads.py --help`（1.7.7.1，默认值以此为准）；anaconda.org bioconda/getorganelle 最新 1.7.7.1
- 本项目：`workflow/rules/mito.smk`、`nuclear.smk`、`scripts/mt_iterate.py`、`nuclear_markers.py`、`collect_fasta.py`、
  `envs/getorganelle.yaml`、`config/config.yaml`、`docs/TESTING.md` §3–4、README Troubleshooting；MJ20260820041 `ISSUES.md`
  #5 #10 #14 #19 #24 #25 #45 #47 #49 #51、`Temp_scripts_dev/wgs_skim_pipeline/TROUBLESHOOTING.md`；`wgs_skim_test/*/benchmarks/`

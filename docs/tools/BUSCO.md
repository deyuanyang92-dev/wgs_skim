# BUSCO（wgs_skim 排错参考）

## 1. 简介与用途

BUSCO 用单拷贝直系同源基因集（OrthoDB lineage）评估组装完整度。在 wgs_skim 中：

- `rule busco`（`workflow/rules/assembly.smk`）：对 SPAdes 基因组组装 `assembly/<s>/<s>.scaffolds.fasta` 做 `-m genome`；
  其 `busco/run_*/full_table.tsv` 还被 `nuc_copies.py` 用来取单拷贝基因的覆盖度，作为 rDNA / 组蛋白拷贝数的分母。
- `rule busco_assembly`（`workflow/rules/transcriptome.smk`）：其他 SPAdes 模式和 Trinity 的组装；rna / rnaviral / trinity 用
  `-m transcriptome`，其余用 `-m genome`。
- `rule busco_lineage`：`busco --download <lineage>` 一次性下载到 `<db_dir>/busco_downloads/lineages/`（`busco_lineage` 为本地目录时跳过）。
- Nextflow：`nextflow/modules/assembly.nf` 的 `BUSCO_LINEAGE`（`storeDir`）与 `BUSCO`。

链接：
- 代码：https://gitlab.com/ezlab/busco （issue 在 GitLab，不在 GitHub）
- 用户手册：https://busco.ezlab.org/busco_userguide.html
- 数据集列表：https://busco-data.ezlab.org/v5/data/file_versions.tsv
- 论文：Manni et al. 2021 MBE, doi:10.1093/molbev/msab199；odb12 / BUSCO 6：Tegenfeldt et al. 2025 NAR 53(D1):D516
  （手册引用此文；DOI 应为 10.1093/nar/gkae987，未在手册页面上核对，标记为待确认）

## 2. 版本

| | 版本 | 日期 |
|---|---|---|
| wgs_skim 固定 | `busco =6.0.0`（`workflow/envs/busco.yaml`） | 2025-07-09 |
| 最新 | 6.1.0 | 2026-06-05 |

6.0.0 主要变化（相对 5.x）：miniprot 后处理改进（#728 #739 #742 #744 #772 #795 #809）、`--plot`、接受 gzip 输入、
可分别下载 odb10 / odb12、作业失败写入 busco.log（#755）、解压到可写位置（#693 #779 #789）。

6.1.0：支持 **odb12.2** 数据集（重新引入 `lengths_cutoff`：Complete 的长度下限改为 median / sMAD 规则，上限 1.5×median），
仍兼容 odb12 / odb10。手册原文："the odb12.2 datasets can only be correctly parsed by v6.1.0 or later"。

是否升级：**暂不必须，但应留意**。6.0.0 + odb12 能正常工作（e2e 实测）。数据服务器上现在同时有 `metazoa_odb12`
（2026-05-22 重打包）与 `metazoa_odb12.2`（2026-05-13）；只要 `busco_lineage` 写的是不带 `.2` 的名字，6.0.0 就没问题。
若以后改用 odb12.2，必须同时升到 6.1.0。另：6.1.0 有一个 open issue #855（odb12.2 下 Fragmented 判定与手册描述不一致），
#852（6.1.0 不产生任何输出），升级前先在小数据上回归。

## 3. wgs_skim 中的调用

```bash
# 下载（rule busco_lineage / BUSCO_LINEAGE）
busco --download metazoa_odb12 --download_path <db_dir>/busco_downloads
# 评估（rule busco；busco_assembly 同理，-o busco_<asm>，-m 由组装类型决定）
rm -rf <out>/busco
busco -i <s>.scaffolds.fasta -l <lineage 目录的完整路径> -m <busco_mode，默认 genome> -o busco --out_path assembly/<s> \
      -c 16 --offline -f
cp assembly/<s>/busco/short_summary.*.txt assembly/<s>/<s>.busco.txt
```

- `-l` 给**完整路径**：手册说明给完整路径时 BUSCO 不再自动管理数据集（不会联网、不在 `./busco_downloads` 里找）。
  Nextflow 中暂存的链接名是 bare name，必须 `$(readlink -f lineage)`（见第 6 节 #17）。
- `--offline`：不下载任何文件（也不检查数据集新版本）。
- `-f`：覆盖已有结果目录。规则里先 `rm -rf` 再 `-f`，两者重复但无害；**不要**改成 `-r`（restart），因为 Snakemake 失败重跑时
  旧目录可能是不完整的。
- 资源：`threads: 16`、`mem_mb: 16000`（Nextflow `16000 × attempt`）。
- 可配置键（`config/config.yaml`）：`busco_lineage`（名称或本地目录；空 = 分类群预设，annelid → metazoa_odb12，mollusc →
  mollusca_odb12，fish → actinopterygii_odb12，crustacean → arthropoda_odb12）、`db_dir`、`assembly`、`steps: busco`、
  `spades_modes`。CLI：`--busco-lineage`。
- `busco_mode`（genome / transcriptome / proteins）：Snakemake `rule busco` 用作 `-m {params.mode}`（工作区 2026-10-29
  的未提交修改中刚加入；`wgs_skim busco -m` 子命令设置它）。**Nextflow 的 `BUSCO` 进程仍写死 `-m genome`**
  （`nextflow/modules/assembly.nf:103`），两个引擎在 `busco_mode` ≠ genome 时结果不同。`busco_assembly` 不看此键，用 `BUSCO_MODE` 字典。

## 4. 重要参数

| 参数 | 含义 | wgs_skim 取值 / 建议 |
|---|---|---|
| `-m genome/transcriptome/proteins` | 模式 | genome；RNA 组装 transcriptome |
| `-l` | lineage 名称或路径 | 完整路径 |
| `--offline` | 不联网 | 总是加 |
| `--download_path` | 数据集存放 | 仅下载时用 |
| `-c` | 线程；不给则只用 1 个 CPU | 16 |
| `--miniprot`（默认，真核 genome） / `--metaeuk` / `--augustus` | 基因预测器 | 默认 miniprot（日志可见 `miniprot_index` / `miniprot_align`）；transcriptome 模式固定用 metaeuk |
| `-f` / `-r` | 覆盖 / 从中断处继续 | `-f` |
| `--skip_bbtools` | 不算组装统计 | 可加，省几秒；wgs_skim 不用 BUSCO 的 N50 |
| `--limit` | 每个 BUSCO 考虑的候选区域数（默认 3） | 默认 |
| `--contig_break` | 多少个 N 视为 contig 断点（默认 10） | 默认 |
| `--plot`（6.0+） | 生成图 | 未用 |

## 5. 已知问题与 issues

| 现象 / 报错 | 原因 | 解决 | 链接 |
|---|---|---|---|
| `Dataset ... does not exist` / 在 `./busco_downloads` 找 lineage | `-l` 给的是名称（或链接名），离线时 BUSCO 在当前目录下查 | `-l` 给 `readlink -f` 后的完整路径 | 手册 "Offline"；本项目 #17 |
| odb12 中没有该 lineage 时离线模式报错不清楚 | 在线时会回退 odb10，离线时不会 | 换存在的 lineage，或下载 odb10 | https://gitlab.com/ezlab/busco/-/issues/839 |
| `UnboundLocalError: ... 'autoconfig'`（--auto-lineage --offline） | 缺少 placement 文件时错误信息不明确 | 不用 auto-lineage（wgs_skim 不用） | https://gitlab.com/ezlab/busco/-/issues/826 |
| miniprot `error -9` / `Job failed with error code -6` / `[morecore] insufficient memory` | miniprot 被 OOM 杀掉 | 提高 `mem_mb`，或 `--metaeuk` | https://gitlab.com/ezlab/busco/-/issues/863 、/843 |
| miniprot 流程内存极高（脊椎动物 388 Mb 基因组 170 GB） | miniprot 对齐后处理 | 6.0.0 已改进；仍高时换 `--metaeuk` | https://gitlab.com/ezlab/busco/-/issues/809 、/739 |
| transcriptome 模式 metaeuk 报错 | metaeuk 输入 / 版本问题 | 检查 `busco.log` 中 metaeuk 命令；用 bioconda 6.0.0 | https://gitlab.com/ezlab/busco/-/issues/806 、/818 |
| 子进程（如 metaeuk）失败被当成"没找到基因"，结果全 Missing | 5.x 行为；6.0 起作业失败写入 busco.log | 结果异常低时先查 busco.log 的 ERROR | https://gitlab.com/ezlab/busco/-/issues/755 |
| 下载的 lineage tar 不完整 | S3 上文件损坏 | 重新下载并核对 file_versions.tsv 中的 md5 | https://gitlab.com/ezlab/busco/-/issues/848 、/854 、/858 |
| odb12.2 结果与 odb12 的 F 比例不同 | 新的 lengths_cutoff 规则（6.1.0） | 同一批次用同一数据集版本比较 | https://gitlab.com/ezlab/busco/-/issues/855 |
| conda 环境中 `busco` 启动失败（`#!/usr/bin/env python3` 解析到别的 python） | PATH 中其他环境的 python3 在前面 | 由 Snakemake/Nextflow 激活环境运行；手动运行时先把环境 bin 放到 PATH 前面 | 本机经验（MJ20260820041） |

## 6. 本项目中遇到的问题

- **TESTING.md #17（Nextflow）**：`-l metazoa_odb12`（Nextflow 暂存的链接名）被当作 lineage 名称 → "Dataset ... does not exist"。
  修复：`-l $(readlink -f lineage)`。
- 数据集日期：本机 `/media/cody/18T_R2/Common_databases/BUSCO_database/metazoa_odb12/dataset.cfg` 为
  `creation_date=2025-02-25`，n=672；服务器上同名数据集现为 2026-05-22。**同名不等于同一版本**，批次之间比较 BUSCO 时
  应记录 `dataset.cfg` 的 creation_date（short_summary 中也有 "dataset ... (eukaryota, 2025-02-25)"）。
- 实测（e2e2，0.3 Gb 小样品）：metazoa_odb12，miniprot，43 s，峰值 RSS 1.7 GB；转录组样品 C 仅 12–13 %（数据量太小所致）。
- 分类群预设：annelid 用 metazoa_odb12（672 个 BUSCO）。服务器上有 `lophotrochozoa_odb12`（2026-05-22，下载核对：
  **1252 个 BUSCO**，75 个物种，taxid 1206795，OrthoDB 12.1，无 `lengths_cutoff` = odb12 格式，6.0.0 应可读——未实跑），
  对环节动物 / 软体动物更近、基因数近两倍，作为 `nuc_copies` 的单拷贝分母也更稳。是否改 annelid 预设需在真实样品上比较后决定。

## 7. 排错清单

1. `logs/busco/<s>.log`（或 `logs/busco_assembly/<s>.<asm>.log`）：搜 `ERROR`、`does not exist`、`error code`、`-9`、`insufficient memory`。
2. `assembly/<s>/busco/logs/` 中各子程序日志（miniprot / metaeuk / hmmsearch）。
3. 日志开头 `Running BUSCO using lineage dataset <名> (<domain>, <日期>)`：确认用的是哪个数据集。
4. `<s>.busco.txt` 内容为 `no assembly`：SPAdes 没组装（数据量 < `spades_min_clean_gb`），不是 BUSCO 的错误。
5. `benchmarks/busco/<s>.tsv` 的 `max_rss`：接近 `mem_mb` 时提高内存。
6. 结果几乎全部 Missing：先看是否子程序失败（第 5 节 #755），再看组装是否只有细胞器 / 很小。

## 8. 资料来源（查阅日期 2026-10-29）

- GitLab releases API：https://gitlab.com/api/v4/projects/ezlab%2Fbusco/releases （6.1.0、6.0.0 说明）
- 用户手册：https://busco.ezlab.org/busco_userguide.html
- 数据集版本：https://busco-data.ezlab.org/v5/data/file_versions.tsv
- Issues：#755 #806 #809 #818 #826 #839 #843 #848 #852 #855 #863（https://gitlab.com/ezlab/busco/-/issues/<号>）
- 本项目：`workflow/rules/assembly.smk`、`transcriptome.smk`、`nextflow/modules/assembly.nf`、`docs/TESTING.md` §4 §9.2、
  `../wgs_skim_test/e2e2/logs/busco/`

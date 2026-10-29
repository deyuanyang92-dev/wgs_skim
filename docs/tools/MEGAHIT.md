# MEGAHIT（wgs_skim 排错参考）

## 1. 简介与用途

MEGAHIT 是基于简洁 de Bruijn 图（SdBG）的短读长组装器，为大型 / 复杂宏基因组设计，速度快、内存省，多 k 迭代
（默认 k 21→141）。在 wgs_skim 中**不被直接调用**：

- 只在 MitoFinder 自己组装 reads 时使用：`config.mitofinder_assembly: megahit`，或为空且未运行 SPAdes 步骤时
  （`workflow/rules/mitofinder.smk::mitofinder_source()`）。规则把 `-1 R1 -2 R2 --megahit -m <mitofinder_mem_gb>` 传给 MitoFinder，
  MitoFinder 再调用 megahit。MitoFinder 本身见 `docs/tools/MitoFinder.md`。
- 注：`mitofinder.smk`、`envs/mitofinder.yaml`（`mitofinder=1.4.1`，或 `mitofinder_sif` 的 1.4.2 镜像）在 2026-10-29 仍为未提交的工作区文件。
- 本机另有独立环境 `~/miniforge3/envs/megahit`（v1.2.9），以及 MitoZ / FastMitoAssembler / captus 环境中各自的 megahit。

链接：
- GitHub：https://github.com/voutcn/megahit （README 即手册；`megahit -h`）
- 论文：Li et al. 2015 *Bioinformatics* 31:1674，doi:10.1093/bioinformatics/btv033；Li et al. 2016 *Methods* 102:3，doi:10.1016/j.ymeth.2016.02.020

## 2. 版本

| | 版本 | 日期 |
|---|---|---|
| wgs_skim 固定 | 无直接 pin；由 MitoFinder 依赖带入（本机解析到 `megahit-1.2.9`，bioconda 多个 build） | — |
| 最新 | **1.2.9**（仍为最新 release） | 2019-10-15 |

1.2.9 修复：零长度序列段错误（#244）、旧 macOS 内存检测（#232）、unitig 图刷新的错误断言（#237）；新增 `--verbose`。
2019 年后无新 release，仓库只有构建相关的合并（2025：aarch64 支持、CMake 最低版本修复）。**无可升级版本；不必操心升级。**

## 3. wgs_skim 中的调用方式（经 MitoFinder）

wgs_skim → MitoFinder：

```bash
mitofinder -j <s> -1 R1 -2 R2 --megahit -r ref.gb -o <genetic_code> -p {threads} -m <mem_GB> -t <trna> --adjust-direction
```

MitoFinder 1.4.2 `runMegahit.py`（本机 `~/mybiosofts/MitoFinder`，1.4.1 是否相同未核对）实际执行：

```bash
megahit -1 R1 -2 R2 -o <s>_megahit --out-prefix <s>_megahit --min-contig-len <MitoFinder -l, 默认 100> -t <p> -m <mem>000000000
```

| 参数 | 含义 |
|---|---|
| `-m <mem>000000000` | 以**字节**给出 SdBG 构建阶段的内存上限（MitoFinder 把 GB 后面补 9 个 0） |
| `-t` | 线程（`rule mitofinder` `threads: 8`） |
| `--min-contig-len` | 最短输出 contig（MEGAHIT 默认 200，MitoFinder 改为 100） |
| k 列表 | 未指定 → 默认 `21,29,39,59,79,99,119,141`，大于"最大读长 + 20"的 k 会被自动去掉（`megahit` 脚本源码） |

可调键：`mitofinder_assembly`（given / spades / megahit / metaspades）、`mitofinder_mem_gb`（默认 32；`mem_mb` 重试 ×1.5）。
MitoFinder 已有 `<s>_megahit` 目录时会跳过 MEGAHIT（"Mitofinder will skip MEGAHIT step"）；wgs_skim 每次先 `rm -rf` 工作目录，所以总是重新组装。

## 4. 重要参数说明（`megahit -h`，v1.2.9）

| 参数 | 默认 | 说明 |
|---|---|---|
| `--k-list` / `--k-min --k-max --k-step` | 21,29,39,59,79,99,119,141 | 奇数，15–255，步长 ≤28；`--k-min/max/step` 会覆盖 `--k-list` |
| `--presets meta-sensitive / meta-large` | — | **会覆盖 `--k-list`**（#313） |
| `--min-count` | 2 | (k_min+1)-mer 最低多重度；低覆盖（线粒体之外的核基因组）不影响 mt，高覆盖 mt 足够 |
| `--no-mercy`、`--prune-level`、`--bubble-level` | — / 2 / 2 | 图清理强度 |
| `-m/--memory` | 0.9 | 字节数，或 0–1 表示占机器总内存比例；**只约束 SdBG 构建阶段**（作者 #175） |
| `--mem-flag` | 1 | SdBG 构建内存模式：0 最少，1 中等，其他 = 用满 `-m` |
| `-t` | 全部逻辑核 | 线程 |
| `--tmp-dir` | 输出目录 | 临时目录；输入 .gz 放在 CIFS 等网络盘上时 `mkfifo` 失败，改到本地盘（#303） |
| `--continue` | — | 从最后检查点继续（需相同 `-o`） |
| `--keep-tmp-files` | 关 | 保留中间文件 |
| `--no-hw-accel` | — | 不用 BMI2/POPCNT（旧 CPU 报非法指令时） |

## 5. 已知问题与 issues

| 现象 / 报错 | 原因 | 解决办法 | 链接 |
|---|---|---|---|
| `Exit code -9` | 进程被 SIGKILL：作者——超出容器 / 调度器内存限制，或被他人杀掉 | 提高内存（`mitofinder_mem_gb`），查调度器日志 | https://github.com/voutcn/megahit/issues/272 、 https://github.com/voutcn/megahit/issues/379 |
| 内存超过 `-m` | `-m` 只限制 SdBG 构建，其他阶段不受限（作者） | 调度器内存要比 `-m` 留余量 | https://github.com/voutcn/megahit/issues/175 |
| `Error occurs when reading inputs` | 一位用户：fastq 文件传输中损坏，重新生成后解决；另有多文库时出现的未解决报告 | `gzip -t`、检查 R1/R2 读数一致 | https://github.com/voutcn/megahit/issues/271 |
| `.gz` 输入在服务器上报错、解压后正常 | 网络文件系统上无法创建命名管道 | `--tmp-dir /本地盘` | https://github.com/voutcn/megahit/issues/303 |
| 指定的 k 不生效 | `--presets` 覆盖了 `--k-list` | 不要同时用 preset | https://github.com/voutcn/megahit/issues/313 、 https://github.com/voutcn/megahit/issues/367 |
| `Exit code -11` / 段错误 | 多种（1.2.8 的零长度序列 #244 已在 1.2.9 修复；其余未定） | 确认 1.2.9；去掉零长度 reads（fastp `-l`） | https://github.com/voutcn/megahit/issues/244 、 https://github.com/voutcn/megahit/issues/310 |
| `Exit code -6` / `-7`（大型组装） | 断言失败 / 未定，作者索要日志 | 附 `log` 文件报 issue | https://github.com/voutcn/megahit/issues/222 、 https://github.com/voutcn/megahit/issues/363 |
| `--continue` 似乎不续跑 | 未有维护者答复 | 重新组装（mt 数据量下 MEGAHIT 很快） | https://github.com/voutcn/megahit/issues/376 |
| `Output directory ... already exists, please change the parameter -o` | MEGAHIT 拒绝覆盖已有输出目录（`megahit` 脚本源码） | 删除目录或用 `--continue` | 源码，非 issue |
| 多线程在 macOS 失败 | 平台问题 | Linux 上无关 | https://github.com/voutcn/megahit/issues/385 |

## 6. 本项目中遇到的问题

- TESTING.md §4 与 README Troubleshooting 中没有 MEGAHIT 的记录；MitoFinder（含 MEGAHIT）的 Snakemake 规则仍在开发，尚未在本仓库端到端测试。
- 注意：经 MitoFinder 调用时，MEGAHIT 的日志在 MitoFinder 工作目录内（`<s>_megahit/log`），而 wgs_skim 在规则结束时 `rm -rf` 整个工作目录——
  失败时日志可能随之消失（规则在 MitoFinder 失败时 `exit 1`，在 `rm -rf` 之前退出，所以失败时目录保留；成功时删除）。
- 对动物 genome skimming 的取舍：MEGAHIT 对整份 1–20 Gb 数据组装核基因组 + 线粒体，内存远低于 SPAdes，但 MitoFinder 只要线粒体 contig；
  已有 SPAdes 组装时应优先 `mitofinder_assembly: spades`（用户规则：不重复组装）。

## 7. 排错清单

1. `logs/mitofinder/<s>.log`（MitoFinder 输出）；失败时 MitoFinder 工作目录 `mitofinder/<s>/run/` 下的 `*_megahit/log`（MEGAHIT 自己的日志，含每个 k 的命令与退出码；确切子路径未实测）。
2. 关键字：`Exit code -9`（被杀，内存）、`Exit code -11`（段错误）、`Error occurs when reading inputs`、`already exists`、
   `Memory used`、`k-list`（实际使用的 k）。
3. `benchmarks/mitofinder/<s>.tsv` 的 `max_rss` 与 `mitofinder_mem_gb` 对比。
4. 输入：`gzip -t` 两个 clean 文件；`fastp.json` 的 R1/R2 读数一致。
5. 版本：`megahit --version`（应为 `MEGAHIT v1.2.9`）。

## 8. 资料来源与查阅日期（2026-10-29）

- Releases：https://github.com/voutcn/megahit/releases ；commits（`gh api repos/voutcn/megahit/commits`）。
- 本机 `megahit -h` 与 `megahit` 脚本源码（`~/miniforge3/envs/megahit`，v1.2.9）。
- MitoFinder 1.4.2 `runMegahit.py`（`~/mybiosofts/MitoFinder`）。
- 上表 issues（2026-10-29 查阅）；wgs_skim `workflow/rules/mitofinder.smk`、`workflow/envs/mitofinder.yaml`、`config/config.yaml`（工作区版本）。

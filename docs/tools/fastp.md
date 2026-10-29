# fastp（wgs_skim 排错参考）

## 1. 简介与用途

fastp 是 FASTQ 质控 / 接头去除 / 过滤一体化工具（C++，多线程）。在 wgs_skim 中：

- `rule fastp`（`workflow/rules/qc.smk`）是所有下游步骤的入口：输出 `qc/<s>/<s>.R1.fq.gz`、`.R2.fq.gz`、`.fastp.html`、`.fastp.json`。
- `fastp.json` 的 `summary.after_filtering.total_bases` 被 `common.smk::clean_gb()` 读取，换算成"clean Gb"，决定
  SPAdes（`spades_gib_per_gb`）、其他 SPAdes 模式、Trinity（`trinity_gib_per_gb`）的内存估计，以及 `spades_min_clean_gb` 是否跳过组装、
  GetOrganelle `--max-reads` 的换算（`mt_iterate.py` 用 fastp 的平均读长）。**json 不对 → 下游内存全部算错。**
- Nextflow：`nextflow/modules/qc.nf` 的 `FASTP`（8 CPU、2000 MB）。

链接：
- GitHub：https://github.com/OpenGene/fastp （README 即手册）
- 论文：Chen et al. 2018 *Bioinformatics* 34(17):i884–i890，doi:10.1093/bioinformatics/bty560；
  Chen 2023 *iMeta* 2:e107，doi:10.1002/imt2.107（fastp 超快版 / 新功能）

## 2. 版本

| | 版本 | 发布日期（GitHub release） |
|---|---|---|
| wgs_skim 固定 | `fastp =0.24.1`（`workflow/envs/fastp.yaml`） | 2025-04-16 |
| 最新 | **1.4.0** | 2026-10-07 |

0.24.1 之后的 release 说明（原文要点，`gh api repos/OpenGene/fastp/releases`）：

| 版本 | 变化 |
|---|---|
| 0.24.2（2025-05-27） | **修复"记录接头序列导致的内存爆炸"（#602）**；`--version` 输出到 stdin(原文如此) |
| 0.24.3 | 多线程 stdout 模式输出可重复 |
| 0.25.0 / 0.26.0 | 并排 HTML 报告；改进接头检测与修剪 |
| 1.0.0 / 1.0.1 | 里程碑版；修复 `--overlapped_out` 段错误（#616） |
| 1.1.0（2026-01） | 修复交错输出可重复性回归；"按修剪后长度识别 adapter dimer"（新增 `--dimer_max_len`） |
| 1.2.0 / 1.3.0（2026-03） | SIMD 加速；并行压缩 |
| 1.3.1–1.3.7 | 一连串修复：merge 回归、线程过多时可能挂起（1.3.2/1.3.3）、BGZF 读取挂起（1.3.5）、限制 worker 线程（1.3.6）、I/O（1.3.7） |
| 1.4.0 | 修复挂起问题 #721（自动接头检测后挂起；与 >32 线程挂起 #694 / PR #723 同源） |

`fastp --help` 0.24.1 与 1.4.0 对比（本机 `~/miniforge3/envs/fastp`、`envs/wgs_fastp` 实测）：选项几乎相同；1.4.0 新增
`--allow_gap_overlap_trimming`、`--dimer_max_len`，`--detect_adapter_for_pe` 有了短名 `-2`，STDIN 模式不做接头自动检测。

**是否升级：建议至少升到 0.24.3（或直接 1.4.0）。** 理由：
1. 0.24.1 带有 #602 的内存问题（短插入片段文库内存可涨到 100+ GB，维护者在 0.24.2 修复）。动物 genome skimming 的
   降解样品 / 博物馆标本常有短片段，而 wgs_skim 只给 fastp 声明 `mem_mb=2000`——这是 Snakemake 调度值，不是硬限制，
   内存真涨上去会挤掉同机其他任务（参见 TESTING.md §4 #16）。
2. 1.3.x 系列多次出现线程/挂起回归（#685 `-w` 不再限制总 CPU，仍 open；#694、#721 挂起，1.4.0 关闭）。若升 1.x，
   直接用 1.4.0，不要用 1.3.x。升级后需重跑 e2e 比较 clean Gb（接头检测在 0.26.0 改过，clean 数据量可能略变）。

## 3. wgs_skim 中的调用方式

`qc: true`（默认）：

```bash
fastp -i {r1} -I {r2} -o qc/<s>/<s>.R1.fq.gz -O qc/<s>/<s>.R2.fq.gz --detect_adapter_for_pe \
      -w {threads} -h <s>.fastp.html -j <s>.fastp.json > logs/fastp/<s>.log 2>&1
```

`qc: false`（reads 已经是 clean，`wgs_skim spades ...` 等单工具命令用）：原始文件被软链接为 qc 输出，fastp 只做计数：

```bash
fastp -i {r1} -I {r2} -A -G -Q -L -w {threads} -h ... -j ...
```

| 参数 | 含义 |
|---|---|
| `-i/-I`、`-o/-O` | 双端输入 / 输出（.gz 自动压缩，默认压缩级别 `-z 4`） |
| `--detect_adapter_for_pe` | 双端也做接头序列自动检测（默认双端只靠 overlap 分析去接头） |
| `-w` | worker 线程数，默认 3；wgs_skim 规则 `threads: 8`（`--set-threads fastp=N` 可改） |
| `-h/-j` | HTML / JSON 报告；JSON 是下游内存估计的依据 |
| `-A -G -Q -L` | 关闭接头修剪 / polyG 修剪 / 质量过滤 / 长度过滤（只统计） |

其余全部用默认：`-q 15 -u 40 -n 5 -l 15`，polyG 修剪对 NovaSeq/NextSeq 数据自动开启，不去重（`-D` 未开）。
没有可调的 fastp 配置键；只有 `qc`（true/false）。改参数需要改规则（或今后加一个 `fastp_extra` 键——目前不存在）。

## 4. 重要参数说明（`fastp --help` 0.24.1 原文摘译）

| 参数 | 默认 | 说明 / 与我们的关系 |
|---|---|---|
| `-w, --thread` | 3 | worker 线程。1.3.2 起 `-w` 不约束总 CPU（#685，open）；>32 线程在 1.3.x 可能挂起（#694，1.4.0 修复） |
| `--detect_adapter_for_pe` | 关 | 1.x 帮助："enable adapter detection for PE data to get ultra-clean data. It takes more time…" |
| `-a / --adapter_fasta` | — | 指定接头序列 / FASTA（已知文库接头时比自动检测稳） |
| `-g / -G` | 自动 | 强制 / 关闭 polyG 修剪；大数据时 polyG 未修剪的报告见 #464、#731（open） |
| `-x, --trim_poly_x` | 关 | polyX 尾修剪 |
| `-l, --length_required` | 15 | 短于此长度丢弃；组装用数据常设 30–50（SPAdes 最小 k=21） |
| `-q / -u / -n` | 15 / 40 / 5 | 合格碱基阈值、不合格碱基比例上限、N 个数上限 |
| `-5 / -3 / -r` 与 `-W -M` | 关；4 / 20 | 滑窗质量修剪（cut_front / cut_tail / cut_right） |
| `-c, --correction` | 关 | 双端 overlap 区碱基纠错 |
| `-D, --dedup` / `--dup_calc_accuracy` | 关 / 1 | 去重；精度等级 1–6 对应 1G…24G 内存（1.4.0 为 32G） |
| `-z` | 4 | gzip 压缩级别，影响磁盘与速度 |
| `--overlap_len_require` | 30 | overlap 检测最小长度（影响接头修剪） |

## 5. 已知问题与 issues

| 现象 / 报错 | 原因 | 解决办法 | 链接 |
|---|---|---|---|
| 内存暴涨到数十～100+ GB（短片段文库） | 记录接头序列造成内存爆炸（0.24.0/0.24.1） | 升级 ≥0.24.2 | https://github.com/OpenGene/fastp/issues/602 |
| 内存不足时不报 OOM 而是一直挂着直到超时（13–27 Gb 单端文件） | 维护者确认要修，仍 open | 给足内存（报告中 12 GB 可跑完）；设置任务时限 | https://github.com/OpenGene/fastp/issues/657 |
| 打印 "No adapter detected for read2" 后挂起，输出 0 字节（1.3.7，`-w 35`） | 线程相关死锁 | 升 1.4.0；或线程 ≤16 | https://github.com/OpenGene/fastp/issues/721 |
| R1/R2 读数不一致时警告 "different read numbers of the pack" 后挂起 | 同上（>32 线程触发） | 1.4.0 修复；先检查 R1/R2 是否截断 | https://github.com/OpenGene/fastp/issues/694 |
| `-w 8` 实际占用远超 8 核（1.3.2） | 1.3.x 新增并行压缩等辅助线程 | 仍 open；调度时多留 CPU，或用 ≤1.0.x | https://github.com/OpenGene/fastp/issues/685 |
| `--overlapped_out` 段错误（1.0.0） | bug | 1.0.1 修复 | https://github.com/OpenGene/fastp/issues/616 |
| 大数据 polyG 尾未被修剪 | 未定（open） | 加 `-g` 或 `--trim_poly_x`；检查 HTML 中 polyG 统计 | https://github.com/OpenGene/fastp/issues/464 、 https://github.com/OpenGene/fastp/issues/731 |
| HTML 去重率直方图消失 / 去重率解读 | 0.22.0 起报告变化 | 看 json `duplication.rate` | https://github.com/OpenGene/fastp/issues/437 |
| `error while loading shared libraries: libisal.so.2` | 直接运行 conda **pkgs 缓存**里的二进制（无依赖库） | 只用激活后的环境（`conda run -p <env> fastp`） | 本机实测，非 issue |

## 6. 本项目中遇到的问题

- TESTING.md §4 中没有 fastp 本身的错误记录。实测峰值内存 0.4 GB（0.6 Gb 样品，TESTING §3.4）。
- 相关：#16（机器内存紧张导致后台任务被停）——fastp 内存 `mem_mb=2000` 只是调度值；若遇 #602 类内存暴涨不会被 Snakemake 限住。
- README Troubleshooting 中与 fastp 有关的只有输入识别（`no paired FASTQ files ... found`、`sample 'X' found twice`），
  那是 `samples.py` 的问题，不是 fastp。

## 7. 排错清单

1. `logs/fastp/<s>.log`：最后应有 `JSON report: ...`、`HTML report: ...` 和 `fastp v0.24.1, time used: N seconds`。没有 → 被杀 / 挂起。
2. 关键字：`ERROR:`、`different read numbers`（R1/R2 不同步，常因上传截断）、`gzip`/`igzip` 错误（文件损坏，先 `gzip -t`）、
   `Segmentation fault`、`Detecting adapter sequence` 后长时间无输出（挂起，见 #721）。
3. `benchmarks/fastp/<s>.tsv`：`max_rss` 若远大于 2 GB → 怀疑 #602（看 log 中接头检测部分、插入片段分布）。
4. `qc/<s>/<s>.fastp.json`：`summary.before_filtering` vs `after_filtering` 的 `total_reads/total_bases`、
   `adapter_cutting`、`duplication.rate`、`insert_size.peak`。`after_filtering.total_bases` 异常小 → 下游内存估计偏低，SPAdes 易 OOM。
5. 版本核对：`conda run -p <env> fastp --version 2>&1`（0.24.2 起改了版本号的输出流，合并 stderr 最稳）。

## 8. 资料来源与查阅日期（2026-10-29）

- GitHub releases（`gh api repos/OpenGene/fastp/releases`）与上表各 issue 页面，2026-10-29 查阅。
- README：https://github.com/OpenGene/fastp#readme
- 本机 `fastp --help`：0.24.1（`~/miniforge3/envs/fastp`）与 1.4.0（`~/miniforge3/envs/wgs_fastp`）。
- wgs_skim：`workflow/rules/qc.smk`、`workflow/envs/fastp.yaml`、`workflow/rules/common.smk`、`docs/TESTING.md`、`README.md`（工作区版本，2026-10-29）。

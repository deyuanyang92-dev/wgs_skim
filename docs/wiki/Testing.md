# 测试

完整的测试记录（数据、每项结果、发现并修复的 16 个错误、尚未测试的部分）：仓库中的
[docs/TESTING.md](https://github.com/deyuanyang92-dev/wgs_skim/blob/main/docs/TESTING.md)。本页说明有哪些测试、如何运行。

## 三层测试

| 测试 | 命令 | 时间 | 需要 |
|---|---|---|---|
| 单元测试（16 项：样品识别、环化、GetOrganelle 轮次选择、核基因比较、BLAST 解析、汇总……） | `python -m pytest -q tests/` | < 1 秒 | Python、Biopython |
| 命令矩阵（33 项：每个子命令、每个 `--steps` 值与 SPAdes 模式的依赖图、错误输入、内存选项、对已完成结果运行 circularize / compare / nuclear） | `tests/cli_matrix.sh 小数据/ 已完成结果/` | 约 1 分钟 | 已安装的 wgs_skim |
| 端到端（2 个以上样品完整运行，逐样品检查 11 项输出 + BLAST、汇总、报告、`wgs_skim nuclear`、单元测试） | `tests/e2e.sh 小数据/ 输出/` | 约 30 分钟 | 已安装的 wgs_skim 和数据库 |

每个脚本的退出码 = 失败项数；`cli_matrix.sh` 写出每项检查的命令、期望、结果、耗时（`cli_matrix.tsv`）。

## 测试数据

不需要大数据：从任意样品取 1 M 对 reads（约 0.3 Gb）即可检查流程是否完整跑通。

```bash
mkdir small
for s in A B; do for r in 1 2; do
    zcat raw/$s*_R$r*.fq.gz | head -n 4000000 | gzip > small/${s}_R$r.fq.gz
done; done
tests/e2e.sh small/ e2e_out/
tests/cli_matrix.sh small/ e2e_out/
```

数据量太小时线粒体不会闭合（`linear_low_reads`）、BUSCO 接近 0 —— 这是预期的；小数据只检查功能，不评价组装质量。

## 最近一次结果（2026-10-29）

| 测试 | 结果 |
|---|---|
| 单元测试 | 16 / 16 通过 |
| 命令矩阵 | 33 / 33 通过 |
| 端到端（pisione + cwsc010，各 0.3 Gb） | 39 / 39 个任务成功，0 项检查失败，26 分钟 |
| 单样品完整流程（HNP1sub，0.6 Gb） | 21 / 21 个任务，28 分钟；rDNA 5 个基因两种方法完全相同 |
| 核基因批量（`wgs_skim nuclear`，Snakemake） | 8 个任务成功；再次运行 "Nothing to be done"（2 秒） |
| 线粒体比较批量（`wgs_skim compare --results`） | 4 个样品 110 秒：identical / consistent / consistent / skipped（缺 MitoFinder） |
| 真实批量 MJ20260820041（22 个样品） | rDNA 注释全部通过 NCBI table2asn（0 错误） |
| 转录组（rnaSPAdes + Trinity + BUSCO，2 个样品各 0.1 Gb） | 见下 |

### 转录组测试

数据：两个 RNA-seq 样品（A1、R1），各约 50 MB × 2。命令 `wgs_skim run -r rna2 -o rna2_out --steps rna`。

| 样品 | rnaSPAdes 转录本 | rnaSPAdes 时间 / 峰值内存 | BUSCO（metazoa，transcriptome） |
|---|---|---|---|
| A1 | 20 573 | 87 s / 1.3 GB | C:12.1% [S:11.9%, D:0.1%], F:22.8%, M:65.2% |
| R1 | 27 620 | 92 s / 1.2 GB | C:13.2% [S:12.9%, D:0.3%], F:29.8%, M:57.0% |

BUSCO 低是因为数据只有约 0.1 Gb（完整转录组通常需要 ≥ 2 Gb）；本测试检查流程。Trinity 结果见 docs/TESTING.md 第 9 节。

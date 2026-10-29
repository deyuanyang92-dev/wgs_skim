# Jellyfish + GenomeScope2（wgs_skim 排错参考）

## 1. 简介与用途

基因组调查（`genome_size: true`，`workflow/rules/genome_size.smk`；Nextflow `JELLYFISH` / `GENOMESCOPE` in
`nextflow/modules/assembly.nf`）：

1. `rule jellyfish`：对 fastp 后的 reads 数 k-mer（canonical），输出直方图 `genome_size/<s>/<s>.k21.histo`，删掉 `.jf`。
2. `rule genomescope`：GenomeScope2 拟合 k-mer 谱，估计单倍体基因组大小、杂合度、重复比例、错误率；摘要复制为
   `genome_size/<s>/<s>.genomescope.txt`，汇总进 `summary/samples.tsv` 的 `genome_size_mb`、`heterozygosity_pct`。

链接：
- Jellyfish：https://github.com/gmarcais/Jellyfish ；Marçais & Kingsford 2011 Bioinformatics 27(6):764, doi:10.1093/bioinformatics/btr011
- GenomeScope2：https://github.com/tbenavi1/genomescope2.0 （**该仓库已关闭 issues**；问题都在 GenomeScope 1 的仓库
  https://github.com/schatzlab/genomescope/issues ）；网页版 http://genomescope.org/genomescope2.0/
- 论文：Ranallo-Benavidez, Jaron & Schatz 2020 Nat Commun 11:1432, doi:10.1038/s41467-020-14998-3；
  GenomeScope 1：Vurture et al. 2017 Bioinformatics, doi:10.1093/bioinformatics/btx153

## 2. 版本

| | 固定（`workflow/envs/genome_size.yaml`） | 最新 |
|---|---|---|
| Jellyfish | `kmer-jellyfish =2.3.1` | 2.3.1（2023-12-07，仍是最新） |
| GenomeScope2 | `genomescope2 =2.1.0` | 2.1.0（2025-10-02） |

- 2.1.0 相对 2.0.1：代码整理（去全局变量）、安装脚本改进、新增 **`--json_report`**（在文本摘要之外写 JSON 报告），完全向后兼容。
  注意：2.1.0 的日志 / 摘要第一行仍写 `GenomeScope version 2.0`（e2e2 实测），不能用它判断安装版本。
- 是否升级：已是最新。建议**利用 `--json_report`**，让 `summary.py` 读 JSON 而不是解析文本表（文本列宽固定、数值带 `%`/`bp`）。
- conda 包名：Jellyfish 在 bioconda 叫 `kmer-jellyfish`；PyPI / conda-forge 上的 `jellyfish` 是一个字符串相似度库，不是它
  （另见 Jellyfish #161 关于 Python 绑定的混乱）。

## 3. wgs_skim 中的调用

```bash
jellyfish count -C -m 21 -s 2G -t 16 -o genome_size/<s>/reads.jf <(zcat R1.fq.gz) <(zcat R2.fq.gz)
jellyfish histo -t 16 genome_size/<s>/reads.jf > genome_size/<s>/<s>.k21.histo
rm -f genome_size/<s>/reads.jf
genomescope2 -i <s>.k21.histo -o genome_size/<s>/genomescope -k 21 -p 2 -n <s> \
  || echo "[wgs_skim] GenomeScope did not converge (low coverage?)"
# 有 <s>_summary.txt 就复制，否则写 "GenomeScope failed, see log"（规则本身不失败）
```

- `-C` canonical（双链），必须；`-m 21` = `kmer_size`；`-s 2G` 初始哈希槽数（不是内存、也不是输入大小，满了会自动翻倍 / 落盘合并）；
  `-t 16` 线程。Jellyfish 不能读 gz，用进程替换喂 `zcat`（官方文档推荐的做法）。
- `jellyfish histo` 默认 `-h 10000`：≥ 10000 次的 k-mer 全部堆在最后一行（e2e2：`10001 94`）。
- `-p 2` = `ploidy`；`-n` 输出文件前缀。资源：jellyfish 16 线程、`mem_mb: 12000`；GenomeScope 1 线程、2000 MB。
- 可配置键：`genome_size`（开关）、`kmer_size`（21）、`ploidy`（2）、`steps: genome_size`。

## 4. 重要参数

| 工具 | 参数 | 含义 / 建议 |
|---|---|---|
| jellyfish count | `-m` | k；官方推荐 21（Schatz：约 90 % 基因组唯一，又不至于让错误 k-mer 和杂合 k-mer 糊在一起） |
| | `-s` | 哈希槽数，约等于预期不同 k-mer 数；人类 reads 3G 足够（#198 回答）。设太大 → 一开始就申请巨大内存 |
| | `-C` | canonical，必须 |
| | `-L` / `-U` | 计数上下限（不建议在 count 时截） |
| jellyfish histo | `-h` | 上限（默认 10000）。重复极高 / 巨大基因组可提高到 1e6，但更易受细胞器 / phiX 影响 |
| genomescope2 | `-k` | 必须与 count 的 `-m` 一致 |
| | `-p` | 倍性（1–6） |
| | `-l` | 初始 λ（杂合 k-mer 覆盖）猜测；低覆盖难收敛时可手动给 |
| | `-m/--max_kmercov` | 排除高频 k-mer（细胞器、phiX、污染） |
| | `--json_report`（2.1.0） | JSON 摘要 |
| | `--fitted_hist`、`--testing` | 输出拟合直方图 / 调试信息 |

## 5. 已知问题与 issues

| 现象 / 报错 | 原因 | 解决 | 链接 |
|---|---|---|---|
| `Failed to converge` / 模型不收敛 | 覆盖太低：**单倍体覆盖至少 15×**，否则拟合差或不收敛（README） | 先确认 `-C`；再试 k=19/17；最终只能加数据 | README FAQ；https://github.com/schatzlab/genomescope/issues/120 、/127 |
| "收敛"但基因组大小离谱（小一半或大一倍） | 低覆盖时模型把纯合峰和杂合峰混淆 | 不要信低覆盖结果 | https://github.com/schatzlab/genomescope/issues/146 |
| 模型 fit 很低 / 拟合怪 | 数据质量、污染、高重复 | 看图；`--max_kmercov` 去高频 | /161 /139 /152 |
| 无法自动判断"结果不可信" | GenomeScope 本身没有低覆盖失败状态（开发者承认，功能请求未实现） | 流程里自己加判据 | https://github.com/schatzlab/genomescope/issues/148 |
| 杂合峰被当作错误峰 | 低覆盖 / 高错误率 | 给 `-l`；加数据 | /131 /173 |
| 倍性不确定 | GenomeScope 不估倍性 | 用 Smudgeplot | /142 /130；README FAQ |
| CLI 与网页结果不同 | 参数默认值不同 | 统一参数 | /172 |
| k 的选择 | 21 是测试最多的；巨大（≫10 Gb）/ 高重复基因组可用更大 k | 保持 21 | /160 |
| jellyfish `terminate called ... large_hash` / 内存不足 | `-s` 太大 | 按不同 k-mer 数估 `-s` | https://github.com/gmarcais/Jellyfish/issues/198 、/177 |
| 内存分配失败后留下孤儿 generator 进程 | Jellyfish bug | 手动杀掉残留 `zcat`/generator | https://github.com/gmarcais/Jellyfish/issues/162 |
| 不能读 .gz | 设计如此 | `zcat` 管道 / 进程替换 | https://github.com/gmarcais/Jellyfish/issues/186 |

## 6. 本项目中遇到的问题

- **低覆盖下给出"收敛"的荒谬结果，并写进 summary**（e2e2 两样品，各 0.3 Gb clean）：
  `Model converged het:0.0372 kcov:35.5 err:0.114 model fit:2.3 len:282719` → 单倍体基因组 0.28 Mb、杂合度 3.7 %、
  错误率 11.4 %。环节动物基因组通常 0.3–2 Gb，0.3 Gb 数据只有 < 1×；kcov 35 的"峰"实际来自线粒体 / rDNA 等多拷贝序列。
  规则里的 `|| echo "did not converge"` 抓不到这种情况（程序正常退出）。**建议**：summary 中对以下情况标记为不可信而不是报数：
  错误率 > 5 %、或估计长度 < 10 Mb、或 `clean_gb × 1e9 / 估计长度 < 15`（覆盖不足 15×），并在 `genome_size_mb` 旁写原因。
- 典型 genome skimming 数据（每样品 5–10 Gb、基因组约 1 Gb）也只有 5–10×，低于 15× 门槛：对多数样品 GenomeScope 只能当参考，
  可考虑默认关闭 `genome_size` 或只在 clean data 足够时运行（类似 `spades_min_clean_gb`）。
- 资源实测（e2e2）：jellyfish 2 min、峰值 6.8 GB（主要是 `-s 2G` 的初始哈希）；GenomeScope 31 s。
- GenomeScope 日志首行 `[1] "已设置北大阿里云镜像"` 来自用户 `~/.Rprofile`，无害；若要干净日志可在规则中设 `R_PROFILE_USER=/dev/null`（未测试）。

## 7. 排错清单

1. `logs/jellyfish/<s>.log`：`terminate called`、`std::bad_alloc`、`Can't open file` → `-s`、内存、输入路径。
2. `<s>.k21.histo`：第一行（count=1 的错误 k-mer）应远大于其他；看是否存在第二个峰。没有峰 = 覆盖太低。
3. `logs/genomescope/<s>.log`：`Model converged ... kcov:X err:Y model fit:Z len:L`；`Failed to converge`。
4. `genome_size/<s>/genomescope/<s>_linear_plot.png`、`_log_plot.png`：拟合线是否穿过真实峰。
5. 判定可信：kcov（杂合 k-mer 覆盖，纯合峰 = 2×kcov）≥ ~7.5（即纯合 ≥ 15×）、err < ~1 %（Illumina）、长度与近缘种 C-value 同量级。
6. 结果离谱时：同样数据降 k 到 19/17 再跑，看是否一致（`--config kmer_size=19`）。

## 8. 资料来源（查阅日期 2026-10-29）

- https://github.com/tbenavi1/genomescope2.0 README 与 releases（v2.1.0、v2.0.1）
- https://github.com/schatzlab/genomescope/issues/120 /127 /131 /146 /148 /152 /160 /161 /172 /173
- https://github.com/gmarcais/Jellyfish releases 与 issues #162 #177 #186 #198
- 本项目：`workflow/rules/genome_size.smk`、`nextflow/modules/assembly.nf`、`../wgs_skim_test/e2e2/genome_size/`、
  `../wgs_skim_test/e2e2/logs/genomescope/`、`benchmarks/jellyfish/`

# 各软件资料：版本、参数、已知问题（排错用）

每页：用途与官方链接 → 固定版本与最新版本的差别 → wgs_skim 中的实际命令 → 重要参数 → 已知问题与 GitHub issues
（现象 / 原因 / 解决 / 链接）→ 本项目遇到的问题 → 排错清单 → 资料来源。查阅日期 2026-10-29。

| 软件 | 页面 | wgs_skim 固定版本 |
|---|---|---|
| fastp | [fastp.md](fastp.md) | 0.24.3 |
| MEANGS | [MEANGS.md](MEANGS.md) | 1.2.1 |
| GetOrganelle | [GetOrganelle.md](GetOrganelle.md) | 1.7.7.1（内部 SPAdes 3.15.5） |
| SPAdes（全部模式） | [SPAdes.md](SPAdes.md) | 4.2.0 |
| Trinity | [Trinity.md](Trinity.md) | 2.15.2（build pl5321h077b44d_3） |
| MEGAHIT（MitoFinder 内部） | [MEGAHIT.md](MEGAHIT.md) | 1.2.9 |
| MitoFinder | [MitoFinder.md](MitoFinder.md) | 1.4.2（Singularity）/ 1.4.1（bioconda） |
| MitoM | [MitoM.md](MitoM.md) | 固定提交（`mitom_commit`、`mitom_nuc_commit`） |
| BUSCO | [BUSCO.md](BUSCO.md) | 6.0.0 |
| BLAST+ | [BLAST.md](BLAST.md) | 2.16.0 |
| Jellyfish + GenomeScope2 | [GenomeScope.md](GenomeScope.md) | 2.3.1 / 2.1.0 |
| Snakemake | [Snakemake.md](Snakemake.md) | 8.18 |
| Nextflow | [Nextflow.md](Nextflow.md) | 25.x |

## 研究中发现、已改进的问题（2026-10-29）

| 问题 | 证据 | 改动 |
|---|---|---|
| Trinity 第二阶段串行运行 | conda build _6 的 ParaFly 忽略 `-CPU`（trinityrnaseq#1473）；8 个测试任务：_6 用 16 s，_3 用 2 s | 固定 build `pl5321h077b44d_3` |
| fastp 0.24.1 内存暴涨 | OpenGene/fastp#602（短插入片段文库可达 100+ GB），0.24.2 修复 | 0.24.3 |
| SPAdes 续跑可能混用旧数据 | `--restart-from` 沿用旧的 reads 与参数（手册：续跑时不能改输入） | 记录 reads + 参数签名：相同才续跑，否则从头；完成后删除签名 |
| rnaSPAdes `-k 21,33,55` | 手册：小 k 产生嵌合转录本，建议让 SPAdes 按读长自选 2 个 k | 默认 `auto`（旧值可配置） |
| MitoFinder 1.4.1 不检查环化 | 源码比对：检查脚本启动时没有传参数；1.4.2 发布说明承认 | 默认用 1.4.2 镜像（`--cleanenv`）；1.4.1 结果注明 |
| MitoFinder 在整个基因组组装上慢、环化检查不运行 | 旧脚本 spades_mitofinder.py：覆盖度预筛后 CHY006 只剩 1 条 contig | 恢复预筛（≥ 1 kb、≥ 5 × 加权中位覆盖度） |
| MitoFinder 自组装用全部 reads | 作者建议 5–10 M 对（MitoFinder#13） | `mitofinder_max_pairs` = 1000 万 |
| MitoM 静默降级 | 环境缺 matplotlib、scikit-learn：不画图、tRNA 退回 MiTFi | 加入 `envs/mitom.yaml` |
| GenomeScope 低覆盖下给出荒谬结果 | e2e2：0.28 Mb、错误率 11 % 仍 "converged" | 汇总表 `genome_size_note` 标为 unreliable |
| Nextflow ≥ 26.04 无法解析 | 严格语法解析器默认开启 | 运行时设 `NXF_SYNTAX_PARSER=v1` |

## 尚待验证（需要真实数据比较后再改）

| 建议 | 来源 | 状态 |
|---|---|---|
| MEANGS `--skipextend` 或用 `*_deep_detected_mito.fas` 作种子（省去最耗内存的延伸） | MEANGS.md，MEANGS#20 | 待与 GetOrganelle 结果比较 |
| GetOrganelle 1.8.0.1（链接 reads、修复 anonym 短线性结果 #239） | GetOrganelle.md | 未发布到 bioconda，有报错（#388、#390） |
| BLAST 跨种 COI 用 `-task dc-megablast` | BLAST.md | 待比较命中 |
| 环节动物 / 软体动物 BUSCO 用 `lophotrochozoa_odb12`（1252 个基因） | BUSCO.md | 待真实样品比较 |
| SPAdes 4.3.0 `--frugal` 作为内存不足时的最后重试 | SPAdes.md | 待比较 N50、BUSCO |
| fastp 1.4.0 | fastp.md | 接头识别改变，需重新核对 clean 数据量 |

# Trinity（wgs_skim 排错参考）

## 1. 简介与用途

Trinity 是 RNA-seq de novo 转录组组装器：Phase 1 = in silico 归一化 → Jellyfish 计 k-mer → Inchworm → Chrysalis 把 reads 分簇；
Phase 2 = 对每个 read 簇单独运行一次 Trinity（Butterfly），由 ParaFly 并行调度；最后汇总（genome-free 模式下再做 salmon 过滤）。
在 wgs_skim 中：

- `rule trinity`（`workflow/rules/transcriptome.smk`）→ `assembly/<s>/<s>.trinity.fasta`（+ `.gene_trans_map`），
  `rule busco_assembly` 用 `-m transcriptome` 评估。`--steps trinity` 或 `--steps rna`（rnaSPAdes + Trinity + BUSCO）。
- 面向 RNA-seq 样品（来自用户旧脚本 `trinity_batch3.1.py`）；不用于 DNA skim 数据。Nextflow 版不支持 Trinity。

链接：
- GitHub：https://github.com/trinityrnaseq/trinityrnaseq ；Wiki：https://github.com/trinityrnaseq/trinityrnaseq/wiki
  （Running-Trinity、minimize_RAM_usage、grid_conf）
- 论文：Grabherr et al. 2011 *Nat Biotechnol* 29:644，doi:10.1038/nbt.1883；Haas et al. 2013 *Nat Protoc* 8:1494，doi:10.1038/nprot.2013.084

## 2. 版本

| | 版本 | 日期 |
|---|---|---|
| wgs_skim 固定 | `trinity =2.15.2`（`workflow/envs/trinity.yaml`；实际解析到 bioconda build `pl5321h077b44d_6`） | 2024-07-31 |
| 最新 | **2.15.2**（无更新的 release；主分支最后提交 2024-08，之后只有构建相关合并） | — |

2.15.2 release 说明：Docker / Apptainer 镜像升到 Ubuntu 24.04；bamsifter 编译修复；辅助脚本小改。
较早的相关变化：2.14.0 把 `--max_memory` 上限设为 200G（`Trinity` 脚本：">200G 时重置为 200G"）、最短 contig 长度绝对下限 100；
2.15.0 归一化改用 DB_File 存 read 名，大幅降低内存。**版本本身无需升级（已是最新）；但 conda build 有问题，见第 5、6 节。**

## 3. wgs_skim 中的调用方式

```bash
mem=$(( mem_mb / 1024 ))G     # mem_mb = max(16, clean_Gb × trinity_gib_per_gb) × 1.5^(attempt-1)
Trinity --seqType fq --left R1 --right R2 --max_memory $mem --CPU {threads} \
        --output assembly/<s>/trinity.Trinity_tmp --full_cleanup {trinity_extra} > logs/trinity/<s>.log 2>&1
mv assembly/<s>/trinity.Trinity_tmp.Trinity.fasta <s>.trinity.fasta
mv assembly/<s>/trinity.Trinity_tmp.Trinity.fasta.gene_trans_map <s>.trinity.fasta.gene_trans_map || true
```

| 参数 / 键 | 含义 |
|---|---|
| `--seqType fq --left/--right` | 双端 FASTQ（fastp 输出） |
| `--max_memory` | "suggested max memory"（Jellyfish、排序等可限制的部分）；>200G 会被重置为 200G |
| `--CPU` | 线程；规则 `threads: 16`。Inchworm 用 min(6, CPU)；Phase 2 由 ParaFly 用 `-CPU` 并行 |
| `--output …trinity.Trinity_tmp` | 目录名必须含 "trinity"（脚本安全检查）——满足 |
| `--full_cleanup` | 只保留 `${output}.Trinity.fasta`（及 gene_trans_map），删除整个工作目录 |
| `trinity_gib_per_gb`（6） | 内存估计；Wiki 经验值 "~1G RAM per ~1M ~76 base Illumina paired reads"，"often much less" |
| `trinity_extra`（""） | 追加选项，如 `--min_contig_length 300`、`--SS_lib_type RF`、`--no_version_check` |

## 4. 重要参数说明（Wiki Running-Trinity 与 2.15.2 `Trinity --show_full_usage_info`）

| 参数 | 默认 | 说明 |
|---|---|---|
| `--SS_lib_type RF/FR` | 无 | 链特异文库方向（dUTP = RF） |
| `--min_contig_length` | 200 | 最短输出 contig |
| `--normalize_max_read_cov` / `--no_normalize_reads` | 200 / 归一化开 | 归一化默认开启 |
| `--min_kmer_cov` | 1 | Inchworm 最低 k-mer 计数；提高可减内存 / 时间（minimize_RAM_usage） |
| `--workdir` | = `--output` | Phase 2 计算目录；可放本地 SSD 减少网络盘小文件 I/O |
| `--no_distributed_trinity_exec` | — | 只生成 Phase 2 命令列表后停止 |
| `--grid_exec "<cmd>"` | — | 用外部程序（HPC GridRunner、或任何"读一个命令文件、全部成功返回 0"的工具）分发 Phase 2 命令 |
| `--max_reads_per_graph` | 200000 | 每个图锚定的最大 reads 数 |
| `--FORCE` | — | 维护者在 #1220 建议：Phase 2 卡在最后一簇时用它收尾（该簇无结果） |
| `--no_salmon` | — | 不做最终 salmon 表达过滤 |
| `--no_version_check` | — | 不做联网版本检查（离线机器可加） |
| `--no_cleanup` | — | 与 `--full_cleanup` 互斥 |

## 5. 已知问题与 issues

| 现象 / 报错 | 原因 | 解决办法 | 链接 |
|---|---|---|---|
| **ParaFly `-CPU 16` 实际一次只跑一个命令，Phase 2 极慢** | 报告者："only one process is used to do the ParaFly Job"（bioconda 环境，2.15.1）；本机复现见第 6 节 | 换 Trinity 构建 / 用 `--grid_exec` 外部并行 / 官方 Docker·Singularity 镜像 | https://github.com/trinityrnaseq/trinityrnaseq/issues/1473 （open） |
| Phase 2 慢（一般情况） | 每簇一次 Trinity，受文件系统 / 机器负载影响 | 维护者：用 `--grid_exec` 分发到集群；Phase 1 只能提高 `--CPU` | https://github.com/trinityrnaseq/trinityrnaseq/issues/739 、 https://github.com/trinityrnaseq/trinityrnaseq/issues/1019 |
| Phase 2 停在 99.99 % 数天 | 某个难组装的簇（常为污染菌 / 内共生菌 / 质粒） | `--FORCE` 收尾，单独处理该簇 | https://github.com/trinityrnaseq/trinityrnaseq/issues/1220 |
| `commands in file: [FailedCommands] failed` | Phase 2 个别命令失败 | 维护者：原命令重跑，只重试失败的命令 | https://github.com/trinityrnaseq/trinityrnaseq/issues/1394 、 https://github.com/trinityrnaseq/trinityrnaseq/issues/1422 |
| **退出码 0 但没有 Trinity.fasta** | 输出目录残留旧的 `__salmon_filt.chkpts/salmon.filt.ok`，写最终结果的命令被跳过（2.13.2 实测，2.15.2 代码相同） | 重跑前删除 `__salmon_filt.chkpts` 和 `Trinity.tmp.fasta`（或整个输出目录） | https://github.com/trinityrnaseq/trinityrnaseq/issues/1477 |
| `tmp_normalized_read` 占满磁盘（66 GB 数据 → 666 GB） | 归一化中间文件 | 已过归一化：缩小大文件后续跑；或保留归一化 fastq 用 `--no_normalize_reads` 重新开始 | https://github.com/trinityrnaseq/trinityrnaseq/issues/1471 |
| 归一化阶段被杀（多样品） | 内存 | Wiki minimize_RAM_usage | https://github.com/trinityrnaseq/trinityrnaseq/issues/1472 |
| `jellyfish count: unrecognized option '--canonical'` | bioconda 环境中 jellyfish 版本不对 | 维护者只支持 Docker/Singularity（不维护 bioconda 包）；推断：conda 中确认 `jellyfish --version` 为 2.x | https://github.com/trinityrnaseq/trinityrnaseq/issues/1383 |
| `EOF marker is absent` / 读取 reads_partition 失败 | 中断后续跑的残留文件 / 非官方镜像 | 删除不完整的中间文件或换官方镜像 | https://github.com/trinityrnaseq/trinityrnaseq/issues/1323 、 https://github.com/trinityrnaseq/trinityrnaseq/issues/1325 |
| salmon 步骤报错 / 卡住 | salmon 安装或版本问题 | `--no_salmon` 或修复 salmon | https://github.com/trinityrnaseq/trinityrnaseq/issues/1444 、 https://github.com/trinityrnaseq/trinityrnaseq/issues/1365 |

rnaSPAdes vs Trinity：rnaSPAdes 论文（Bushmanova 2019）在小鼠等数据上 rnaSPAdes 组装出更多 50 %/95 % 完整基因，但作者也说
"not possible to detect the absolute leader according to all quality metrics"；独立比较 Hölzer & Marz 2019 *GigaScience* 8:giz039
（doi:10.1093/gigascience/giz039）认为 Trinity、SPAdes、Trans-ABySS 总体较好，且没有一个工具在所有数据上最好。
本项目 0.1 Gb 测试：rnaSPAdes ~90 s；Trinity 同数据 > 2 h（主要因下述 ParaFly 问题）。

## 6. 本项目中遇到的问题

- TESTING §9.2：两个约 0.1 Gb 的 RNA-seq 样品，Trinity Phase 2 有 13 961 个簇（A1），"在机器被其他 SPAdes 任务占满时很慢（20 分钟约 7 %）"。
- **2026-10-29 实测原因（本文档作者）**：运行中的 ParaFly 进程（`-CPU 16`）`/proc/<pid>/status` 显示 `Threads: 1`，任何时刻每个样品只有 1 个
  子 Trinity 进程；运行 1 h 47 min 后 A1 30 %、R1 22 %；工作目录 18 291 个文件、1.6 GB。用 8 条 `sleep 2` 测试：
  - wgs_skim 环境（`trinity-2.15.2-pl5321h077b44d_6`，libgomp 16.2.0）的 ParaFly `-CPU 8`：**16 s（串行）**；换成旧 libgomp 仍 16 s →
    问题在该 build 的 ParaFly 二进制本身；`OMP_NUM_THREADS`、`OMP_DYNAMIC` 无影响。
  - `~/miniforge3/envs/trinity`（build `pl5321h077b44d_3`）与 Trinity 2.8.5 的 ParaFly：**2 s（并行）**。
  - 结论：机器忙不是主因；当前 conda 解析出的 build _6 的 ParaFly 不并行。建议在 `trinity.yaml` 固定 `trinity =2.15.2=pl5321h077b44d_3`
    （本机验证并行），或在规则里用 `--grid_exec` 交给 GNU parallel 等工具。修改前需在干净环境复测。
- `--full_cleanup` 与续跑：Snakemake 重试时工作目录 `trinity.Trinity_tmp/` 仍在，Trinity 会按检查点续跑（`*.ok`、`recursive_trinity.cmds.completed`），
  这是好的；但若输入 reads 已变（fastp 重跑），旧检查点会被沿用，且有 #1477 的"空结果退出 0"风险——规则中 `mv` 会因文件不存在而失败，不会静默成功。
- 小文件：Phase 2 在 `read_partitions/Fb_*/CBin_*/` 下每簇写若干文件，`--full_cleanup` 在成功结束时全部删除；失败时需手动清理（记录到 CLEANUP 日志）。

- **2026-10-30 salmon 2.x**：固定 build _3 后重建的环境解析出 salmon 2.8.0（Rust 重写版）。Trinity 2.15.2 的
  `salmon_runner.pl` 调用 `salmon quant --minAssignedFrags 1 --validateMappings`，2.x 已删除这些选项：
  `error: unexpected argument '--minAssignedFrags' found`，Trinity 在表达量过滤一步退出（tools2 trinity-busco：A1、R1 都失败）。
  处理：`trinity.yaml` 固定 `salmon =1.10.3`（最后的 C++ 版）。

## 7. 排错清单

1. `logs/trinity/<s>.log`：阶段标记 `Trinity Phase 1`、`Inchworm`、`Chrysalis`、`Trinity Phase 2: Assembling Clusters of Reads`、
   ParaFly 进度 `succeeded(N) X% completed`；结尾应有 `All commands completed successfully` 与 `Trinity.fasta`。
2. 关键字：`Error`、`FailedCommands`、`died with ret`、`Killed`、`std::bad_alloc`、`unrecognized option`、`No space left`、`EOF marker`。
3. Phase 2 是否真的并行：`ps -eo args | grep -c "trinity.reads.fa --output"`（应接近 `--CPU`）；`grep Threads /proc/$(pgrep -f "ParaFly -c recursive")/status`。
4. 进度：`wc -l trinity.Trinity_tmp/recursive_trinity.cmds` 与 `wc -l recursive_trinity.cmds.completed`。
5. 磁盘 / inode：`du -sh`、`find trinity.Trinity_tmp | wc -l`、`df -i`。
6. 重跑前确认工作目录里没有不属于本次输入的旧检查点（#1477）。

## 8. 资料来源与查阅日期（2026-10-29）

- Releases：https://github.com/trinityrnaseq/trinityrnaseq/releases ；Wiki Running-Trinity；2.15.2 `Trinity` 脚本（本机 conda pkgs）。
- 上表 issues 及维护者回复（`gh api`，2026-10-29）。
- 本机实测：`wgs_skim_test/rna2_out`（运行中）、ParaFly 计时测试（scratchpad）。
- wgs_skim `workflow/rules/transcriptome.smk`、`workflow/envs/trinity.yaml`、`config/config.yaml`、`docs/TESTING.md`。

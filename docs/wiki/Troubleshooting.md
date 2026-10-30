# 故障排查

## 先看哪里

| 文件 | 内容 |
|---|---|
| `pipeline_info/logs/<任务>/<样品>.log` | 该任务所调用工具自己的输出（错误信息通常在最后几行） |
| `pipeline_info/benchmarks/<任务>/<样品>.tsv` | 时间（s）与峰值内存（max_rss，MB）—— 判断是否内存不足 |
| `.snakemake/log/*.log` | 每次运行的 Snakemake 日志：哪些任务失败（`Error in rule …`） |
| `wgs_skim run ... -n -p` | 不运行，列出每个任务的完整 shell 命令（可复制出来单独运行调试） |
| `<marker>/getorganelle/<s>/rounds.tsv` | rDNA / 组蛋白每轮的数据量、结果、基因覆盖、是否采用 |
| `mitogenome/<s>/<s>.mt.json` | 线粒体来源、成环与否、迭代停止原因 |

一个任务失败，其他样品照常进行（`keep-going`），最后的汇总等待；修正问题后运行**同一条命令**即可续跑。

## 常见问题

| 信息 / 现象 | 原因 | 处理 |
|---|---|---|
| `no paired FASTQ files ... found` | 文件名无法识别 | `wgs_skim samples -r DIR` 看识别结果；否则用样品表 `-s`（见 [输入](Input.md)） |
| `sample 'X' found twice` | 同一样品两个 lane / 两种命名 | 合并 lane，或用样品表 |
| `Directory cannot be locked` / `IncompleteFilesException` | 上次运行被中断 | `wgs_skim run ... -- --unlock`，再运行同一条命令（不完整的文件会重做） |
| `CreateCondaEnvironmentException`，`wgs_skim install` 中止 | 建环境时网络 / 频道问题 | 再运行 `wgs_skim install`（已建好的环境保留）；`--conda-prefix` 放到剩余 > 20 GB 的磁盘 |
| `unknown step 'x'` / `unknown SPAdes mode 'x'` | 拼写错误 | 信息中列出全部可用值 |
| `paths must contain only ASCII characters` | 结果目录或输入路径含中文等非 ASCII 字符；SPAdes 会拒绝（`path contains non-ASCII characters`），wgs_skim 在开始前就检查 | 用英文 / 拼音命名目录，如 `results_2026` |
| Trinity：`unexpected argument '--minAssignedFrags'` | 环境里是 salmon 2.x，Trinity 2.15.2 需要 1.x | 已固定 `salmon =1.10.3`；旧环境删掉后重建（`wgs_skim install`） |
| MitoFinder：`MitoFinder is not installed ... Aborting` | bioconda 1.4.1 环境缺 `install.sh.ok` | 用默认的 1.4.2 镜像（`mitofinder_sif`）；wgs_skim 把它报为失败，不当作"没找到" |
| GetOrganelle：`Too few seed reads found` | 种子与样品相距太远（其他科），或 reads 太少 | 组蛋白已用局部比对招募；自己给种子（`histone_seed`、`nr_seed`），或 `nuclear_max_gb: 10,20` |
| GetOrganelle：`No space left on device` | 它把用到的 reads 解压到输出目录（最多数据量 1.5 倍） | 结果目录放到大磁盘；`nuclear_max_gb`、`mt_iter_max_data_gb` 限制了数据量 |
| 线粒体没有成环（`mt_circular: no`） | 重复迭代仍未闭合：控制区重复、覆盖度低 | `topology/` 说明末端能否连接；`mt_iter_max_data_gb: 0` 用全部 reads；用 `wgs_skim compare` 对照 MitoFinder 结果 |
| SPAdes：`<jemalloc>: ... out of memory` 或被杀 | 内存估计（`spades_gib_per_gb` × clean 数据）偏低 | 自动重试给 1.5 倍内存并从检查点继续；或提高 `spades_min_mem_gb` |
| SPAdes 没有运行（`clean data ... < spades_min_clean_gb`） | 数据量低于阈值 | 降低 `spades_min_clean_gb`（0 = 全部组装） |
| Trinity 很慢（第二阶段 `succeeded(n) x% completed` 增长缓慢） | 第二阶段对每个 read 簇运行一次 Trinity（小数据也有上万个），大量小文件，磁盘 / CPU 被其他任务占用时很慢 | 机器空闲时运行；给 Trinity 更多线程 `-- --set-threads trinity=32` |
| `MitoFinder needs a reference mitogenome GenBank` | MitoFinder 必须有参考 | `--ref related.gb`（`wgs_skim mitofinder`）或 `--mitofinder-ref`（`run`）；与 `meangs-getorganelle` 一起运行时自动用样品自己的 MitoM 注释 |
| `MitoFinder found no mitochondrial contig`（结果为空文件） | 组装中没有与参考相似的 contig：参考太远、或线粒体覆盖度太低没有组装出来 | 换更近的参考；看 `pipeline_info/logs/mitofinder_<来源>/<s>.log`；用 GetOrganelle 组装 |
| MitoFinder：`singularity: command not found` / 镜像不存在 | 配置了 `mitofinder_sif` 但本机没有 Singularity 或文件 | 修正路径，或 `--config mitofinder_sif=` 改用 bioconda MitoFinder 1.4.1 |
| `no FASTA files ... in` | `-i` 目录中没有 .fasta / .fa / .fna / .fas | 检查扩展名；或直接列出文件 |
| `sample 'X' found twice`（`-i`） | 两个文件的名称（去掉扩展名后）相同 | 改名，或分两次运行 |
| BUSCO：`lineage ... not found` / 下载失败 | 无网络或名称错误 | `--busco-lineage /path/to/metazoa_odb12`（已下载的谱系目录） |
| BLAST 要几个小时 | 每次搜索都读一遍完整 nt | 一次性 `wgs_skim build-blastdb`（后生动物 ≤ 50 kb：68 GB，每次 1–2 分钟），然后 `--blast-db` |
| `blast_topn.tsv` 中 `(efetch failed)` | NCBI Entrez 未返回这些命中的元数据 | 比对结果仍在；以后 `-- --forcerun blast` 重做 |
| `mitom_db`：`git clone` 失败 | 无网络，无法获取固定版本的 MitoM / MITOS2 参考 | 配置 `mitom_dir: ~/mybiosofts/Mitom` 和 `mitos_refdir:`（不再需要下载） |
| 后台运行消失 | 机器内存不足，进程被系统（OOM killer / 会话管理）停止 | 用 `nohup … &` 或 `tmux` 运行；`--mem-gb` 让本流程只占它的份额 |
| 报告中 `rdna_copies` / `histones` 为空或 `x?` | 没有 SPAdes 组装或 BUSCO 单拷贝基因可作标准 | 需要 `spades` + `busco`（`nuclear` 步骤会运行它们）；数据太少时无法估计 |
| `wgs_skim compare --results` 中样品被 `skipped` | 该样品缺少 GetOrganelle 或 MitoFinder 结果 | `all.compare.tsv` 末尾列出缺少的文件 |
| `wgs_skim nuclear`：`no known results layout` | 结果目录既不是本软件（`qc/`、`assembly/`）也不是旧流程（`01.fastp/`、`03.spades/`）的结构 | 指向正确的结果目录；或对原始数据用 `wgs_skim run --steps nr,nuclear` |

## 仍然解决不了

1. `wgs_skim run ... -n -p` 找到失败任务的完整命令；
2. 激活该任务的 conda 环境（日志开头 `Activating conda environment: …`）单独运行它；
3. 在 GitHub 提 issue，附上 `pipeline_info/logs/<任务>/<样品>.log` 和 `benchmarks/…tsv`。

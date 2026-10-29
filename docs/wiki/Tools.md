# 单个工具

每个工具都可以**单独、批量**运行：一条命令处理目录中的全部样品，由 Snakemake 按核数和内存调度，每个样品一个日志，
中断后运行同一条命令即可续跑。所有工具共用这些选项：

| 选项 | 含义 |
|---|---|
| `-r DIR` / `-s sheet.tsv` | reads 输入（reads 类工具）：目录（样品名自动识别）或样品表，见 [输入数据](Input.md) |
| `-i FILE\|DIR …` | FASTA 输入（BUSCO、MitoFinder、MitoM、BLAST）：文件或目录，样品名 = 文件名；`--split-by-id` 拆分多样品合并的 FASTA |
| `-o DIR` | 结果目录（默认 `./<工具>_results`） |
| `-j N`、`--mem-gb N` | 全部任务合计最多用的核数、内存（默认全部核、本机内存 90%） |
| `-t TAXON` | 类群预设：annelid（默认）/ mollusc / fish / crustacean —— 遗传密码、MitoM 类群、BUSCO 谱系、MEANGS 类群 |
| `--qc` | reads 类工具：先用 fastp 清洗（默认直接用给出的 reads，视为已清洗） |
| `-n` | 只显示将运行的任务 |
| `--config 键=值` | 任何配置项，见 [配置参数](Configuration.md) |
| `-- …` | 之后的内容原样传给 Snakemake，如 `-- --set-threads spades=48` |

| 工具 | 命令 | 输入 | 做什么 |
|---|---|---|---|
| [fastp](Tool-fastp.md) | `wgs_skim fastp` | reads | 质控、去接头 |
| [MEANGS](Tool-MEANGS.md) | `wgs_skim meangs` | reads | 线粒体种子序列 |
| [GetOrganelle](Tool-GetOrganelle.md) | `wgs_skim getorganelle` | reads | 线粒体（未成环迭代）、rDNA、组蛋白 |
| [SPAdes](Tool-SPAdes.md) | `wgs_skim spades` | reads | 组装：genome / rna / meta / plasmid … 全部模式 |
| [Trinity](Tool-Trinity.md) | `wgs_skim trinity` | reads | 转录组组装 |
| [BUSCO](Tool-BUSCO.md) | `wgs_skim busco` | FASTA | 组装 / 转录组 / 蛋白完整度 |
| [MitoFinder](Tool-MitoFinder.md) | `wgs_skim mitofinder` | FASTA 或 reads | 从组装中找出并注释线粒体 |
| [MitoM](Tool-MitoM.md) | `wgs_skim mitom` | FASTA | 线粒体注释（GenBank，NCBI 检查） |
| [BLAST](Tool-BLAST.md) | `wgs_skim blast` | FASTA | 比对 nt，Top 命中 + GenBank 物种分类 |
| [GenomeScope](Tool-GenomeScope.md) | `wgs_skim genome-size` | reads | 基因组大小、杂合度 |

每个工具的版本、官方参数、已知问题（GitHub issues）与排错方法：仓库 [docs/tools/](../tools/)。

多个工具串起来的推荐用法见 [组合流程](Combinations.md)。

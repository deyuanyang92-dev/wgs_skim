# wgs_skim 使用手册

wgs_skim 是动物基因组浅测序（genome skimming）与转录组数据的批量分析软件。手册分三部分：

1. **[单个工具](Tools.md)**：fastp、MEANGS、GetOrganelle、SPAdes、Trinity、BUSCO、MitoFinder、MitoM、BLAST、GenomeScope ——
   每个工具一条命令（`wgs_skim fastp`、`wgs_skim spades`、`wgs_skim busco` …），对目录中全部样品批量运行；
   输入可以是原始 / clean reads，也可以是已有的组装或序列。
2. **[组合流程](Combinations.md)**：推荐的工具组合，一条命令从原始数据到结果：
   - 线粒体：**MEANGS-GetOrganelle**（未成环自动迭代），或 SPAdes-MitoFinder；
   - 全基因组：**SPAdes-BUSCO**；
   - 转录组：**rnaSPAdes-BUSCO**（推荐）或 Trinity-BUSCO；
   - 核 rDNA / 组蛋白、物种鉴定（BLAST）、全流程。
3. **参考**：资源调度与续跑、全部配置项、输出文件、报错排查、测试记录、旧脚本对照，以及每个软件的版本、参数和已知问题
   （[docs/tools/](../tools/)）。

所有批量任务都由 **Snakemake** 管理：按 CPU 与内存调度、每个任务单独日志、中断后运行同一条命令断点续跑；
也可以用 **Nextflow** 运行同一流程（`wgs_skim run --engine nextflow`，见 [资源调度与续跑](Resources-and-resume.md)）。

```bash
wgs_skim busco -i assemblies/ -o busco_out/                       # 一个工具
wgs_skim run -r raw/ -o res/ --steps spades-busco                 # 一个组合
wgs_skim run -r raw/ -o res/                                      # 全部
```

源码：<https://github.com/deyuanyang92-dev/wgs_skim> · 作者：De-Yuan Yang（YZZ Team）

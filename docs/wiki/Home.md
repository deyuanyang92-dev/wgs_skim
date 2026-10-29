# wgs_skim 使用手册

wgs_skim 是动物基因组浅测序（genome skimming）的批量分析软件：给一个原始数据目录，自动完成质控、线粒体组装与注释、
核多拷贝基因（18S-ITS-28S rDNA、组蛋白 H1–H4）组装与注释、SPAdes 全基因组组装（含 rnaSPAdes、meta、plasmid 等各种模式）
与 Trinity 转录组组装、BUSCO 评估、基因组大小、BLAST 物种鉴定，最后生成汇总表与报告。

所有批量任务都由 **Snakemake** 管理（按 CPU 与内存调度、每个任务单独日志、中断后断点续跑），也可以用 **Nextflow** 运行同一流程（`--engine nextflow`，结果与 Snakemake 完全相同，见 [资源调度与续跑](Resources-and-resume.md)）。

| 页面 | 内容 |
|---|---|
| [安装](Installation.md) | 安装软件与数据库（一次） |
| [快速开始](Quick-start.md) | 三条命令跑完一个项目 |
| [输入数据](Input.md) | 文件命名、样品表、已有组装 |
| [步骤详解](Steps.md) | 每个步骤：做什么、命令、输出、实测时间与内存 |
| [线粒体](Mitogenome.md) | MEANGS → GetOrganelle → 迭代 → SPAdes 兜底、MitoM 注释、环 / 线性验证、无 reads 环化、组装比较 |
| [核多拷贝基因](Nuclear-markers.md) | rDNA 与组蛋白：GetOrganelle 组装、MitoM 注释（.gbf）、与 SPAdes 比较 |
| [SPAdes 各模式与转录组](Assembly-modes.md) | genome / rna / meta / plasmid … 、Trinity、BUSCO |
| [BLAST 物种鉴定](BLAST.md) | 本地小库、在线 BLAST、查询基因 |
| [资源调度与续跑](Resources-and-resume.md) | Snakemake：CPU 与内存、进度、续跑、集群 |
| [配置参数](Configuration.md) | 全部配置项 |
| [输出文件](Output.md) | 结果目录结构、汇总表各列 |
| [报错排查](Troubleshooting.md) | 常见报错、原因、解决办法 |
| [测试](Testing.md) | 测了什么、怎么测、结果 |
| [旧脚本对照](From-old-scripts.md) | Temp_scripts 中各脚本的功能在哪里 |

源码：<https://github.com/deyuanyang92-dev/wgs_skim> · 作者：De-Yuan Yang（YZZ Team）

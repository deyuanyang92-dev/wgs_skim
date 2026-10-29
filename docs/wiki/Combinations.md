# 组合流程

单个工具见 [单个工具](Tools.md)。这里是把工具串起来的推荐用法：每个组合一条命令，从原始 reads 开始（先 fastp），
中间结果自动传给下一个工具，已完成的部分不重做。

| 目标 | 推荐组合 | 命令 |
|---|---|---|
| 线粒体基因组 | **MEANGS-GetOrganelle**（未成环自动迭代） | `wgs_skim run -r raw/ -o res/ --steps meangs-getorganelle` |
| 线粒体（另一种方法 / 对照） | SPAdes-MitoFinder | `wgs_skim run -r raw/ -o res/ --steps spades-mitofinder --mitofinder-ref ref.gb` |
| 全基因组（二代数据） | **SPAdes-BUSCO** | `wgs_skim run -r raw/ -o res/ --steps spades-busco` |
| 转录组 | **rnaSPAdes-BUSCO**（推荐）、Trinity-BUSCO | `--steps rnaspades-busco` / `--steps trinity-busco` / 两者 `--steps rna` |
| 转录组中的线粒体基因 | rnaSPAdes-MitoFinder、Trinity-MitoFinder（多样品按基因汇总） | `--steps rnaspades-mitofinder` / `--steps trinity-mitofinder` |
| 核 rDNA、组蛋白 | GetOrganelle-MitoM，与 SPAdes-MitoM 互相验证 | `--steps nr` / `--steps nuclear` |
| 物种鉴定 | MitoM / mitom nuc → BLAST | `--steps blast` |
| 全部 | 以上全部 + 基因组大小 + 汇总报告 | `wgs_skim run -r raw/ -o res/` |

组合可以一起写：`--steps meangs-getorganelle,spades-mitofinder,spades-busco`。

| 页面 | 内容 |
|---|---|
| [线粒体](Combo-mitogenome.md) | MEANGS-GetOrganelle + 迭代、SPAdes-MitoFinder、两者比较、环 / 线性验证、无 reads 环化、注释 |
| [全基因组](Combo-genome.md) | SPAdes-BUSCO、基因组大小 |
| [转录组](Combo-transcriptome.md) | rnaSPAdes-BUSCO、Trinity-BUSCO；转录组中的线粒体基因（rnaSPAdes / Trinity-MitoFinder） |
| [核多拷贝基因](Combo-nuclear.md) | rDNA 与组蛋白：GetOrganelle、MitoM .gbf、与 SPAdes 比较 |
| [物种鉴定](Combo-species.md) | COX1、18S、ITS、28S 的 BLAST |
| [全流程](Combo-all.md) | 一条命令完成全部分析，汇总表与报告 |

# 转录组：rnaSPAdes-BUSCO / Trinity-BUSCO

两种组装方法，都在组装后用 BUSCO（transcriptome 模式）评估。**推荐 rnaSPAdes-BUSCO**：速度快、内存少；Trinity 可作对照。

```bash
wgs_skim run -r rnaseq/ -o res/ --steps rnaspades-busco     # 推荐：fastp → rnaSPAdes → BUSCO
wgs_skim run -r rnaseq/ -o res/ --steps trinity-busco       # fastp → Trinity → BUSCO
wgs_skim run -r rnaseq/ -o res/ --steps rna                 # 两者都做，各自 BUSCO
```

| 步骤 | 工具 | 输出 `spades/<s>/` |
|---|---|---|
| 质控 | [fastp](Tool-fastp.md) | — |
| 组装 | [SPAdes](Tool-SPAdes.md) `--rna -k 21,33,55` | `<s>.rna.fasta` |
| 组装 | [Trinity](Tool-Trinity.md) `--full_cleanup` | `<s>.trinity.fasta` |
| 评估 | [BUSCO](Tool-BUSCO.md) `-m transcriptome` | `<s>.rna.busco.txt`、`<s>.trinity.busco.txt` |

实测（2 个样品，各约 0.1 Gb RNA-seq）：

| 样品 | rnaSPAdes 转录本 | 时间 / 内存 | BUSCO（metazoa） |
|---|---|---|---|
| A1 | 20 573 | 87 s / 1.3 GB | C:12.1% |
| R1 | 27 620 | 92 s / 1.2 GB | C:13.2% |

BUSCO 低是数据量太小（完整转录组通常 ≥ 2 Gb）；Trinity 同样数据第二阶段在繁忙机器上需要数小时。

已有转录本只想评估：`wgs_skim busco -i transcripts/ -m transcriptome`。

## 从转录组获取线粒体基因：rnaSPAdes-MitoFinder / Trinity-MitoFinder

只要线粒体基因（例如多个物种的 COX1、CYTB、ND 基因用于建树）时，不必组装线粒体全长：转录组中线粒体基因的表达量很高，
转录本里就有。MitoFinder 在转录本中找出并注释它们，再按基因汇总所有样品。

```bash
wgs_skim run -r rnaseq/ -o res/ --steps rnaspades-mitofinder --mitofinder-ref related.gb
wgs_skim run -r rnaseq/ -o res/ --steps trinity-mitofinder --mitofinder-ref related.gb
wgs_skim run -r rnaseq/ -o res/ --steps rnaspades-mitofinder,trinity-mitofinder,rnaspades-busco --mitofinder-ref related.gb
wgs_skim mitofinder -i transcripts/ -o mf/ --ref related.gb        # 已有转录本
```

流程：fastp → rnaSPAdes / Trinity → MitoFinder（每个样品）→ 按基因汇总所有样品 → **BLAST 污染检查**（每个基因比对 nt，
最佳命中是别的科且高度相似时标记）。

结果在 `mitofinder_rnaspades/`（rnaSPAdes）和 `mitofinder_trinity/`（Trinity）；污染检查在 `mitofinder_<来源>/blast/samples_check.tsv`；`mitofinder_summary.tsv`（每个样品每个基因的长度）、
`genes/<基因>.fasta`（所有样品的同一基因，可直接比对），详见 [MitoFinder](Tool-MitoFinder.md)。

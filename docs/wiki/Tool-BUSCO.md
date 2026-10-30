# BUSCO：完整度评估

对任意组装、转录组或蛋白序列评估单拷贝直系同源基因的完整度。

```bash
wgs_skim busco -i assemblies/ -o busco_out/                            # 目录中每个 FASTA = 一个样品
wgs_skim busco -i A.fasta B.fasta -o busco_out/
wgs_skim busco -i transcripts/ -o busco_out/ -m transcriptome
wgs_skim busco -i proteins/ -o busco_out/ -m proteins -l arthropoda_odb12
```

- `-m genome | transcriptome | proteins`；`-l` 谱系名称（第一次自动下载）或已下载的谱系目录；默认按 `--taxon`（环节动物 metazoa_odb12）。
- BUSCO 6.0.0，`--offline`；每个任务 16 线程、16 GB。
- 输出：`busco/<s>/<s>.busco.txt`（short summary：C/S/D/F/M）、`spades/<s>/busco/`（完整结果）。
  输入文件以链接的形式放在 `spades/<s>/<s>.scaffolds.fasta`。
- 实测：44 Mb 组装 40 s / 1.7 GB；0.1 Gb RNA-seq 的转录本 1–2 min。

在组合流程中 BUSCO 自动评估每个组装：`spades-busco`、`rnaspades-busco`、`trinity-busco`。

深入资料：[docs/tools/BUSCO.md](../tools/BUSCO.md)。

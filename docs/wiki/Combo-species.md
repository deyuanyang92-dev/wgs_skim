# 物种鉴定：BLAST

每个样品的 COX1（来自 MitoM 注释）和 18S、ITS、28S（来自 rDNA 注释）比对 nt，取 Top N 命中，附 GenBank 物种与分类。

```bash
wgs_skim run -r raw/ -o res/ --steps blast                                    # 组装 + 注释 + BLAST
wgs_skim run -r raw/ -o res/ --steps blast --blast-mt cox1,16s,cytb --blast-nr 18S,ITS,28S
wgs_skim run -r raw/ -o res/ --steps blast --blast-mt full --blast-nr none     # 只比整条线粒体
wgs_skim blast -i seqs/ -o blast_out/                                         # 已有序列：单独的 BLAST 工具
```

查询：`--blast-mt` = `full` 或基因（cox1 = COI、cox2、cox3、cytb、nad1–6、nad4l、atp6、atp8、16s、12s）；`--blast-nr` = `full`、
`18S`、`ITS`（ITS1-5.8S-ITS2）、`ITS1`、`5.8S`、`ITS2`、`28S`；`none` 关闭。没有该基因的样品用整条序列代替（日志注明）。

数据库与速度见 [BLAST](Tool-BLAST.md)。输出 `blast/blast_topn.tsv`；汇总表 `blast_mt` / `blast_nr` 列为每个查询的最佳命中。

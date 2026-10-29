# BLAST 物种鉴定

每个样品的 COX1（来自 MitoM 注释）和 18S、ITS、28S（来自 MitoM nuc 注释）比对 nt，取 Top N 个命中，附 GenBank 的物种名与分类。

```bash
wgs_skim run -r raw/ -o res/ --steps blast --blast-db /db/blast_marker/metazoa_marker   # 本地库（快）
wgs_skim run -r raw/ -o res/ --steps blast --email you@example.com                      # NCBI 在线（慢）
wgs_skim run -r raw/ -o res/ --steps blast --blast-mt cox1,16s,cytb --blast-nr 18S,ITS,28S
wgs_skim run -r raw/ -o res/ --steps blast --blast-mt full                              # 整条线粒体
```

查询：`--blast-mt` 可用 `full` 或基因（cox1 = COI、cox2、cox3、cytb、nad1–6、nad4l、atp6、atp8、16s、12s）；
`--blast-nr` 可用 `full`、`18S`、`ITS`（ITS1-5.8S-ITS2）、`ITS1`、`5.8S`、`ITS2`、`28S`。没有该基因的样品用整条序列代替（日志中注明）。

## 快速本地库：`wgs_skim build-blastdb`

全 nt（约 1.2 TB）每次搜索都要读一遍，实测一次 COX1 + rDNA 搜索 42–52 分钟。条形码 / 标记基因只会命中短记录，
所以从本地 nt 提取后生动物 ≤ 50 kb 的记录建小库：68 GB，每次搜索 1–2 分钟，命中与全 nt 相同。

```bash
wgs_skim build-blastdb --nt /db/blast/nt --out /db/blast_marker/metazoa_marker            # 一次，数小时
wgs_skim build-blastdb --nt /db/blast/nt --out /db/x --taxid 6340 --max-length 30000      # 只要环节动物
```

把 `blast_db: /db/blast_marker/metazoa_marker` 写进 `~/.wgs_skim/config.yaml`，以后每次自动使用。

输出 `blast/blast_topn.tsv`：查询、命中、一致度、覆盖度、物种、分类、标题、长度；NCBI 没有返回元数据的命中写
`(efetch failed)`（比对结果本身保留）。汇总表 `blast_mt` / `blast_nr` 列为每个查询的最佳命中。

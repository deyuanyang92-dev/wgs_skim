# BLAST：物种鉴定

把序列比对到 nt，每个查询取 Top N 个命中，附 GenBank 的物种名与分类。

```bash
wgs_skim blast -i seqs/ -o blast_out/                                   # 整条序列（默认 --genes full）
wgs_skim blast -i mitogenomes/ -o blast_out/ --genes cox1,16s           # 先 MitoM 注释，再取基因比对
wgs_skim blast -i seqs/ -o blast_out/ --db /db/blast_marker/metazoa_marker --topn 10
```

- 数据库：`--db` 本地库（快）；不给则 NCBI 在线 BLAST（慢，无需下载）。个人配置 `blast_db` 设置后每次自动使用。
- 快速本地库：`wgs_skim build-blastdb --nt /db/nt --out /db/metazoa_marker`（后生动物 ≤ 50 kb，68 GB，每次搜索 1–2 min，
  命中与完整 nt 相同；完整 nt 每次 40–50 min）。
- 输出：`blast/blast_topn.tsv`：查询、命中、一致度、覆盖度、物种、分类、标题、长度；NCBI 未返回元数据时写 `(efetch failed)`。
- 内存：本地库被映射进内存，实测峰值 47.5 GB（任务申请 48 GB）。

在组合流程中自动使用 COX1（MitoM 注释）和 18S、ITS、28S（rDNA 注释）：`--blast-mt`、`--blast-nr`，见 [物种鉴定](Combo-species.md)。

深入资料：[docs/tools/BLAST.md](../tools/BLAST.md)。

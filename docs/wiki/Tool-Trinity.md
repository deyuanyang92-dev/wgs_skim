# Trinity：转录组组装

```bash
wgs_skim trinity -r clean_rnaseq/ -o trinity_out/
wgs_skim trinity -r rnaseq/ -o trinity_out/ --qc --config trinity_extra="--min_contig_length 300"
```

- Trinity 2.15.2，`--seqType fq --full_cleanup`，`--max_memory` 按数据量估计（`trinity_gib_per_gb` = 6 GiB/Gb，至少 16 GB），16 线程。
- 输出：`assembly/<s>/<s>.trinity.fasta`、`<s>.trinity.fasta.gene_trans_map`。
- 注意：第二阶段为每个 read 簇单独运行一次 Trinity（小数据也有上万个簇），机器繁忙时很慢；可给更多线程：
  `-- --set-threads trinity=32`。

推荐组合：[Trinity-BUSCO](Combo-transcriptome.md)。深入资料：[docs/tools/Trinity.md](../tools/Trinity.md)。

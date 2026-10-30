# SPAdes：组装（全部模式）

```bash
wgs_skim spades -r clean_reads/ -o spades_out/                          # 全基因组（genome）
wgs_skim spades -r clean_reads/ -o spades_out/ --mode rna               # 转录组 rnaSPAdes
wgs_skim spades -r clean_reads/ -o spades_out/ --mode meta,plasmid      # 多个模式：各自一个组装
wgs_skim spades -r raw/ -o spades_out/ --qc --mode isolate --config spades_careful=true
```

| `--mode` | SPAdes 参数 | 输出 `spades/<s>/` |
|---|---|---|
| `genome`（默认） | `--only-assembler` | `<s>.scaffolds.fasta` |
| `isolate` / `sc` / `meta` | `--isolate` / `--sc` / `--meta` | `<s>.<mode>.fasta`（scaffolds） |
| `rna` / `rnaviral` | `--rna` / `--rnaviral`，`-k 21,33,55` | `<s>.rna.fasta`（transcripts） |
| `plasmid` / `metaplasmid` / `metaviral` | 同名参数 | `<s>.<mode>.fasta` |
| `bio` | `--bio` | `<s>.bio.fasta`（gene_clusters） |
| `sewage` / `corona` | 同名参数 | `<s>.<mode>.fasta` |

- 内存：按 clean 数据量估计（genome：`spades_gib_per_gb` = 6.5 GiB/Gb，至少 32 GB；其他模式 5 GiB/Gb，至少 16 GB），
  同时作为 SPAdes 的 `-m`；内存不足自动重试一次（1.5 倍），从最后的检查点继续（`--restart-from last`）。
- 其他配置：`spades_kmers`、`rnaspades_kmers`、`spades_careful`、`spades_only_assembler`、`spades_min_clean_gb`。
- 中间文件（K*/、tmp/、misc/、corrected/）在组装完成后删除。
- 实测（0.3 Gb）：genome 6–8 min / 2.8 GB；rna（0.1 Gb RNA-seq）90 s / 1.3 GB。

推荐组合：全基因组 [SPAdes-BUSCO](Combo-genome.md)；转录组 [rnaSPAdes-BUSCO](Combo-transcriptome.md)；
线粒体 [SPAdes-MitoFinder](Combo-mitogenome.md)。

深入资料：[docs/tools/SPAdes.md](../tools/SPAdes.md)。

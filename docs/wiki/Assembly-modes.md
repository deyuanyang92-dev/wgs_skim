# SPAdes 各模式与转录组

## SPAdes：`--spades-modes`

与旧脚本 `spadesv3.5.5.py` 相同，支持 SPAdes 的全部模式，可同时跑多个，每个模式每个样品一个任务，各自 BUSCO 评估：

| 模式 | SPAdes 参数 | 结果文件 | BUSCO 模式 |
|---|---|---|---|
| `genome`（默认） | （`--only-assembler`） | `assembly/<s>/<s>.scaffolds.fasta` | genome |
| `isolate` | `--isolate` | `<s>.isolate.fasta` | genome |
| `sc` | `--sc`（单细胞 MDA） | `<s>.sc.fasta` | genome |
| `meta` | `--meta` | `<s>.meta.fasta` | genome |
| `rna` | `--rna`（rnaSPAdes） | `<s>.rna.fasta`（transcripts） | transcriptome |
| `rnaviral` | `--rnaviral` | `<s>.rnaviral.fasta` | transcriptome |
| `plasmid` / `metaplasmid` | `--plasmid` / `--metaplasmid` | `<s>.plasmid.fasta` … | genome |
| `metaviral` | `--metaviral` | `<s>.metaviral.fasta` | genome |
| `bio` | `--bio` | `<s>.bio.fasta`（gene_clusters） | genome |
| `sewage` / `corona` | `--sewage` / `--corona` | … | genome |

`genome` 模式的结果供 BUSCO、核基因、线粒体兜底使用；其他模式是独立的组装。

```bash
wgs_skim run -r raw/ -o res/ --steps spades,busco --spades-modes genome,meta,plasmid
wgs_skim run -r raw/ -o res/ --steps spades --spades-modes isolate --config spades_careful=true spades_kmers=21,33,55
```

相关配置：`spades_kmers`（非 genome 模式的 -k）、`rnaspades_kmers`（rna / rnaviral，默认 21,33,55）、`spades_careful`
（isolate / sc / plasmid）、`spades_only_assembler`、内存估计 `spades_gib_per_gb`（genome）/ `spades_mode_gib_per_gb`（其他）。
中断后重跑会从 SPAdes 最后的检查点继续（`--restart-from last`）。

## 转录组：rnaSPAdes 与 Trinity

```bash
wgs_skim run -r rnaseq/ -o res/ --steps rna        # rnaSPAdes + Trinity + 两者 BUSCO（-m transcriptome）
wgs_skim run -r rnaseq/ -o res/ --steps trinity    # 只跑 Trinity
wgs_skim run -r rnaseq/ -o res/ --steps rna --config trinity_extra="--min_contig_length 300"
```

- Trinity 2.15.2（与旧脚本 `trinity_batch3.1.py` 使用的版本相同），`--full_cleanup`，`--max_memory` 按数据量估计（`trinity_gib_per_gb`）。
- 输出：`assembly/<s>/<s>.rna.fasta`、`<s>.trinity.fasta`（+ `.gene_trans_map`）、`<s>.rna.busco.txt`、`<s>.trinity.busco.txt`。
- BUSCO 谱系由 `--taxon` 预设或 `--busco-lineage` 指定。

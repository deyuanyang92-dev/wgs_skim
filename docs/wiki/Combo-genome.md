# 全基因组：SPAdes-BUSCO

二代数据的全基因组组装推荐 SPAdes，组装后用 BUSCO 评估完整度。

```bash
wgs_skim run -r raw/ -o res/ --steps spades-busco                        # fastp → SPAdes → BUSCO
wgs_skim run -r raw/ -o res/ --steps spades-busco -j 96 --mem-gb 600 -- --set-threads spades=48
wgs_skim run -r raw/ -o res/ --steps spades-busco,genome_size            # + 基因组大小、杂合度
```

已有组装、只想评估：`wgs_skim busco -i assemblies/`；已有 clean reads、只想组装：`wgs_skim spades -r clean/`。

| 步骤 | 工具 | 输出 |
|---|---|---|
| 质控 | [fastp](Tool-fastp.md) | `qc/<s>/` |
| 组装 | [SPAdes](Tool-SPAdes.md) genome 模式，`--only-assembler` | `assembly/<s>/<s>.scaffolds.fasta` |
| 评估 | [BUSCO](Tool-BUSCO.md) genome 模式 | `assembly/<s>/<s>.busco.txt` |
| 基因组大小（可选） | [Jellyfish + GenomeScope2](Tool-GenomeScope.md) | `genome_size/<s>/` |

内存：SPAdes 按数据量估计（6.5 GiB / Gb clean 数据；环节动物实测 5.3–6.1），`--mem-gb` 之内调度，
大任务不会同时启动而互相挤爆内存；不足时自动重试（1.5 倍）并从检查点继续。

其他模式（宏基因组 meta、质粒 plasmid、单细胞 sc …）：`--spades-modes meta,plasmid`，每个模式各自组装并 BUSCO。

SPAdes 组装还用于：线粒体兜底与 [SPAdes-MitoFinder](Combo-mitogenome.md)、[核多拷贝基因](Combo-nuclear.md) 的拷贝数。

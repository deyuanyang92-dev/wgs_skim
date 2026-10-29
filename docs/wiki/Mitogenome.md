# 线粒体

## 组装策略（步骤 `mito`）

1. **MEANGS** 从 200 万条 reads 中找线粒体种子（`meangs_reads`）；找不到时 GetOrganelle 用自带的 animal_mt 种子。
2. **GetOrganelle 第 0 轮**：`-F animal_mt -R 10`，GetOrganelle 自己的读段子集。
3. **迭代**（`mt_iterate`）：第 0 轮没成环，就用上一轮结果作种子、用全部 clean reads（最多 `mt_iter_max_data_gb` = 20 Gb）、
   `-R 20` 再组，最多 `mt_iter_rounds` = 3 轮；某一轮成环、或结果与上一轮完全相同（再组也没用）时停止。
4. **SPAdes 兜底**（`mt_select`，需要 `spades`）：GetOrganelle 没成环、而 SPAdes 中对应的 scaffold 能闭合（图遍历 / 首尾重叠）时，
   用 SPAdes 的结果；否则保留 GetOrganelle 的。`.mt.json` 记录来源与原因。

```bash
wgs_skim run -r raw/ -o res/ --steps mito          # -> mito/<s>/<s>.mt.fasta, .mt.json
wgs_skim run -r raw/ -o res/ --steps getorganelle  # 不用 MEANGS
wgs_skim run -r raw/ -o res/ --steps mito --config mt_iter_max_data_gb=0   # 迭代用全部 reads
```

`mito/<s>/<s>.mt.json` 的 `stop` 说明为什么停止：`circular`、`no change after roundN`、`rounds exhausted`。

## 注释（步骤 `annotate`）

MitoM（<https://github.com/deyuanyang92-dev/Mitom>）：蛋白基因、rRNA、tRNA，以 cox1 为起点，带 NCBI table2asn 检查。
遗传密码、类群由 `--taxon` 预设（环节动物：密码表 5，Annelida-segmented-worms）。输出 `annotation/<s>/<s>.mt.gb`。

## 环还是线性？（步骤 `topology`）

只用 reads 才能判断环是真是假：在 nad5 中间切开，看有多少 reads 跨过原来的首尾连接点。

| 结论 | 含义 |
|---|---|
| `circular_validated` | 连接点被 reads 支持，环是真的 |
| `circular_weak` / `circular_repeat_junction` | 支持弱 / 连接点在重复区，长度可能不准 |
| `circular_rejected` | 首尾重叠但没有 reads 跨越 |
| `linear_gap_circular` | 首尾没接上，但有 read 对跨过缺口（实际是环，有缺口） |
| `linear_unjoined` | 末端没有接上；`mt_end_diagnosis` 说明原因（重复序列、同聚物、覆盖度下降） |
| `*_low_reads` | reads 太少无法判断 |

默认 `topology: auto` 只检查线性的线粒体和 SPAdes 闭合的环；`yes` 检查全部（包括 rDNA 单元）。

## 无 reads 环化：`wgs_skim circularize`

MitoFinder / SPAdes 的线粒体 contig 常常没闭合。不需要 reads，用 SPAdes 组装图闭合：

```bash
wgs_skim circularize -i <样品>_all_mt.fasta --spades spades_dir/ -o circ/ --prefix <样品>
```

四种方法（由强到弱）：组装图证明首尾相连 → 图中走一条线粒体深度的路径补上缺口（串联重复拷贝数按深度估计）→
首尾序列重叠 → 多条 contig 首尾相接成环。

## 两种组装是否一致：`wgs_skim compare`

GetOrganelle 结果作参考（它的组装图能标出重复区），把 SPAdes-MitoFinder（或任何其他组装）先环化、统一方向与起点、
全局比对，每个差异标注所在基因、是否在 GetOrganelle 的重复区或同聚物中；同时做 QUAST 式的结构检查（倒位、易位、缺失基因、移码）。

```bash
# 一个样品
wgs_skim compare -a getorganelle_dir/ -b mitofinder=<s>_all_mt.fasta --spades spades_dir/ --gb <s>.gb -o cmp/ --prefix <s>
# 一批（Snakemake 管理）：旧流程结果目录或样品表
wgs_skim compare --results results/ -o compare/
wgs_skim compare --sheet sheet.tsv -o compare/      # 列：sample ref other [spades gb gff gfa]
```

结论：`identical` 完全相同；`consistent` 差异都是二代数据无法确定的（重复区、线性末端、同聚物长度）；`incomplete`
没有错误但缺基因或缺片段；`differ` 有其他差异，需要人工检查。输出 `all.compare.tsv`、`all.diffs.tsv`、`all.genes.tsv`、
`all.structure.tsv`、`all.compare.html`（含每个样品的比较图 SVG）。

实测（MJ20260820041）：CHY006 identical、CHY007 consistent（MitoFinder 多出控制区 910 bp）、CHY008 consistent（两处 poly-T 长度）；
4 个样品 110 秒。

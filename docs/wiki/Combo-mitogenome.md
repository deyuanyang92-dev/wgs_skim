# 线粒体基因组

两种方法，推荐第一种；两种都做可以互相验证。

## 推荐：MEANGS-GetOrganelle（未成环自动迭代）

```bash
wgs_skim run -r raw/ -o res/ --steps meangs-getorganelle            # fastp → MEANGS → GetOrganelle → 迭代
wgs_skim run -r raw/ -o res/ --steps meangs-getorganelle,annotate   # + MitoM 注释
wgs_skim run -r raw/ -o res/ --steps meangs-getorganelle --config mt_iter_max_data_gb=0   # 迭代用全部 reads
```

1. **MEANGS** 从 200 万条 reads 找到样品自己的线粒体种子（找不到时用 GetOrganelle 自带的 animal_mt 种子）。
2. **GetOrganelle 第 0 轮**：`-F animal_mt -R 10`。
3. **迭代**：没有成环，以上一轮结果作种子、用全部 clean reads（≤ `mt_iter_max_data_gb` = 20 Gb）、`-R 20` 再组装，最多
   `mt_iter_rounds` = 3 轮；成环或结果不再变化时停止。
4. **SPAdes 兜底**（同时运行 SPAdes 时）：GetOrganelle 仍未成环、而 SPAdes 中对应的 scaffold 能闭合时，用 SPAdes 的结果。

输出 `mitogenome/<s>/<s>.mt.fasta`；`<s>.mt.json` 记录来源、轮次与停止原因（`circular`、`no change after roundN`、`rounds exhausted`）。

## SPAdes-MitoFinder

全基因组组装中通常也有完整的线粒体 contig。MitoFinder 从 SPAdes 组装中找出它并注释。

```bash
wgs_skim run -r raw/ -o res/ --steps spades-mitofinder --mitofinder-ref related.gb
# 与 MEANGS-GetOrganelle 一起：参考默认用该样品自己的 MitoM 注释（不必给 --mitofinder-ref）
wgs_skim run -r raw/ -o res/ --steps meangs-getorganelle,spades-mitofinder
# 已有 SPAdes 组装：直接用 MitoFinder 工具
wgs_skim mitofinder -i spades_assemblies/ -o mf/ --ref related.gb
```

输出 `mitofinder_spades/<s>/<s>.mitofinder.fasta`、`.gb`，`mitofinder_spades/mitofinder_summary.tsv`（每个样品每个基因的长度），
`mitofinder_spades/genes/<基因>.fasta`。详见 [MitoFinder](Tool-MitoFinder.md)。

## 两种结果是否一致：`wgs_skim compare`

以 GetOrganelle 结果为参考，把 MitoFinder（或任何其他组装）先环化、统一方向与起点、全局比对；每个差异标注所在基因、
是否位于 GetOrganelle 的重复区或同聚物中，并做结构检查（倒位、易位、缺失基因、移码）。

```bash
wgs_skim compare -a res/mito/S/round1 -b mitofinder=res/mitofinder_spades/S/S.mitofinder.fasta \
    --gb res/annotation/S/S.mt.gb -o cmp/ --prefix S                   # 一个样品
wgs_skim compare --results results/ -o compare/                       # 一批（Snakemake 管理）
wgs_skim compare --sheet sheet.tsv -o compare/                        # 列：sample ref other [spades gb gff gfa]
```

结论：`identical` 完全相同；`consistent` 差异都在二代数据无法确定的地方（重复区、线性末端、同聚物长度）；
`incomplete` 无错误但缺基因 / 片段；`differ` 需要人工检查。实测（MJ20260820041）：CHY006 identical、CHY007 consistent、
CHY008 consistent，4 个样品 110 s。

## 注释：MitoM

`--steps annotate`（或单独 `wgs_skim mitom -i`）：蛋白基因、rRNA、tRNA，cox1 起点，NCBI 检查，输出 `mitom/<s>/<s>.mt.gb`。见 [MitoM](Tool-MitoM.md)。

## 环还是线性？reads 验证

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
wgs_skim circularize -i res/mitofinder_spades/S/S.mitofinder.fasta --spades res/assembly/S/spades -o circ/ --prefix S
```

四种方法（由强到弱）：组装图证明首尾相连 → 图中走一条线粒体深度的路径补上缺口（串联重复拷贝数按深度估计）→
首尾序列重叠 → 多条 contig 首尾相接成环。

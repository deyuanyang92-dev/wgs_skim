# 核多拷贝基因：rDNA 与组蛋白

## 两种方法，互相验证

| 方法 | 怎么做 | 优点 |
|---|---|---|
| GetOrganelle（步骤 `nr`） | 从 reads 直接组装串联重复单元（`-F anonym`），rDNA 用 `resources/rdna` 种子，组蛋白用 `resources/histone`（环节动物 H1–H4 基因簇拷贝 + Hesionidae H3） | 得到完整单元：基因顺序、间隔区 |
| SPAdes-MitoM（步骤 `nuclear`） | MitoM 在 SPAdes 全基因组组装中找出全部拷贝，每个基因合并成一条一致序列 | 拷贝数、与基因组整体比较 |

同一个基因组两种方法应该给出相同的基因；`nuc_compare.tsv` 逐基因比较。

## 数据量：多拷贝基因只需少量数据

- 第 1 轮最多用 **5 Gb** reads。
- 单元**不完整**时，第 2 轮用 **10 Gb**，并以第 1 轮的结果作种子（迭代）。
- "完整" = 闭合成环，或核心基因（rDNA：18S、5.8S、28S；组蛋白：H2A、H2B、H3、H4）都被标签序列完整覆盖。
  线性但基因齐全的 rDNA 单元算完整，不再跑第 2 轮（实测：第 2 轮对这类样品没有任何改进，每个白费 15–40 分钟）。
- 每轮的结果按目标基因覆盖的碱基数打分；第 2 轮如果"闭合"成一个不含目标基因的短重复环（测试中出现过 155 bp 的简单重复），不会被采用。
- `<marker>/getorganelle/<s>/rounds.tsv` 记录每一轮：数据量、种子、结果、每个基因的覆盖度、是否闭合、是否被采用。
- 修改：`--max-gb 5,10`（`wgs_skim nuclear`）或配置 `nuclear_max_gb: "5,10"`。

组蛋白的种子 reads 用 bowtie2 `--very-sensitive-local` 招募：其他科（如 *Pisione*）的组蛋白 DNA 与自带种子差异较大，
默认的端到端比对几乎找不到种子 reads。

## 组蛋白的类别

MitoM 把每个组蛋白拷贝分类，每个基因的代表序列优先取基因簇拷贝：

| 类别 | 判断 | 说明 |
|---|---|---|
| `replication-dependent` | H3 第 31/87/89/90 位 = A S V M，或覆盖度 ≥ 5 倍单拷贝 | 串联基因簇拷贝，系统发育常用的 H3 / H4 |
| `replacement` | H3.3（S A I G），或含内含子 | 替换变体，另一个旁系同源基因 |
| `single-copy` | 覆盖度 < 2 倍单拷贝（如 H2A.X） | 不在基因簇里 |

## 在流程中运行

```bash
wgs_skim run -r raw/ -o res/ --steps nr        # GetOrganelle + MitoM -> rdna/<s>/<s>.rdna.gbf、histone/<s>/<s>.histone.gbf
wgs_skim run -r raw/ -o res/ --steps nuclear   # + SPAdes-MitoM、拷贝数、比较 -> nuclear_markers/<s>/
```

## 对已有结果运行：`wgs_skim nuclear`

旧流程或本流程的结果目录都可以（自动识别结构），已有的 rDNA 组装算作第 1 轮，SPAdes 组装直接使用。由 Snakemake 管理：

```bash
wgs_skim nuclear --results results/ -o nuclear/                  # 全部样品
wgs_skim nuclear --results results/ -o nuclear/ --samples A,B    # 部分样品
wgs_skim nuclear --results results/ -o nuclear/ -j 64 --mem-gb 200 -n   # 干运行
wgs_skim nuclear --status -o nuclear/                            # 进度（运行中也可以）
```

输出：

```
nuclear/
├── summary.tsv                    每个样品 x 每个基因：两种方法的长度、是否不完整、组蛋白类别、错配、插入缺失、一致度、结论
├── rdna/ histone/
│   ├── getorganelle/<s>/          组装结果、组装图、rounds.tsv
│   └── mitom/<s>/            <s>.rdna.gbf / <s>.histone.gbf、nuc_features.tsv、consensus/、genes/、submission/（NCBI）
├── spades_mitom/<s>/              MitoM 在 SPAdes 组装上的结果（两种标记一次完成）
├── compare/<s>.tsv                每个样品的比较
└── logs/ benchmarks/              每个任务的日志、时间、内存
```

`summary.tsv` 的结论：`identical` 相同；`contained` 较短（不完整）的一条完全包含在另一条中；`differ` 有差异（组蛋白不同拷贝之间可能有少量差异）；
`no overlap` 几乎没有可比部分；`only GetOrganelle` / `only SPAdes-MitoM` 只有一方找到。

实测：小样品 rDNA 5 个基因两种方法完全相同；组蛋白 H2B、H3、H4 相同、H2A 1 个碱基差异、H1 只在 SPAdes 中（单拷贝）。
MJ20260820041 22 个样品 rDNA 注释全部通过 NCBI table2asn 校验（0 错误）。

# MitoM 排错参考

## 1. 简介与用途

MitoM（命令 `mitom`，主程序 `batch_mitoz.py`，文件名沿用早期版本、已不调用 MitoZ）是用户自己的动物线粒体注释流程：
原生注释核心 `mito_core.py`（blastx + exonerate 蛋白基因、Infernal/HMMER rRNA、MiTFi 等 tRNA）→ cox1 定向 → 定向后有问题
才第二轮注释 → 补搜 / 内含子 → 端点校正与 NCBI 命名 → table2asn 校验。`mitom nuc`（`mito_nuc.py`）注释 contig 中的核 rDNA
（18S、ITS1、5.8S、ITS2、28S）和组蛋白（H1、H2A、H2B、H3、H4）：BLAST 投影 NCBI 参考的边界，再用 Rfam 模型精修 rRNA 端点。

**在 wgs_skim 中的作用**（取代 MitoZ，ISSUES #46）：

| 规则 | 命令 | 输入 → 输出 |
|---|---|---|
| `mitom_annotate`（annotate.smk） | `mitom run` | 最终线粒体 `mito/<s>/<s>.mt.fasta` → `annotation/<s>/<s>.mt.gb`、`<s>.mitom.tsv` |
| `nr_annotate`（nuclear.smk） | `mitom nuc annotate --genes rdna` | GetOrganelle rDNA 单元 → `rdna/<s>/<s>.rdna.gbf`、`nr_features.tsv` |
| `histone_annotate` | `nuclear_markers.py annotate`（内部调 `mitom nuc annotate --genes histone`） | 组蛋白簇 → `histone/<s>/<s>.histone.gbf` |
| `nuc_markers` | `mitom nuc annotate --genes rdna,histone` | SPAdes `scaffolds.fasta` → 每个拷贝 + 每基因共识（`nuc_features.tsv`、`consensus.tsv`） |

MitoM 的 cox1 注释也是 BLAST 物种核对的默认查询（`blast_mt: cox1`），`mitom nuc` 的 18S/ITS/28S 是 rDNA 查询。
MitoFinder 规则在没有 `mitofinder_ref` 时用 MitoM 的 `.mt.gb` 作参考。

- GitHub：https://github.com/deyuanyang92-dev/Mitom（本地 `~/mybiosofts/Mitom`，main = origin/main = 6f75f73）
- 文档：仓库 `README.md`（中英，1041 行）、`docs/benchmark_*.md`、`knowledge/topics/nuclear_markers.md`
- 方法依据：Yang et al. 2024, *Zoosystematics and Evolution* 100(4): 1297–1314, doi:10.3897/zse.100.127201（COX1 起点；
  tRNA 优先级 tRNAscan-SE > MiTFi > ARWEN 是论文时的做法，MitoM 现在用 MiTFi / MitoM-ML，见下）；MITOS2 参考数据
  refseq89m（Zenodo 4284483，md5 bb6325e2…）。无独立论文。

## 2. 版本（git commit）

| 用途 | 配置键 | 固定 commit | 日期 |
|---|---|---|---|
| 线粒体注释 `mitom run` | `mitom_commit` | **27264d2**（Duplicated gene blocks, gene-order ORF fallback） | 2026-10-23 |
| 核基因 `mitom nuc` | `mitom_nuc_commit` | **6f75f73**（PR #2 合并，= 当前 main） | 2026-10-28 |

两次克隆同一仓库（rule `mitom_db` → `MITOM_DB/Mitom` + refseq89m；rule `mitom_nuc_db` → `MITOM_NUC_DB/Mitom`）；
`mitom_dir` 非空时两者都用该本地目录。

27264d2 之后、6f75f73 之前的 15 个 commit 与线粒体注释有关的有：`mito_nuc` 首版（8a0debe）；rRNA 链检查的 HMM 库、第二拷贝、
被序列末端切断的基因（5368d07）；MitoM-ML tRNA（4e9b99e、f0e815d）；**NCBI 可提交输出、MiTFi 按 MITOS2 方式运行、tRNA 第二意见
（6575782）**；`--rrna_method hmm` 默认（f1bc9b0，450 条留一法：找回率与 cm 相同 98.5 %，更快更省内存）；去掉 MitoZ、原生圈图
（95a2579）；PR #1（组蛋白拷贝分类）、PR #2（contig 末端无比对 rRNA 片段的保留条件）。

**是否升级 `mitom_commit`：值得，但要先补环境再测。** 6f75f73 的 `env/mitotools.yml` 比 wgs_skim `envs/mitom.yaml` 多了
`matplotlib`（圈图）、`scikit-learn=1.6.1`（MitoM-ML tRNA）、`bwa`、`defusedxml`、`pytest`。缺这些时 MitoM **不会报错**：
圈图记 `[PLOT] mito_plot.py not usable` 并跳过，tRNA 记 `MitoM-ML not installed ...: tRNAs by MiTFi` 退回 MiTFi
（源码 `batch_mitoz.py` run_plots、`mito_core.py` 第 928 行）。所以现在设 `mitom_dir: ~/mybiosofts/Mitom` 时，线粒体注释是在
“缺依赖的 HEAD”上跑的，与 MitoM 自己的基准不是同一配置。建议：`envs/mitom.yaml` 加 matplotlib、scikit-learn=1.6.1（bwa 可选），
在 HNP1sub + CHY006/007/008 上对比 27264d2 与 6f75f73 的 `summary.tsv` / 基因边界，通过后把两个 commit 合成一个（只克隆一次）。

## 3. wgs_skim 中的调用方式

```bash
# 线粒体（8 线程，4 GB）；ID 先改成 [A-Za-z0-9_.-]："S|mt|1" -> "S_mt_1"
python <Mitom>/mitom run -i input.fasta -o out --clade Annelida-segmented-worms --genetic_code 5 \
    --threads 8 --max_tasks 1 --id_maxlen 100 --mitos_refdir <db>/mitos_refs
cp out/final_annotation_all.gb annotation/<s>/<s>.mt.gb; cp out/summary.tsv annotation/<s>/<s>.mitom.tsv

# rDNA 单元（4 线程）/ SPAdes 组装（16 线程，8 GB）
python <Mitom>/mitom nuc annotate -i input.fasta -o out --genes rdna --threads 4
python <Mitom>/mitom nuc annotate -i scaffolds.fasta -o mitom_nuc --genes rdna,histone --threads 16
```

| 参数 | 含义 | 配置键 |
|---|---|---|
| `--clade` | 蛋白库（Annelida-segmented-worms、Mollusca、Chordata、Arthropoda…、Animal） | `clade`（空 = `taxon` 预设） |
| `--genetic_code` | 线粒体遗传密码 | `genetic_code` |
| `--mitos_refdir` | MITOS2 refseq89m 目录 | `mitos_refdir`（空 = 下载一次） |
| `--max_tasks 1` | 一次一条记录（每个样品一个 Snakemake 任务） | 固定 |
| `--id_maxlen 100` | 内部 ID 最大长度（≥ 6） | 固定 |
| `--genes` | rdna / histone / rdna,histone | 按规则 |

## 4. 重要参数说明

| 参数 | 默认 | 说明 |
|---|---|---|
| `--clade` | auto | FASTA 输入没有分类信息时**沿用 Annelida 并警告**——非环节动物必须显式给（wgs_skim 总是给） |
| `--genetic_code` | 随 clade（Chordata 2、Echinodermata 9、Porifera 4，其余 5；Cnidaria / Platyhelminthes / Chaetognatha / Animal 必须给） | 显式值优先 |
| `--which_gene_first` / `--fallback_genes` | cox1 / cox2,cox3,cob,nad1 | 定向基因；线性序列从不旋转 |
| `--if_reannotation_2` | auto | 定向后只重注检查未通过的记录（`reports/reannotation_decision.tsv`） |
| `--rrna_method`（6f75f73） | hmm | `cm` = 只用 MITOS2 协方差模型 |
| `--rescue` | auto（= 不运行） | `yes` 需另装 MITOS2 |
| `--plot`（6f75f73） | yes | 需要 matplotlib；`--infq` 加深度环 |
| `--keep_intermediate` | no | 排错时 `yes`（中间目录移到 `work/`）；问题样本核心日志总在 `reports/failed_logs/` |
| `nuc --consensus_identity` | 0.99 | 与共识一致度低于此的片段（他种 / 污染 / 另一类 rDNA）不进共识 |
| `nuc --no_cm` | 关 | 只用参考投影，不用 Rfam 模型精修端点 |
| `nuc --lib` | 随附 `models/nuc/`（NCBI 2026-10-23） | `mitom nuc fetch` / `build` 可加入自己的注释（如 Geneious GenBank） |
| `nuc --organism` / `--src` | — | 写进 NCBI 提交文件（`submission/nuc.sqn`） |

## 5. 已知问题（仓库 issue 区目前只有两个 PR，#1 #2，均已合并）

| 现象 | 原因 | 解决 | 链接 |
|---|---|---|---|
| ITS1 被截短：contig 末端 20 nt 被注成“5′ 不完整 18S” | 无比对支持的末端片段被保留 | 只有与参考基因末端相似 ≥ 70 % 才保留 | [PR #2](https://github.com/deyuanyang92-dev/Mitom/pull/2)（TESTING #4，CHY007） |
| H3 共识实为单拷贝 H3.3，H2A 实为 H2A.X | 组蛋白拷贝未分类，共识取了完整的单拷贝变体 | 拷贝分为 replication-dependent / replacement / single-copy，共识优先取基因簇拷贝 | [PR #1](https://github.com/deyuanyang92-dev/Mitom/pull/1)（TESTING #5） |
| 非环节动物结果差 | FASTA 无分类信息时默认 Annelida 蛋白库 | 显式 `--clade` `--genetic_code` | README 注意事项 1–2 |
| 圈图缺失 / tRNA 只用 MiTFi | 环境缺 matplotlib / scikit-learn（静默降级） | 补依赖（第 2 节） | README 依赖；源码 |
| SPAdes 组蛋白簇只剩 < 150 nt 的碎片 | SPAdes 把高拷贝簇拆碎 | 报告簇碎片（标不完整），不拿单拷贝变体顶替；完整单元看 GetOrganelle 组蛋白 | README nuc 节 |
| 参考稀少类群（海绵、刺胞、扁形）rDNA 边界多靠外推 | 库中参考少 | `mitom nuc build` 加入本类群注释 | README nuc 节 |
| 线性序列 | 定向后仍保持 linear，不强制 circular | — | README 注意事项 4 |

## 6. 本项目中遇到的问题

1. **MitoZ 3.6 → MitoM（ISSUES #46）**：MitoZ 的 tbl2asn 不接受 `HNP1sub|mt|1` 这样的 ID（"Unable to make sequence identifier"），
   不生成 .gbf；MitoM 在同一条 15 662 bp 线粒体上 15 s，13 PCG / 22 tRNA / 2 rRNA，NCBI 0 错误。规则仍先把 ID 中 `|` 等字符换成 `_`
   （topology 按序列内容对应注释）。
2. **TESTING #4、#5** 见上表（MitoM PR #2、#1）；#12：组装更新后旧注释被复用 → 注释标记中加入输入文件大小与时间。
3. **rDNA 批量（TESTING 3.5）**：MJ20260820041 22 个 `.rdna.gbf` 全部 table2asn 0 错误 0 警告；HNP1sub rDNA 5 个基因 GetOrganelle 与
   SPAdes-MitoM 完全相同（18S 1882、ITS1 377、5.8S 158、ITS2 226、28S 3773 nt）。
4. **`mitom_db` 需要网络**：`git clone` 失败时用 `mitom_dir: ~/mybiosofts/Mitom` + `mitos_refdir:`（README Troubleshooting）；
   但见第 2 节：本地 main 与 `mitom_commit` 不是同一版本。
5. 资源：`mitom run` 一条线粒体 14–15 s / 0.2 GB（CWSC010 140 s，多 contig 输入）；`nuc annotate` rDNA 单元 2 s / 0.13 GB；
   SPAdes 组装 15–118 s / 0.9–2.2 GB；MitoM-ML 首次运行时在 `models/trna_ml/cache/` 训练并缓存模型（仓库目录需可写）。

## 7. 排错清单

- 日志：`logs/mitom/<s>.log`、`logs/nr_annotate/<s>.log`、`logs/histone_annotate/<s>.log`、`logs/nuc_markers/<s>.log`；
  完整过程在工作目录 `annotation/<s>/mitom/out/run.log`，问题记录的核心日志在 `out/reports/failed_logs/`。
- 空输出：日志只有 `no mitogenome` / `no rDNA` / `no histone cluster` / `no assembly` = 上游没有序列（先查 GetOrganelle / SPAdes）。
- 关键字：`[Step k/11]`（卡在哪一步）、`WARNING ... Annelida`（clade 回退）、`[PLOT] mito_plot.py not usable`、
  `MitoM-ML not installed`、`round 2`（`summary.tsv` annotation_used）、`oriented_by_fallback.txt`（没有 cox1）、
  `submission_summary.tsv` ERROR 数；nuc：`consensus.tsv` 的 dropped 原因、`nuc_features.tsv` 的 `evidence`（hmm/cm/aligned/projected）
  与 `aligned_nt`（越少越要人工核对）、`class`（组蛋白）。
- 环境自检：`python <Mitom>/mitom doctor`（所有必需程序应为 OK）；版本：`git -C <Mitom> rev-parse --short HEAD`。
- 复现：在 mitom 环境里对 `annotation/<s>/mitom/input.fasta` 重跑规则命令，加 `--keep_intermediate yes`。

## 8. 资料来源与查阅日期（2026-10-29）

- 本地仓库 `~/mybiosofts/Mitom`：`README.md`（快速上手、关键参数、输出目录、核基因、注意事项）、`git log`（27264d2..6f75f73）、
  `env/mitotools.yml` 差异、`batch_mitoz.py`（run_plots）、`mito_core.py`、`trna_ml.py`、`mitom -h`、`mitom nuc annotate -h`、
  `docs/benchmark_six_tools_2026-10-22.md`
- GitHub PR #1、#2（gh api；issue 区无 issue）
- 本项目：`workflow/rules/annotate.smk`、`nuclear.smk`、`scripts/nuclear_markers.py`、`envs/mitom.yaml`、`config/config.yaml`、
  `docs/TESTING.md` §3.2、3.5、4、README Troubleshooting；MJ20260820041 `ISSUES.md` #46；`wgs_skim_test/*/benchmarks/`

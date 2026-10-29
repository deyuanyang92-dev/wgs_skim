# MitoFinder 排错参考

## 1. 简介与用途

MitoFinder 用参考线粒体（GenBank，可多条）在组装中找线粒体 contig（BLAST），再按参考注释蛋白基因 / rRNA（相似性 +
找起止密码子）、用 MiTFi / ARWEN / tRNAscan-SE 注释 tRNA；也可先用 MEGAHIT / metaSPAdes / IDBA-UD 组装 reads。

**在 wgs_skim 中的作用**（rule `mitofinder`，`workflow/rules/mitofinder.smk`，正在加入）：从一个组装中找出并注释线粒体，
作为 GetOrganelle 之外的第二条独立证据（用户 2024 年 *Perinereis* 论文的做法：多组装器 + MitoFinder 跑在 de novo 组装上）。
组装来源 `mitofinder_assembly`：样品表给的组装 / SPAdes 全基因组组装（`-a`），或让 MitoFinder 自己组装 reads
（`--megahit` / `--metaspades`）。结果再交给 `wgs_skim circularize`（无 reads 环化）和 `wgs_skim compare`（与 GetOrganelle 比较）。

- GitHub：https://github.com/RemiAllio/MitoFinder；容器：https://github.com/RemiAllio/MitoFinder_container
- 论文：Allio R, Schomaker-Bastos A, Romiguier J, Prosdocimi F, Nabholz B, Delsuc F (2020) MitoFinder: Efficient
  automated large-scale extraction of mitogenomic data in target enrichment phylogenomics. *Mol Ecol Resour* 20: 892–905.
  doi:10.1111/1755-0998.13160
- 用户旧脚本：`Temp_scripts/mitofinder-batch/mitofinder_pipeline_v2_7.py`（annotation2summary / assembly2summary / summary）

## 2. 版本

| | 版本 | 说明 |
|---|---|---|
| bioconda | **1.4.1**（2021-10；py27；依赖 blast、megahit、spades、idba、trnascan-se、openjdk、biopython、pandas，**均未固定版本**） | `envs/mitofinder.yaml` 默认 |
| 最新 / Singularity | **1.4.2**（2024-04-11），`singularity pull --arch amd64 library://remiallio/default/mitofinder:v1.4.2` | 本机 `~/mybiosofts/mitofinder_v1.4.2.sif`（986 MB）；配置 `mitofinder_sif` |

1.4.1 → 1.4.2（GitHub compare，源码逐行确认）：

1. **环化检查修复**：1.4.1 里 `Popen(args1, ..., shell=True)` 把参数列表交给 shell，只执行了 `circularizationCheck.py`
   本身、没有参数 → 环化检查形同虚设（作者原话 "deprecated (involuntarily)"）。1.4.2 去掉 `shell=True`。
2. 每条 contig 生成注释示意图 `*_mtDNA_contig.png`（需 pillow 6.2.2；bioconda 1.4.1 的依赖里没有 pillow）。
3. 所有脚本 shebang 改为 `#!/usr/bin/python2.7`；COX1/2/3 名称识别多容忍 "-"。

**是否升级：应当默认用 1.4.2 镜像。** bioconda 1.4.1 的 “Circularization: Yes/Not found” 不可信，环化后的首尾修剪也
不会发生；作者在 README 和多个 issue 中只推荐 Singularity（#47 #48 #63）。conda 1.4.1 只适合没有 Singularity 的机器，
且环化结论一律以 `wgs_skim circularize` / topology 为准。2026-10-29 对 `mitofinder=1.4.1` 做 mamba dry-run 可解：
python 2.7.15、blast 2.16.0、megahit 1.2.9、spades 3.15.5、trnascan-se 2.0.13、openjdk 11（本机旧环境 `envs/mitofinder`
只有 python 2.7，没装上 mitofinder）。

## 3. wgs_skim 中的调用方式

```bash
# 组装已有（given / spades）：只做查找 + 注释
mitofinder -j <s> -a <scaffolds.fasta> -r <ref.gb> -o <genetic_code> -p 8 -m <GB> -t mitfi --adjust-direction
# 让 MitoFinder 组装 reads
mitofinder -j <s> -1 R1.clean.fq.gz -2 R2.clean.fq.gz --megahit -r <ref.gb> -o <code> -p 8 -m 32 -t mitfi --adjust-direction
# Singularity（mitofinder_sif 非空时）：
singularity exec -B <输入目录...> -B <工作目录> mitofinder_v1.4.2.sif mitofinder ...
```

规则在 `mitofinder/<s>/run/` 里运行，完成后把 `<s>/*_Final_Results` 移到 `mitofinder/<s>/Final_Results`，合并
`*_mtDNA_contig*.fasta`（不含 `_genes_`）为 `<s>.mitofinder.fasta`、`.gb` 为 `<s>.mitofinder.gb`；自组装时保留
`<s>_link_<assembler>.scafSeq` → `<s>.<assembler>.contigs.fasta`，然后删 `run/`。

| 参数 | 含义 | 配置键 |
|---|---|---|
| `-j` | 样品名（输出目录 / 文件前缀） | 样品名 |
| `-a` / `-1 -2 --megahit|--metaspades` | 已有组装 / 自己组装 | `mitofinder_assembly`（given / spades / megahit / metaspades；空 = 自动） |
| `-r` | 参考线粒体 GenBank | `mitofinder_ref`；空 = 本样品 GetOrganelle 线粒体的 MitoM 注释 |
| `-o` | NCBI 遗传密码（5 无脊椎线粒体、2 脊椎、9 棘皮/扁形、13 海鞘…） | `genetic_code`（空 = 按 `taxon` 预设） |
| `-t` | tRNA 工具 mitfi / arwen / trnascan | `mitofinder_trna`（默认 mitfi） |
| `-p` / `-m` | 线程 / 组装内存（GB，只对 MEGAHIT/metaSPAdes 有意义） | rule threads / `mitofinder_mem_gb`（-a 模式规则只申请 8 GB） |
| `--adjust-direction` | 按参考调整 contig 方向 | 固定 |

## 4. 重要参数说明

| 参数 | 默认 | 说明 |
|---|---|---|
| `-r` | — | 必须带注释（CDS/rRNA）；只有序列没有特征时报 "didn't found any nucleotide sequence in the reference"（#48）。可给多条，MitoFinder 每个基因选最佳参考（#10） |
| `--min-contig-size` | 1000 | 小于此长度的 contig 不考虑（捕获 / 低深度数据要调低） |
| `--max-contig-size` | 25000 | 大于此长度不考虑（README：长线粒体或嵌合 scaffold 会被丢掉） |
| `--max-contig` | 0（不限） | 只保留前 N 条匹配 contig；**只有找到 1 条（或 `--max-contig 1`）时才做环化检查** |
| `--circular-size` / `--circular-offset` | 45 / 200 | 首尾 200 bp 内找 ≥ 45 bp 的自身 BLAST 重叠才判环（#58 #62） |
| `--blast-identity-nucl` / `-prot` / `--blast-size` / `-e` | 50 / 40 / 30 % / 1e-5 | 远缘参考找不到 contig 时调低（#23） |
| `--new-genes` / `--ignore` | 关 | 非标准基因 |
| `--allow-intron` / `--numt` | 关 | 需近缘参考，建议先找 contig 再第二步 `-a` 注释 |
| `--override` | 关 | 重跑自组装时必须删旧目录或加此项，否则跳过组装 |
| `-l` | 100 | MEGAHIT 最短 contig |

组装器：README 写 MEGAHIT 默认更快、metaSPAdes “a bit slower but more efficient”（不支持单端）。

## 5. 已知问题与 issues

| 现象 / 报错 | 原因 | 解决 | 链接 |
|---|---|---|---|
| `Evidences of circularization could not be found, but everyother step was successful` | 1.4.1 检查失效；或末端停在控制区重复，首尾重叠不在 200 bp 窗口 / 外侧有不一致 | 用 1.4.2；换 metaSPAdes；用 wgs_skim circularize / topology 判断 | [#36](https://github.com/RemiAllio/MitoFinder/issues/36)、[#62](https://github.com/RemiAllio/MitoFinder/issues/62)、[#58](https://github.com/RemiAllio/MitoFinder/issues/58) |
| 环化检查时报 Traceback（subprocess） | 容器内 SPAdes 3.13 / 环境问题 | 附 log + geneChecker.log | [#26](https://github.com/RemiAllio/MitoFinder/issues/26) |
| `invalid syntax`（print 语句） | 用 python3 运行了 python2.7 代码 | `python2.7 mitofinder` 或用容器 | [#54](https://github.com/RemiAllio/MitoFinder/issues/54)、[#19](https://github.com/RemiAllio/MitoFinder/issues/19)、[#55](https://github.com/RemiAllio/MitoFinder/issues/55) |
| 容器里 python 库报错 | 宿主 conda 的 PYTHONPATH 等泄漏进容器 | `singularity exec --cleanenv`，不在激活的 conda 环境里调用 | [#67](https://github.com/RemiAllio/MitoFinder/issues/67) |
| MEGAHIT `No such file or directory: /tmp/k29/...` | conda 版 MEGAHIT 临时目录问题 | 换 Singularity 后成功 | [#47](https://github.com/RemiAllio/MitoFinder/issues/47) |
| `ERROR: MetaSPAdes didn't run well` / metaSPAdes 卡住 | 原因未统一（内存 / reads / 容器内 SPAdes 3.13，#17 作者也无解） | 改用自己的 SPAdes 组装再 `-a`（作者在 #17 的建议） | [#42](https://github.com/RemiAllio/MitoFinder/issues/42)、[#17](https://github.com/RemiAllio/MitoFinder/issues/17) |
| 组装 > 18 h | 全部 WGS reads 交给 MEGAHIT | 作者建议抽 5–10 M reads；已有组装就用 `-a` | [#13](https://github.com/RemiAllio/MitoFinder/issues/13) |
| BLAST 有命中但 Final_Results 为空 | 阈值（identity / size / e 值） | 调低 `--blast-*` | [#23](https://github.com/RemiAllio/MitoFinder/issues/23) |
| 缺某基因（如 COX1），blastn e 值为 0.0 | e 值取整 / 多 contig 选择 | 只给候选 contig 再跑 | [#60](https://github.com/RemiAllio/MitoFinder/issues/60) |
| 长度不一、缺基因、基因找到多次 | 片段化、NUMT、污染（日志 WARNING "genes were found more than once"） | 逐 contig 检查 | [#51](https://github.com/RemiAllio/MitoFinder/issues/51) |
| 没有 `_final_genes_NT.fasta` / png | geneChecker 出错；1.4.1 无 png | 看 geneChecker(_error).log；用容器 | [#43](https://github.com/RemiAllio/MitoFinder/issues/43)、[#63](https://github.com/RemiAllio/MitoFinder/issues/63) |
| `makeblastdb is not executable` / Ubuntu 22.04 blast 问题 | 自带 BLAST 二进制 | 容器或系统 BLAST | [#61](https://github.com/RemiAllio/MitoFinder/issues/61)、[#44](https://github.com/RemiAllio/MitoFinder/issues/44) |
| `Feature references another sequence` / `NameError: featureName` | 多半是参考 GenBank 的特征写法（跨记录 join、非标准名）——只看了标题，未核实 | 换参考或清理特征 | [#66](https://github.com/RemiAllio/MitoFinder/issues/66)、[#50](https://github.com/RemiAllio/MitoFinder/issues/50) |
| tmp.gb 里只有最后一条参考 | 正常（逐条处理） | — | [#10](https://github.com/RemiAllio/MitoFinder/issues/10) |

## 6. 本项目中遇到的问题

1. **MitoFinder 判环不可靠（ISSUES #36）**：只认首尾 200 bp 内 ≥ 45 bp 的重叠。CHY006 SPAdes contig NODE_1248（16 828 bp）
   两端 58 bp 重叠前 24 bp 不一致 → `Circularization: Not found`；而那 58 bp 本身落在 94 bp 串联重复里，按重叠闭环会丢 2 个
   重复单元（GetOrganelle 17 052 bp，3 拷贝）。CHY007 报 `Yes` 并把 scaffold 从 17 305 截到 17 040 bp（ISSUES #50）。
   → 环 / 线性一律看 `wgs_skim circularize`（图游走 + 深度补拷贝）和 topology（reads 跨接）。结果：CHY006 17 052 bp（与
   GetOrganelle 100 %）、CHY008 17 048 bp、CHY007 17 228 bp（重复拷贝数不确定）。
2. **整个 SPAdes 组装太大**：旧流程（`spades_mitofinder.py`）先预筛 scaffold：长度 ≥ 1000 且 k-mer 覆盖度 ≥ 5 × 长度加权
   中位数，CHY006 剩 3871 条（10 MB），MitoFinder 5 min 08 s、只找到 1 条 contig（于是才做了环化检查）。注释说明
   “对整个 2 GB 组装做 BLAST 很慢”。当前 rule 直接给完整 `scaffolds.fasta`，未预筛（见第 7 节建议）。
3. **参考用本样品的 MitoM 注释**：旧运行 `-r CHY006_0729_locus_mt_1.gb`（MitoM `final_annotation`），15 个基因（13 PCG +
   rrnL/rrnS）全部识别，tRNA 由 MiTFi 完成 → MitoM 的基因命名与 MitoFinder 兼容。但这是“自参考”：适合找同一线粒体在另一
   组装中的版本，不能当独立的物种参考；没有 GetOrganelle 结果时需给近缘种 GenBank（MitoM 基准：缺近缘参考的 2/24 条
   MitoFinder 失败；环节动物基准中只用带全部基因的记录作参考）。
4. 旧脚本 v2.7 已处理：多样品合并 FASTA 按 header 第一列拆分（重复 header 改名 `<s>_contigN`）；运行前检查 `mitofinder`
   可执行（否则几十个样品同时 `No such file or directory`）；缺 Biopython 时报 `Biopython 未安装`（TROUBLESHOOTING）。
5. 准确度参考（MitoM 基准，17 条共同记录）：MitoFinder 1.4.2 蛋白基因两端完全一致 59.2 %、缺基因 8、table2asn ERROR 17；
   中位 204 s、238 MB。wgs_skim 用 MitoM 做最终注释，MitoFinder 只用于“找 contig”。

## 7. 排错清单

- 日志：`logs/mitofinder/<s>.log`（MitoFinder 标准输出）；运行目录中途失败时保留在 `mitofinder/<s>/run/<s>/`：
  `<s>_MitoFinder.log`、`geneChecker.log`、`geneChecker_error.log`、`<s>_blast_out.txt`。
- 关键字：
  - `MitoFinder found a single mitochondrial contig` / `found N contigs matching` + `Did not check for circularization`（多条 = 不判环，
    多为 NUMT 或片段）。
  - `Evidences of circularization were found!`（会按重叠修剪）/ `could not be found` / `.infos` 里 `Circularization: Yes|Not found`。
  - `dit not found any mitochondrial contig with a size lower than 25000 bp` —— 被 `--max-contig-size` 挡掉。
  - `didn't found any nucleotide sequence in the reference` —— 参考没有特征。
  - `genes were found more than once` —— 片段化 / NUMT / 污染。
  - `tRNA annotation with MitFi run well` —— tRNA 步正常。
  - `[wgs_skim] MitoFinder found no mitochondrial contig` / `empty assembly` / `MitoFinder failed`。
- 版本：`mitofinder -v`（1.4.1 = 环化检查失效）；容器内路径 `/opt/MitoFinder/mitofinder`。
- 绑定目录：容器只看得到 `-B` 的目录；输入是软链接时要绑定 `realpath` 所在目录（rule 的 `bind_dirs` 已这样做）。

**建议（会改变 wgs_skim 调用方式，需实测后改，旧行为放开关后）**：
1. 默认 1.4.2 镜像；conda 1.4.1 时在日志 / 汇总中把 MitoFinder 的 Circularization 标为“不可用”。
2. `-a spades` 前恢复旧流程的预筛（≥ 1 kb、覆盖度 ≥ 5 × 中位数），提速并让 MitoFinder 只找到 1 条 contig。
3. `--megahit` 模式先抽样到 5–10 M 对 reads（作者 #13 建议；当前把全部 clean reads 交给 MEGAHIT）。
4. Singularity 调用加 `--cleanenv`（#67）。

## 8. 资料来源与查阅日期（2026-10-29）

- GitHub README（v1.4.2）、Releases v1.4.1 / v1.4.2、compare v1.4.1...v1.4.2：https://github.com/RemiAllio/MitoFinder
- Issues #10 #13 #17 #19 #23 #26 #36 #42 #43 #44 #47 #48 #50 #51 #54 #55 #58 #60 #61 #62 #63 #66 #67（gh api）
- anaconda.org bioconda/mitofinder（1.4、1.4.1；依赖列表）；mamba dry-run（2026-10-29）
- 本地：`~/mybiosofts/MitoFinder/mitofinder`（1.4.2 源码，环化逻辑第 905–990 行）、`mitofinder_v1.4.2.sif`（`singularity inspect`）
- 本项目：`workflow/rules/mitofinder.smk`、`envs/mitofinder.yaml`、`config/config.yaml`；`Temp_scripts/mitofinder-batch/`；
  MJ20260820041 `results/10.mitofinder/`（CHY006–008 日志、`mitofinder_summary.tsv`）、`ISSUES.md` #36 #50、
  `Temp_scripts_dev/wgs_skim_pipeline/scripts/spades_mitofinder.py`、TROUBLESHOOTING.md；Mitom `docs/benchmark_six_tools_2026-10-22.md`

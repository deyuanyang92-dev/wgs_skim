# MEANGS 排错参考

## 1. 简介与用途

MEANGS（MitoDNA Extending Assembler from NGS data）是“无种子”的动物线粒体组装器：先用 nhmmer + 各门线粒体
profile HMM 从 reads 中挑出线粒体 reads，组装成初始 scaffold（自发现种子），再（`--deepin`）用这些种子在全部输入
reads 上做 trie 搜索延伸。

**在 wgs_skim 中的作用**：只负责给 GetOrganelle 提供**本样品自己的种子**（`mito/<s>/meangs/<s>.seed.fasta`）。
GetOrganelle 内置 animal_mt 种子对多毛类等类群很远，常出现 "Too few seed reads"（旧流程 README 已记录）；MEANGS
失败或种子为空时，GetOrganelle 退回 animal_mt 种子，流程不中断。MEANGS 自己的组装结果不进入最终结果
（README 比较表：MEANGS 4–11 条 contig、只覆盖 60–80 %）。

- GitHub：https://github.com/YanCCscu/MEANGS
- conda：`conda install -c yccscucib meangs`（不在 bioconda）；Docker：`bioinfodocker/meangs`
- 论文：Song M-H, Yan C, Li J-T (2022) MEANGS: an efficient seed-free tool for de novo assembling animal
  mitochondrial genome using whole genome NGS data. *Briefings in Bioinformatics* 23(1): bbab538,
  doi:10.1093/bib/bbab538（DOI 经 Crossref 检索得到，作者 / 卷期按记忆，未逐字核对）
- 鱼类评估：Xu et al. 2024, *Brief Bioinform* bbae041（doi:10.1093/bib/bbae041）

## 2. 版本

| | 版本 | 说明 |
|---|---|---|
| wgs_skim 固定 | **1.2.1**（`workflow/envs/meangs.yaml`，频道 yccscucib，python 3.9） | 2023-04-04；"fix some bugs in filtering fasta files" |
| GitHub 最新 | **1.3.1**（2026-03-17） | conda 频道 yccscucib 截至 2026-10-29 仍只有 1.2.1 |

1.2.1 → 1.3.1 的全部改动（GitHub compare，README / meangs.py / ReadsIsolation.py / detercirc.py）：

1. 支持单端数据（只给 `-1`）；`-1` 也接受 `.fa/.fas`。
2. `QC_Convert` 的 `trimfq -q` 改为参数（默认值不变，0.01）。
3. `detercirc.py`（`--clip` 用）修复格式串错误：多个候选切点时原来会崩溃。
4. 版本号字符串从 "1.0" 改为 "1.3.1"（1.2.1 的 `meangs.py` 里写的仍是 VERSION="1.0"，日志里看到 1.0 不代表装错）。

**是否升级：不必。** wgs_skim 只用双端数据、不用 `--clip`，两项改动都碰不到；而且 1.3.1 要从源码装（需 gcc ≥ 8.3.1、
pcre ≥ 8.41，见 issue #10 / #4），conda 可复现性反而下降。若以后要支持单端数据再考虑。

## 3. wgs_skim 中的调用方式

`workflow/rules/mito.smk` rule `meangs`（8 线程，mem_mb 12000）：

```bash
cd mito/<s>/meangs
meangs.py --silence -1 R1.clean.fq.gz -2 R2.clean.fq.gz -o meangs -t 8 -i 350 \
    --species_class A-worms -n 2000000 --deepin > logs/meangs/<s>.log 2>&1 || echo "[wgs_skim] MEANGS failed ..."
seed=$(find . -name 'scaffold_seeds.fa*' | head -n 1)     # 有且非空 -> <s>.seed.fasta，否则写空文件
rm -rf meangs
```

| 参数 | 含义 | 配置键 |
|---|---|---|
| `-1/-2` | fastp 之后的 clean reads（软链接先 `readlink -f`） | — |
| `-o meangs` | 输出前缀兼目录名 | — |
| `-t` | 线程（nhmmer、ReadsIsolation） | rule threads |
| `-i 350` | 文库插入片段长度，ReadsIsolation 用来配对 | `insert_size` |
| `--species_class` | 选哪一门的 HMM（A-worms = 环节动物，N-worms = 纽形动物） | `meangs_class`（空 = 按 `taxon` 预设） |
| `-n 2000000` | 每个文件只取 QC 后的**前** N 条 reads（`head`，不是随机抽样） | `meangs_reads` |
| `--deepin` | 生成 `scaffold_seeds.fas` 并在全部输入 reads 上延伸 | 固定 |
| `--silence` | 屏幕输出写进日志 | 固定 |

开关：`meangs: false` 不跑 MEANGS（`--steps getorganelle` 同样不用 MEANGS 种子）。

## 4. 重要参数说明

| 参数 | 默认 | 说明 / 建议 |
|---|---|---|
| `-n` | 0（全部） | 作者建议一定要设，以控制时间和内存；结果随 -n 不单调变化（#18）。种子用途 2 M 足够：MJ20260820041 22 个样品全部得到种子（1–21 条、7.5–13.6 kb） |
| `-i` | 无 | 插入片段；不知道时问测序公司（#14）。fastp JSON 里有 `insert_size.peak`，可据此自动给（未实现） |
| `--species_class` | 无 | 门外类群：先不加 `--deepin` 试几个 class，有 `*_hmmout_tbl` 命中再跑 deepin（#16） |
| `--deepin` | 关 | 关闭时只有初始 scaffold，不写 `scaffold_seeds.fas`（wgs_skim 需要它） |
| `--skipextend` | 关 | `--deepin` 时跳过全量延伸；`scaffold_seeds.fas` 在延伸**之前**就已写出（见第 6 节） |
| `--clip` | 关 | 按重复 k-mer 找环切点；作者建议首次运行不要用（#9） |
| `--keepIntMed` | 关 | 保留 `*.input.fas` 等中间文件（排错时用） |

## 5. 已知问题与 issues

| 现象 / 报错 | 原因 | 解决 | 链接 |
|---|---|---|---|
| `Killed300` / `Killed400`，停在 `tools/assembler ... -s scaffold_seeds.fas` | deepin 延伸阶段内存不足 | 加内存、减少并行、减小 `-n` | [#20](https://github.com/YanCCscu/MEANGS/issues/20) |
| 结果多条 scaffold、同一基因出现多次 | 深度不足 / 片段化 | 加大 `-n` 或不用 `-n`；不加 `--deepin` 可得到 13 个 PCG；用 MEANGS 结果作种子交给种子型软件 | [#13](https://github.com/YanCCscu/MEANGS/issues/13) |
| 线粒体偏短（13 kb vs 15.6 kb），PCG 完整 | 控制区（D-loop）重复组装不出 | 选一段非末端序列作种子，用其他种子型软件重组 | [#19](https://github.com/YanCCscu/MEANGS/issues/19) |
| 增加 reads 后 scaffold 数、基因重复数反而波动 | 算法对数据量敏感；完整结果只认 `*_deep_detected_mito.fas` | 多试几个 -n；注释用 MITOS2 | [#18](https://github.com/YanCCscu/MEANGS/issues/18) |
| 输出没有 circular 标记 | MEANGS 不判断环；`--clip` 只找切点 | 用其他方法（wgs_skim topology / circularize） | [#9](https://github.com/YanCCscu/MEANGS/issues/9) |
| 深度模式不组装线粒体 reads | 早期版本问题，作者要求屏幕日志 | 不用 `--silence` 看完整输出 | [#1](https://github.com/YanCCscu/MEANGS/issues/1) |
| `libpcre.so.1: cannot open shared object file` | 源码编译的 assembler 找不到 pcre | 用 conda 包，或装 pcre 8.41 并设 LD_LIBRARY_PATH | [#10](https://github.com/YanCCscu/MEANGS/issues/10) |
| 编译 / 系统要求问题 | gcc ≥ 8.3.1、pcre ≥ 8.41 | 用 conda / Docker | [#4](https://github.com/YanCCscu/MEANGS/issues/4)、[#22](https://github.com/YanCCscu/MEANGS/issues/22) |
| `-i` 不知道填多少 | — | 文库插入片段，问测序公司 | [#14](https://github.com/YanCCscu/MEANGS/issues/14) |
| 列表外的门 / 非动物 | 只有 9 个 class 的 HMM | 先无 deepin 试各 class | [#16](https://github.com/YanCCscu/MEANGS/issues/16) |

## 6. 本项目中遇到的问题

1. **种子路径错（ISSUES #45，已修复）**：`scaffold_seeds.fas` 写在**当前工作目录**（源码 `mitoSeeds='scaffold_seeds.fas'`
   是相对路径），不在 `-o` 目录里。规则原来只在 `meangs/` 下找 → 种子为空 → GetOrganelle 用 animal_mt 种子
   （种子 reads 只有 25 KB）。现在 `cd` 到规则目录并 `find .`。
2. **`rm -rf meangs` 也删掉了 MEANGS 自己的日志**：`--silence` 时详细输出在 `meangs/` 里；wgs_skim 把 stdout 重定向到
   `logs/meangs/<s>.log`（实测 4.3 MB，含 assembler 进度），一般够用；要看中间文件时手工加 `--keepIntMed` 并去掉 `rm -rf`。
3. **资源**：HNP1sub（0.6 Gb）338–432 s、峰值 7.6 GB；e2e2 两个小样品 239–289 s、5.8 GB；旧流程 2 M reads 每样品约 7 min。
   MEANGS 比 GetOrganelle 第 0 轮还慢（254–276 s）。
4. **`--deepin` 的全量延伸对 wgs_skim 是白算（源码确认，未实测提速）**：1.2.1 的流程是 初始组装 → `FindMitoScaf` →
   `scaffold2seed.py` 写 `scaffold_seeds.fas`（去掉含 N 的片段、≥150 bp）→ 再在全部输入 reads 上延伸
   （`*_deep_scaffolds.fa`）→ 规则只拷贝 `scaffold_seeds.fas`、删除其余。延伸是最占内存的一步（#20 的 Killed 就在这里）。
   两个可选改进（需在 HNP1sub 和 2–3 个 MJ 样品上对比 GetOrganelle 结果后再改，旧行为保留在开关后）：
   a) 加 `--skipextend`：种子不变，省去延伸（之后 `FindMitoScaf` 找不到 deep 文件可能报错，规则的 `|| echo` 能兜住）；
   b) 反过来用 `*_deep_detected_mito.fas`（更长）作种子。
5. 种子质量：旧批次 22 个样品种子 1–21 条、7.5–13.6 kb；CWSC009 有 21 条碎片，仍能带起 GetOrganelle。

## 7. 排错清单

- 日志：`logs/meangs/<s>.log`；基准：`benchmarks/meangs/<s>.tsv`；结果：`mito/<s>/meangs/<s>.seed.fasta`（0 字节 = 无种子）。
- 关键字：
  - `there are N sequences pairs` —— HMM 挑出的线粒体 read 对数（HNP1sub：1776 对）；N 很小 → class 不对或深度太低。
  - `[wgs_skim] MEANGS failed (rc=...)` —— MEANGS 非零退出，GetOrganelle 会用 animal_mt 种子。
  - `Killed` / `Killed300` —— 内存；`Error occured when running command:` 后面是失败的子命令。
  - `gene [] not exists in assembly scaffold` —— FindMitoScaf 的提示，HNP1sub 正常运行时也出现，不是错误。
  - `libpcre.so.1` —— assembler 动态库。
- 判断种子是否起作用：看 GetOrganelle 日志 `Seed reads made: ... (N bytes)`，几十 KB 说明种子太远（MEANGS 种子
  常为几百 KB，如 pisione 794 KB）。
- 复现：`conda run -p <meangs env> meangs.py -1 .. -2 .. -o t -t 8 -i 350 --species_class A-worms -n 2000000 --deepin --keepIntMed`，
  不加 `--silence`，在空目录里运行（种子写在当前目录）。

## 8. 资料来源与查阅日期（2026-10-29）

- GitHub README、Releases（1.0.1 / 1.2.1 / 1.3.1）、compare 1.2.1...1.3.1：https://github.com/YanCCscu/MEANGS
- Issues #1 #4 #9 #10 #13 #14 #16 #18 #19 #20 #22（gh api 读取）
- anaconda.org：yccscucib/meangs 只有 1.2.1；bioconda 无此包
- 本地源码：wgs_skim conda 环境 / `~/miniforge3/envs/meangs/bin/meangs.py`（1.2.1）、`tools/scaffold2seed.py`
- 本项目：`workflow/rules/mito.smk`、`config/config.yaml`、`docs/TESTING.md`、`README.md`、
  MJ20260820041 `ISSUES.md` #45、`results/02.organelle/01.MEANGS_results/`、`wgs_skim_test/*/benchmarks/meangs/`

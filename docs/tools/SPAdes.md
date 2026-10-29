# SPAdes / rnaSPAdes（wgs_skim 排错参考）

## 1. 简介与用途

SPAdes 是 de Bruijn 图短读长组装器，同一个 `spades.py` 通过模式开关提供 isolate / sc / meta / plasmid / metaviral /
metaplasmid / bio / sewage / corona / rna（rnaSPAdes）/ rnaviral。在 wgs_skim 中：

- `rule spades`（`workflow/rules/assembly.smk`，"genome" 模式）：全基因组 skim 组装 → `assembly/<s>/<s>.scaffolds.fasta`，
  供 BUSCO、核标记拷贝数（rDNA / 组蛋白）、`mt_from_spades` 线粒体后备、MitoFinder（`mitofinder_assembly: spades`）使用。
- `rule spades_mode`（`workflow/rules/transcriptome.smk`）：`spades_modes` 中除 genome 以外的模式 → `<s>.<mode>.fasta`；
  `--steps rna` = rnaSPAdes + Trinity + BUSCO（`-m transcriptome`）。
- GetOrganelle 内部调用 SPAdes：`workflow/envs/getorganelle.yaml` 固定 `spades =3.15.5`（与 4.2.0 环境分开）。
- Nextflow：`nextflow/modules/assembly.nf` 的 `SPADES`（只有 genome 模式；Nextflow 版不接受其他 SPAdes 模式 / Trinity）。

链接：
- GitHub：https://github.com/ablab/spades ；手册：https://ablab.github.io/spades/ （运行选项 running.html，rnaSPAdes rna.html）
- 论文：Bankevich et al. 2012 *J Comput Biol* 19:455，doi:10.1089/cmb.2012.0021；Prjibelski et al. 2020 *Curr Protoc
  Bioinformatics* 70:e102，doi:10.1002/cpbi.102；rnaSPAdes：Bushmanova et al. 2019 *GigaScience* 8:giz100，doi:10.1093/gigascience/giz100

## 2. 版本

| | 版本 | 日期 |
|---|---|---|
| wgs_skim 固定（组装 / 各模式） | `spades =4.2.0`（`workflow/envs/spades.yaml`） | 2025-05-06 |
| GetOrganelle 内部 | `spades =3.15.5`（`workflow/envs/getorganelle.yaml`） | 2022-07-16 |
| 最新 | **4.3.0** | 2026-06-14（另有 2025-08 的 continuous 预发布） |

4.2.0 本身：新增 hpcSPAdes（集群）、修 BinSPreader GFA 1.2、Python 小修。4.3.0 release 原文要点：
- 新的省内存模式 `--frugal`（手册："sacrifice some assembly result contiguity in exchange for reduced memory consumption…
  affects the assembly results in an unpredictable way"）
- "Improved overall RAM consumption for all pipelines"
- 正式的错误码体系（清晰的退出码和错误信息），修复错误处理；hpcSPAdes / grid engine 修复；提高若干依赖最低版本。
- `spades.py -h` 对比（本机 `envs/spades` 4.2.0 vs `envs/wgs_spades` 4.3.0）：唯一新增选项是 `--frugal`。

**是否升级：值得在测试后升到 4.3.0，但不紧急。** 理由：我们的主要失败模式是内存（README Troubleshooting），4.3.0 全流程
降内存 + `--frugal` 可作为最后一次重试的兜底；清晰的退出码便于判断是否该重试。风险：结果连续性可能变化，需要用 HNP1sub
与 MJ20260820041 的几个样品比较 N50 / BUSCO 后再改 pin。bioconda 已有 `spades-4.3.0-hde4eca7_2`（本机缓存）。
GetOrganelle 的 3.15.5 不要跟着改（GetOrganelle 1.7.7.1 与 4.x 的兼容性未验证）。

## 3. wgs_skim 中的调用方式

genome 模式（`rule spades`，`threads: 32`）：

```bash
mem=$(( mem_mb / 1024 ))                       # mem_mb = max(spades_min_mem_gb, clean_Gb × spades_gib_per_gb) × 1.5^(attempt-1)
if [ -f assembly/<s>/spades/params.txt ]; then # 已有运行 → 续跑
    spades.py -o assembly/<s>/spades --restart-from last -t {threads} -m $mem
else
    spades.py -1 R1 -2 R2 -o assembly/<s>/spades -t {threads} -m $mem --only-assembler
fi
cp spades/scaffolds.fasta <s>.scaffolds.fasta
rm -rf spades/K* spades/tmp spades/misc spades/before_rr.fasta   # 中间文件
```

其他模式（`rule spades_mode`，`threads: 16`）：`spades.py --<mode> -1 R1 -2 R2 -o spades_<mode> -t -m [-k ...] [--careful] [--only-assembler]`，
同样的 `params.txt` 续跑逻辑；结果文件 rna/rnaviral = `transcripts.fasta`，bio = `gene_clusters.fasta`，其余 `scaffolds.fasta`；
结束后删除 `K* tmp misc corrected`。

| 配置键（`config/config.yaml`） | 默认 | 作用 |
|---|---|---|
| `spades_gib_per_gb` / `spades_min_mem_gb` | 6.5 / 32 | genome 模式 `-m` 估计（注释：环节动物实测 5.3–6.1 GiB/Gb） |
| `spades_min_clean_gb` | 0 | clean 数据少于此值不组装（输出空文件） |
| `spades_modes` | `genome` | 逗号分隔的模式列表 |
| `spades_kmers` | `""` | 非 genome 模式的 `-k`（空 = SPAdes 自选） |
| `rnaspades_kmers` | `21,33,55` | rna / rnaviral 的 `-k`；`auto` = SPAdes 自选（**见第 5 节，建议改默认**） |
| `spades_careful` | false | `--careful`（仅 isolate / sc / plasmid） |
| `spades_only_assembler` | false | 非 genome 模式的 `--only-assembler`（genome 模式总是开） |
| `spades_mode_gib_per_gb` | 5 | 非 genome 模式内存估计（≥16 GB） |
| 线程 | — | `--set-threads spades=48` 或 profile；profile `retries` 决定 OOM 后 ×1.5 重试 |

## 4. 重要参数说明（手册 running.html / `spades.py -h` 4.2.0）

| 参数 | 说明 |
|---|---|
| `-m <int>` | 内存上限（Gb），默认 250。"SPAdes terminates if it reaches this limit"；实际用量低于此值。设成机器真实可用量 |
| `-t <int>` | 线程，默认 16 |
| `-k` | 奇数、<128、升序；默认 auto：读长 ≥150 时 `21,33,55,77`，250 bp 用到 127（源码 `options_storage.py`） |
| `--only-assembler` | 跳过 BayesHammer 纠错（省时间省内存，wgs_skim genome 模式固定使用） |
| `--restart-from <cp>` | `ec` / `as` / `k<int>` / `mc` / `last`；可改选项，**但不能改基本选项、任何输入数据选项、`--only-error-correction`、`--only-assembler`** |
| `--continue` | 从最后检查点继续，只允许 `-o` |
| `--checkpoints last/all` | 非必需；会多占时间和"significant amount of disk space" |
| `--tmp-dir` | 纠错临时文件目录（默认 `<out>/tmp`）；网络盘上卡住时移到本地盘（#1444 维护者建议） |
| `--careful` | 只推荐小基因组；"strongly recommend not to use it for large and medium-size eukaryotic genomes" |
| `--isolate` | 推荐用于高覆盖度 isolate / 多细胞数据（动物低覆盖 skim 不适用，维持默认模式） |
| `--meta` | 只支持一个双端短读长文库，不支持 careful 和 `--cov-cutoff` |
| `--rna` + `--ss rf/fr` | rnaSPAdes；链特异文库需给 `--ss`；不能与 `--careful`、`--only-error-correction` 同用 |
| `--frugal`（4.3.0） | 省内存，结果不可预测地变差 |
| `--gfa11` | 输出旧版 GFA 1.1（4.0 起默认 GFA 1.2，旧工具读不了时用） |

rnaSPAdes k 值：手册 rna.html："By default rnaSPAdes uses 2 k-mer sizes, which are automatically detected using read length.
We recommend not to change this parameter because smaller k-mer sizes typically result in multiple chimeric (misassembled)
transcripts." 4.2.0 源码 `rna_k_values()`：上 k = 读长/2−1，下 k = min(max(读长/3, 29), 55)，即 150 bp 读长 → **k = 49,73**。
输出：`transcripts.fasta`（推荐）、`hard_filtered_transcripts.fasta`、`soft_filtered_transcripts.fasta`。

## 5. 已知问题与 issues

| 现象 / 报错 | 原因 | 解决办法 | 链接 |
|---|---|---|---|
| `<jemalloc>: Error in malloc(): out of memory` | 操作系统未能满足分配（不是 SPAdes 能恢复的）；`-m` 只是上限 | 提高 `-m` 与实际可用内存；减少同机任务；wgs_skim 重试 ×1.5 | https://github.com/ablab/spades/issues/494 、 https://github.com/ablab/spades/issues/19 |
| 2 TB 机器上仍 OOM | 实际 *空闲* 内存不足 / 其他进程占用 | 检查空闲内存而非总内存 | https://github.com/ablab/spades/issues/871 |
| 实际用量超过 `-m` 被 HPC 杀（hard limit） | `-m` 是 SPAdes 自身的上限，调度器的限制需留余量 | 调度内存 > `-m` | https://github.com/ablab/spades/issues/1536 |
| `finished abnormally, OS return value: -9` / `err code -9` | 进程被 SIGKILL（OOM killer / 调度器） | 同上；查 `dmesg`/调度日志 | https://github.com/ablab/spades/issues/76 、 https://github.com/ablab/spades/issues/1338 （rnaSPAdes） |
| `OS return value: -11`（段错误，4.2.0 报告） | 未定（open） | 换 k 值 / 版本，附 spades.log 报 issue | https://github.com/ablab/spades/issues/1531 、 https://github.com/ablab/spades/issues/1521 |
| `Assertion nucls.size() >= index.k() failed`，restart 后同样 | 维护者："likely indicates some memory corruption" | 从头重跑 / 换机器 | https://github.com/ablab/spades/issues/1540 |
| 卡在 "Collecting perfect loops"（4.2.0，K33，7 天无输出） | 未解决（open） | 记录 K 阶段，考虑换 k 或 `--restart-from k<下一个>` | https://github.com/ablab/spades/issues/1539 |
| 日志无报错但长时间不动（纠错阶段） | 维护者："Likely the problem is around I/O"；临时目录在网络盘 | `--tmp-dir` 放本地盘 | https://github.com/ablab/spades/issues/1444 |
| `Cannot create temporary file` / `OS return value: 122` | 磁盘满 / 配额 / 目录不可写 | 腾空间；`--tmp-dir` | https://github.com/ablab/spades/issues/998 |
| 换服务器后 `--restart-from last` 找错 reads 路径 | 续跑使用 params 中记录的原始路径 | 保持 reads 绝对路径不变 | https://github.com/ablab/spades/issues/666 |
| 4.2.0 只看前 10 000 条 reads 估计读长，按长度排序的数据 k 被错误剔除 | `get_max_reads_length()` 抽样（已关闭） | 读长混杂时手动给 `-k` | https://github.com/ablab/spades/issues/1591 |
| 需要 ~1.1 TB 内存的超大数据 | 数据量 | 先单独纠错再组装、降数据量（讨论中） | https://github.com/ablab/spades/issues/1373 |
| 3.15.3：`module 'collections' has no attribute 'Hashable'` | 旧 pyyaml 不兼容 Python 3.10 | 用 3.15.4+（我们 GetOrganelle 环境是 3.15.5 + python 3.10，没问题） | https://github.com/ablab/spades/issues/873 |
| Python 3.12 `SyntaxWarning: invalid escape sequence` | 警告，非错误 | 忽略 | https://github.com/ablab/spades/issues/1457 |

## 6. 本项目中遇到的问题与代码观察

- README Troubleshooting：`<jemalloc>: ... out of memory` 或被杀 → 内存估计偏低，profile `retries: 1` 给 1.5 倍并从最后检查点续跑；
  `clean data ... < spades_min_clean_gb` → 不组装。TESTING §3.4：中间文件清理 2.1 GB → 817 MB；0.6 Gb 样品峰值 6.9 GB。
- TESTING §9.2：rnaSPAdes（`--rna -k 21,33,55`）两个 0.1 Gb RNA-seq 样品 87–92 s、1.2–1.3 GB，正常。
- **代码观察 1（未实测，按手册推断）**：续跑条件只是 `params.txt` 存在。手册规定 `--restart-from` 不能改输入数据；所以若 fastp
  重跑（clean reads 变了）或改了 `spades_kmers` / `rnaspades_kmers` / `spades_careful`，Snakemake 重跑该规则时，SPAdes 会沿用旧
  reads 和旧参数（`-k` 等新值不传给续跑命令）。建议：只在 `attempt > 1` 时续跑，否则先 `rm -rf` 输出目录；或在目录里记录输入/参数指纹。
- **代码观察 2（未实测）**：成功后删除了 `K*`，但 `params.txt` 仍在；`--forcerun spades` 时会走 `--restart-from last`，行为未验证。
- **代码观察 3**：`rnaspades_kmers: "21,33,55"` 与手册建议（保持自动、小 k 易嵌合）相反；150 bp 读长时自动值为 49,73。建议默认改 `auto`，
  旧值保留为可选（来自用户旧 rnaspades snakefile）。
- 磁盘：SPAdes 中间文件约 2 GB/Gb clean 数据（规则注释）；20 Gb 样品峰值可能 ~40 GB，结果目录放大盘。

## 7. 排错清单

1. `logs/spades/<s>.log`（或 `logs/spades_mode/<s>.<mode>.log`）末尾的 `== Error ==` 段：`system call for: [... spades-core / spades-hammer ...]
   finished abnormally, OS return value: N`。N=-9 被杀（内存）；-11 段错误；-6 断言；-7 SIGBUS（推断：多与磁盘 / 内存映射有关，issues 中无维护者定论）；122 磁盘配额。
2. `assembly/<s>/spades/spades.log`（完整日志）、`params.txt`（实际命令、版本）、`warnings.log`。
3. 关键字：`jemalloc`、`out of memory`、`Memory limit`、`Cannot create temporary file`、`No space left`、`Assertion`、
   `K values to be used`（rnaSPAdes 实际 k）、`Default k-mer sizes were set to`。
4. `benchmarks/spades/<s>.tsv` 的 `max_rss` 与 `-m`（log 第一行附近 `Memory limit`）对比；接近 → 提高 `spades_gib_per_gb`。
5. 续跑后结果异常 → 检查 `params.txt` 中的 reads 路径 / k 是否是当前的（见第 6 节观察 1）。
6. 磁盘：`du -sh assembly/<s>/spades`；`df -h` 结果盘与 tmp。

## 8. 资料来源与查阅日期（2026-10-29）

- Releases：https://github.com/ablab/spades/releases （`gh api` 读取 4.2.0 / 4.3.0 说明）
- 手册：https://ablab.github.io/spades/running.html 、https://ablab.github.io/spades/rna.html
- 本机 `spades.py -h`（4.2.0 `~/miniforge3/envs/spades`，4.3.0 `~/miniforge3/envs/wgs_spades`）与 4.2.0 源码
  `share/spades/spades_pipeline/options_storage.py`、`stages/spades_stage.py`。
- 上表 issues（2026-10-29 查阅）；wgs_skim `workflow/rules/assembly.smk`、`transcriptome.smk`、`common.smk`、`config/config.yaml`、
  `docs/TESTING.md`、`README.md`。

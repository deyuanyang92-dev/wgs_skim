# 开发说明：模块结构与扩展方法

wgs_skim 功能多，所以代码按"一处实现、到处复用"组织。新增功能前先看这里：要用的函数很可能已经有了。

## 结构

```
src/wgs_skim/                Python 包（命令行与 Python 接口）
  cli.py                     命令行：只解析参数，调用下面的模块
  paths.py                   位置：workflow/、~/.wgs_skim、软件环境目录、机器内存
  engine.py                  引擎：由一组工作流选项生成并运行 Snakemake / Nextflow 命令
  tools.py                   工具登记表 TOOLS：每个工具的输入、步骤、专有选项及其配置（wgs_skim <工具>）
  api.py                     Python 接口：wgs_skim.run(...)、wgs_skim.tool(...)、wgs_skim.command(...)
  lib.py                     workflow/scripts 作为库：from wgs_skim.lib import seqtools
workflow/Snakefile           include rules/*.smk
workflow/rules/*.smk         每类工具一个规则文件：输入、输出、资源（threads、mem_mb）、conda 环境、命令
workflow/scripts/*.py        真正的计算；每个脚本同时是：可 import 的模块 + 命令行工具 和 / 或 Snakemake 脚本
workflow/envs/*.yaml         每个工具一个 conda 环境，版本固定（见 docs/tools/<工具>.md 的理由）
nextflow/                    同一流程的 Nextflow 版本，调用同一批 scripts（bin/smk_shim.py）
```

## 脚本的约定（workflow/scripts）

1. **import 时不执行任何操作**：工作放在 `main()`（命令行）或 `main(sm)`（Snakemake `script:`）中，文件末尾：

   ```python
   if __name__ == "__main__" and "snakemake" in globals():
       main(snakemake)  # noqa: F821
   ```

2. **import 同目录的公共模块**，开头加：

   ```python
   sys.path.insert(0, str(Path(__file__).resolve().parent))   # the shared modules next to this script
   if "snakemake" in globals():                                  # a Snakemake script: run from a temporary copy
       sys.path.insert(0, snakemake.scriptdir)  # noqa: F821
   from seqtools import read_fasta  # noqa: E402
   ```

3. **不要再写一份已有的函数**。公共函数：

| 模块 | 函数 | 用途 |
|---|---|---|
| `seqtools` | `iter_fasta(path, first_word, upper)`、`read_fasta`、`parse_fasta(text)`、`write_fasta(recs, 文件或路径, width)`、`rc`、`fasta_lengths`、`n50`、`is_circular(header)`、`genbank_seqs`、`find_circle` | 序列：FASTA 读写（含 .gz、逐条读取大组装）、反向互补、长度与 N50、GenBank 序列（不需 Biopython）、首尾重叠 |
| `tables` | `write_tsv(path, rows, cols, append)`、`read_tsv` | 所有结果表的读写 |
| `runtools` | `run(cmd, log, env, check, cwd)`、`sh`、`output`、`tool_env`、`blastn(query, target, fields, evalue, threads, db)`、`aligned_bp` | 运行外部程序（日志、环境）；BLAST（大目标自动建临时库） |
| `samples` | `find_pairs`、`read_sheet`、`sheet_columns`、`fasta_samples`、`split_fasta`、`by_sample` | 样品识别：reads 目录、样品表、FASTA 输入（含按序列名拆分）、按样品目录取文件 |
| `collect_fasta` | `getorganelle_cmd`、`result(path)`、`result_fasta`、`best`、`clean`、`renamed` | GetOrganelle：命令、结果选择（含组装图）、改名、临时文件清理 |
| `mitofinder_run` | `run_mitofinder`、`collect`、`mitofinder_exe`、`first_pairs` | 一次 MitoFinder 运行（组装或 reads；Singularity 或 bioconda），只用标准库 |
| `mitofinder_prefilter` | `select`、`ref_hits`、`weighted_median` | MitoFinder 输入：深覆盖 scaffold + 参考序列命中的 scaffold |
| `circularize` | `circularize(recs, spades)`、`rc` | 无 reads 环化（组装图 / 重叠 / 拼接） |
| `blast_taxonomy` | `blast_topn(query, out_tsv, db, topn, threads, email)`、`mt_genes`、`nr_regions` | BLAST + GenBank 物种分类（blast2metadata 的流程） |
| `mitofinder_summary` | `sample_result` | MitoFinder 结果汇总 |
| `mitofinder_check` | `check(rows, qids, min_pident)`、`family` | 基因级污染判断 |
| `nuclear_markers` | `getorganelle`、`annotate`、`compare`、`is_complete`、`layout` | 核多拷贝基因（rDNA、组蛋白） |
| `assembly_check` | `check`、`read_gff` | QUAST 式结构比较 |
| `summary` | `genomescope`、`genome_size_check`、`busco_line`、`n50` | 汇总表 |

4. 外部程序用 `runtools.run`（日志、退出码）与 `runtools.blastn`，表格用 `tables`，FASTA 用 `seqtools`，
   GetOrganelle 用 `collect_fasta`；不要在脚本里再写 `subprocess.run`、`csv.DictWriter`、FASTA 解析。
   在旧工具环境（如 Python 2.7 的 bioconda MitoFinder）里运行的脚本只用标准库，并在规则中用 `{PY}`（Snakemake 的 Python 3）启动。

5. 每个新函数在 `tests/test_scripts.py` 有单元测试（秒级，不需要外部软件）；断言用真实数据中观察到的情况，注释写明来源。

## 新增一个工具（例：X）

1. `workflow/envs/x.yaml`：固定版本；在 `docs/tools/X.md` 写明版本理由与已知 issues。
2. `workflow/rules/x.smk`：规则写 `threads`、`resources: mem_mb`（按数据量估计时写函数，重试时增加）、`log`、
   `benchmark`、`conda`；在 `workflow/Snakefile` 中 include。计算放在 `workflow/scripts/x.py`（遵守上面的约定）。
3. `workflow/rules/common.smk`：
   - `STEPS_ALL` 加 `x`；需要时在 `COMBOS` 加组合（如 `"spades-x": ["spades", "x"]`）；
   - `final_outputs()` 的 `targets` 加 `"x": [rules.x.output.…]`；
   - 输入可以来自其他工具时，用输入函数选择来源（参考 `rules/mitofinder.smk` 的 `mitofinder_sources`）。
4. `config/config.yaml` + `workflow/schemas/config.schema.yaml`：新配置键（带注释与默认值）；
   `python docs/wiki/tools/config_page.py > docs/wiki/Configuration.md` 重新生成配置页。
5. `src/wgs_skim/tools.py`：在 `TOOLS` 加一个 `Tool(name, inputs, help, step, column, options, configure)`，即得到
   `wgs_skim x` 批量命令与 `wgs_skim.tool("x", ...)`（reads 或 FASTA 输入、-j、--mem-gb、-n 都是共用的）；
   专有选项写在 `options(parser)`，它们怎样变成工作流配置写在 `configure(args, conf, use_fasta)`。
6. 测试：
   - `tests/cli_matrix.sh`：`tool "tool x" <应出现的规则> -- x …`，组合用 `dry "combo …"`（干运行，1 分钟内）；
   - `tests/tools_e2e.sh`：一项真实运行检查（`check NAME "输出文件" 命令`）；
   - 结果写进 `docs/TESTING.md`。
7. 文档：`docs/wiki/Tool-X.md`（做什么、命令、输入、输出、实测时间）、`Tools.md` 表格、`_Sidebar.md`；
   组合写进对应的 `Combo-*.md`；常见错误写进 `Troubleshooting.md`。
8. Nextflow：暂不支持的步骤列在 `engine.py` 的 `NEXTFLOW_MISSING` 中（直接提示，不会只跑一半）。

## 从 Python 使用

```python
import wgs_skim
wgs_skim.run(reads_dir="rawdata/", outdir="res/", steps="spades-busco", cores=32)       # = wgs_skim run ...
wgs_skim.tool("mitofinder", reads_dir="rnaseq/", ref="related.gb", assembler="rnaspades")
from wgs_skim.lib import seqtools, samples                                             # 库函数
pairs = samples.find_pairs("rawdata/")
```
关键字参数就是命令的选项名（`-` 换成 `_`），与命令行完全一致（`tests/test_package.py`）。

## 检查清单（提交前）

```bash
python -m pytest -q tests/                                 # 单元测试
ruff check --select F,E9 --line-length 120 workflow/scripts src
tests/cli_matrix.sh small/ small_results/                  # 每个命令、每个步骤、每个组合（干运行）
tests/tools_e2e.sh small/ small_results/ out/ "" rna_small/   # 每个工具与组合的真实运行
```

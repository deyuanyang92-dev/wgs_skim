# Nextflow（wgs_skim 排错参考）

## 1. 简介与用途

`wgs_skim run --engine nextflow` 用 `nextflow/main.nf`（模块在 `nextflow/modules/*.nf`）跑与 Snakemake 相同的流程（不含其他
SPAdes 模式与 Trinity）。配置默认值直接读 `config/config.yaml`（两个引擎一个来源），脚本经 `nextflow/bin/smk_shim.py`
以 Snakemake 的 `snakemake` 对象接口运行 `workflow/scripts/*.py`。数据库任务（GetOrganelle 库、MitoM 库、BUSCO lineage）用
`storeDir`，存在即跳过（"Stored"）。

链接：
- 代码：https://github.com/nextflow-io/nextflow ；文档：https://docs.seqera.io/nextflow/
- 缓存与续跑：https://docs.seqera.io/nextflow/cache-and-resume ；严格语法：https://docs.seqera.io/nextflow/strict-syntax
- 26.04 迁移说明：https://docs.seqera.io/nextflow/migrations/26-04
- 论文：Di Tommaso et al. 2017 Nat Biotechnol 35:316, doi:10.1038/nbt.3820

## 2. 版本

| | 版本 |
|---|---|
| wgs_skim 测试用 | 25.04.2（`/usr/local/bin/nextflow`，TESTING.md §10）；`manifest.nextflowVersion = '>=24.04.0'`（无上限） |
| 最新稳定 | 26.10.4（2026-10-26）；同时维护 26.04.9、25.10.9 |

主要变化：
- 25.04：`resourceLimits`（wgs_skim 用它把请求封顶到 `--max_cpus` / `--max_memory`）在 24.04 起已有；25.04 起可选用严格语法解析器 v2。
- 25.10：工作流 outputs 转正、typed processes 预览、Java 25。
- **26.04（BREAKING）**：**严格语法解析器默认开启**（退回用 `NXF_SYNTAX_PARSER=v1`）；map 类型输入的哈希改为键 + 值
  → 升级后第一次 `-resume` 会重跑所有带 map 输入的进程；`stageAs` 参数顺序在 typed 语法中变化。

**升级风险已实测**（本机 25.04.2，`NXF_SYNTAX_PARSER=v2 nextflow config nextflow/` 与 `nextflow lint`，2026-10-29）：

```
Error nextflow.config:5:1: Unexpected input: 'new'
│   5 | new groovy.yaml.YamlSlurper().parse(new File("${projectDir}/../config/
ERROR ~ Config parsing failed
```

`nextflow lint`：9 个文件 81 个错误，根源只有两处——
1. `nextflow.config`：顶层语句（`new ...YamlSlurper().parse(...).each { params[k] = v }`、`def stamp = ...`）。严格配置只允许
   赋值、配置块、`includeConfig`；需改成闭包内的动态 include，或在 `main.nf` 中读 YAML。
2. `modules/common.nf:102`：`def out = new StringBuilder(), err = new StringBuilder()`（一行声明两个变量）→ 拆成两行。
   其余 72 个 "`X` is not defined" 都是 common.nf 无法解析造成的连锁错误。

结论：**在 Nextflow ≥ 26.04 上当前代码无法启动**。在修改前，`src/wgs_skim/engine.py` 的 `nextflow()` 应设置 `NXF_SYNTAX_PARSER=v1`
（或 `manifest.nextflowVersion = '>=24.04.0, <26.04'`）；长期应修这两处让 v1/v2 都能解析。`nextflow self-update` 会直接踩到这个问题。

## 3. wgs_skim 中的调用

```bash
cd <outdir>      # .nextflow/（续跑缓存）放在结果目录
nextflow run <repo>/nextflow/main.nf -profile conda -resume -work-dir <outdir>/work \
  -params-file <outdir>/nextflow_params.yaml --max_cpus <-j> --max_memory '<--mem-gb> GB' \
  --python <sys.executable> --conda_prefix <prefix> --<key> <value> ...
```

- `-params-file`：`~/.wgs_skim/config.yaml` 与 `--configfile` 合并（Nextflow 只接受一个）；命令行 `--key` 优先。
- `-resume` 总是开启；`-work-dir` 在结果目录下；`publish_mode = 'link'`（硬链接；work/ 与结果不在同一文件系统时改 `copy`）。
- `nextflow.config`：`process.shell = ['/bin/bash','-euo','pipefail']`；`errorStrategy` 第 1 次失败 retry、之后 `finish`；
  `maxRetries 1`；内存 `16000 × task.attempt` 等；`executor.$local.cpus/memory`；`conda.cacheDir = WGS_SKIM_CONDA_PREFIX
  或 params.conda_prefix`；`conda.createTimeout = '2 h'`；trace / report / timeline 写 `pipeline_info/`。
- `condaEnv(name)`（common.nf）：若 Snakemake 的 conda 前缀里有内容完全相同的 `<hash>_.yaml` 且 `env_setup_done` 存在，
  直接用该环境目录；否则交给 Nextflow 从 YAML 新建。
- profile：`conda`（默认）、`slurm`（`slurm_queue`、`slurm_options`，queueSize 50、submitRateLimit 10/1min）、`test`。

## 4. 重要参数

| 参数 / 设置 | 含义 | wgs_skim |
|---|---|---|
| `-resume [session]` | 复用缓存 | 总是 |
| `-work-dir` | 任务目录 | `<outdir>/work` |
| `-profile conda` | `conda.enabled = true` | 默认 |
| `-stub-run` | 只跑 stub，测 DAG | 测试用 |
| `-dump-hashes [json]` | 打印每个任务的缓存哈希组成 | 续跑失效时用 |
| `process.resourceLimits` | 请求上限封顶 | cpus/memory |
| `cache 'lenient'` | 只用路径 + 大小（忽略时间戳） | 未用；NFS 上可用 |
| `storeDir` | 输出永久存放；存在即跳过，不看哈希 | 数据库任务 |
| `NXF_SYNTAX_PARSER` | v1 / v2 解析器 | 见第 2 节 |

## 5. 已知问题与 issues

| 现象 / 报错 | 原因 | 解决 | 链接 |
|---|---|---|---|
| `Config parsing failed` / `Unexpected input: 'new'`（26.04+） | 严格语法默认开启 | `NXF_SYNTAX_PARSER=v1`；修 config / common.nf | 26.04 迁移说明；本项目实测 |
| 续跑全部重跑 | 任务哈希 = 会话、进程名、容器 / **conda 环境路径**、输入（路径 + mtime + 大小）、脚本、引用的全局变量、`bin/` 脚本、stub | `-dump-hashes` 对比两次运行 | cache-and-resume 文档 |
| 改了 `conda.cacheDir` 后无法续跑 | conda 环境路径计入哈希 | 保持 `--conda-prefix` 不变 | https://github.com/nextflow-io/nextflow/issues/6613 |
| 重试成功的任务每次 `-resume` 都重跑 | 该次运行曾越过残留 work 目录续跑 | 已关闭（2026-10，修复版本未核实）；删掉残留失败目录 | https://github.com/nextflow-io/nextflow/issues/7787 |
| `WARN: Unrecognized config option 'executor.$local.cpus'` | ConfigValidator 误报（选项实际生效） | 忽略 | https://github.com/nextflow-io/nextflow/issues/7551 |
| storeDir 进程失败时 `cp: '.command.out' and '.command.out' are the same file` | storeDir 失败路径 bug | 看真正的 `.command.err` | https://github.com/nextflow-io/nextflow/issues/5774 、/6311 |
| storeDir 仍暂存输入 | 设计如此 | 无 | https://github.com/nextflow-io/nextflow/issues/5468 |
| storeDir 下的数据库版本过期不会更新 | storeDir 只看文件是否存在 | 手动删除 `<db_dir>/busco_downloads/lineages/<lineage>`（记录在案）后重跑 | 文档 |
| NFS 上续跑失效 | 时间戳不一致 | `process.cache = 'lenient'` | cache-and-resume 文档 |
| `WARN: It appears you have never run this project before -- Option -resume is ignored` | 第一次在该目录运行 | 正常 | — |
| 单个任务内存超限被杀（exit 137） | 本地 executor 不强制内存，只用于调度；真实超用由 OOM killer 处理 | 提高该进程 memory；看 trace 的 peak_rss | https://github.com/nextflow-io/nextflow/issues/6797 |

## 6. 本项目中遇到的问题（TESTING.md §10，Nextflow 25.04.2）

- **#17**：BUSCO 6 把暂存链接名 `metazoa_odb12` 当作 lineage 名称 → `-l $(readlink -f lineage)`。
- **#18**：topology 的 `stageAs: 'fasta/*'` 文件 `.name` 已含 `fasta/` → 路径成了 `fasta/fasta/x` → 传暂存路径本身（`names()`）。
  注意 26.04 typed 语法中 `stageAs` 参数顺序改变，升级时要复查这里。
- **#19**：summary.py 不写日志，Nextflow 缺少声明的输出 `logs/summary.log` 而失败 → `smk_shim.py` 预建输出 / 日志。
- 续跑粒度实测：只改 topology 输入写法 → 只有 TOPOLOGY 及下游重跑；只改 `bin/smk_shim.py` → 用 shim 的任务全部重跑（bin/ 计入哈希）。
- `-resume` 两次：38 Cached + 3 Stored，8 s。
- conda 环境复用的副作用（推断）：若先跑 Nextflow（环境建在 cacheDir），之后再 `wgs_skim install`，`condaEnv()` 改为返回
  Snakemake 环境路径 → 哈希改变 → 下一次 `-resume` 全部重跑（同 #6613 机制）。

## 7. 排错清单

1. 终端报错给出任务 work 目录：看 `.command.sh`（实际命令）、`.command.err`、`.command.log`、`.exitcode`、`logs/<step>/<s>.log`。
2. `<outdir>/.nextflow.log`：搜 `ERROR`、`Caused by`、`WARN`；`nextflow log` 列出历次运行与 session id。
3. `pipeline_info/trace_*.tsv`：`status`、`exit`、`peak_rss`、`attempt`（判断是否重试过、内存是否到顶）。
4. 续跑失效：两次运行都加 `-dump-hashes json -log x.log`，`diff` 对应任务的 `cache hash` 组成。
5. 先 `-stub-run` 确认 DAG（数秒），再实跑。
6. 版本：`nextflow -version`；≥ 26.04 时先确认 `NXF_SYNTAX_PARSER`。

## 8. 资料来源（查阅日期 2026-10-29）

- https://github.com/nextflow-io/nextflow/releases （26.10.4、26.04.x、25.10.x 说明）
- https://docs.seqera.io/nextflow/migrations/26-04 ；https://docs.seqera.io/nextflow/strict-syntax ；https://docs.seqera.io/nextflow/cache-and-resume
- Issues：#5468 #5774 #6311 #6613 #6797 #7551 #7787
- 本机实测：`NXF_SYNTAX_PARSER=v2 nextflow config`、`nextflow lint`（Nextflow 25.04.2）
- 本项目：`nextflow/nextflow.config`、`nextflow/modules/common.nf`、`assembly.nf`、`src/wgs_skim/engine.py`（nextflow）、
  `docs/TESTING.md` §4 §10、`../wgs_skim_test/nf_e2e*.log`

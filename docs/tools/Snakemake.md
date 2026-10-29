# Snakemake（wgs_skim 排错参考）

## 1. 简介与用途

wgs_skim 的默认引擎。`wgs_skim run` / `nuclear` / `compare` / `install` 都是拼好命令再调用 Snakemake
（`src/wgs_skim/cli.py`）：`workflow/Snakefile`（主流程）、`workflow/nuclear_batch.smk`、`workflow/compare_batch.smk`。
每个规则声明 `threads` 与 `resources: mem_mb`，Snakemake 按 `--cores` 和 `--resources mem_mb` 调度；软件环境用
`workflow/envs/*.yaml` 由 conda 建在共享的 `--conda-prefix` 下。

链接：
- 代码：https://github.com/snakemake/snakemake ；文档：https://snakemake.readthedocs.io/
- CLI：https://snakemake.readthedocs.io/en/stable/executing/cli.html ；Profiles：同页 "Profiles" 一节
- 论文：Mölder et al. 2021 F1000Research 10:33, doi:10.12688/f1000research.29032.2

## 2. 版本

| | 版本 |
|---|---|
| wgs_skim 固定 | `environment.yml`: `snakemake =8.18`（本机 `~/miniforge3/envs/snakemake` = 8.18.2）；`pyproject.toml`: `snakemake>=8,<10` |
| 最新 | 9.27.1（2026-10-21） |

- 9.0.0（2025-03-14）唯一标为 BREAKING 的是日志重构（LoggerPluginInterface，#3107）；其余为新功能（`--report-after-run`、
  `--replace-workflow-config` 等）和修复（如 "Re-evaluation of free resources" #3399）。
- 9.x 文档：profile 文件默认名为 `profile.v9+.yaml`，`config.yaml` 仍"为向后兼容"而接受——wgs_skim 的
  `workflow/profiles/default/config.yaml` 在 9.x 下应仍可用（未在 9.x 上实跑）。
- 是否升级：**暂不需要**。8.18 在本项目所有测试中工作正常；升到 9.x 前应跑一次 `tests/e2e.sh` + `tests/cli_matrix.sh`
  （日志格式变化可能影响依赖 Snakemake 输出文本的地方）。pyproject 允许 9.x，而 environment.yml 固定 8.18：两处不一致，
  pip 安装到其他环境时可能拿到 9.x。
- 本机陷阱：`/usr/bin/snakemake` 是 **6.15.1**（系统包），不支持 `--sdm`、workflow profile。`wgs_skim` 优先用与自身 Python
  同目录的 snakemake（`snakemake_exe()`），所以 `wgs_skim run` 不受影响；但手动敲 `snakemake ...` 时可能用到 6.15.1。

## 3. wgs_skim 中的调用

```bash
snakemake -s <repo>/workflow/Snakefile --sdm conda --conda-prefix <prefix> --config <KEY=VALUE ...> \
  -d <outdir> --cores <-j> --resources mem_mb=<--mem-gb × 1024> [-n] [-- 之后的参数原样追加]
# install：-d <tmp> --cores 1 --conda-create-envs-only，再建 mitom_db / mitom_nuc_db / getorganelle_db
```

自动加载的 workflow profile `workflow/profiles/default/config.yaml`：

```yaml
software-deployment-method: conda
conda-frontend: conda
printshellcmds: true
rerun-incomplete: true
keep-going: true
retries: 1          # SPAdes 重试时 1.5× 内存（resources 函数用 attempt）
```

本机 8.18.2 实测：`-s wf/Snakefile -d run` 时打印 `Using workflow specific profile wf/profiles/default`，即 profile 相对
**Snakefile 所在目录**查找，`-d` 改工作目录不影响。`--workflow-profile none` 可关闭。

可调：`-j/--cores`、`--mem-gb`（默认机器内存 90 %）、`--conda-prefix`（默认共享目录）、`-- --set-threads spades=48`、
`-- --rerun-triggers mtime`、`-- --unlock`、`-- --forcerun <rule>`。

## 4. 重要参数

| 参数 / 键 | 含义 | wgs_skim 用法 |
|---|---|---|
| `--resources mem_mb=N` | 全局内存池；规则 `mem_mb` 之和不超过 N | `--mem-gb` |
| 规则 `resources: mem_mb` | 每个任务声明；可为函数 `(wildcards, attempt)` | SPAdes 按数据量估算 × 1.5^attempt |
| `--default-resources` | 未声明时默认 mem_mb = 2×输入大小，限 1000–8000 | 所有规则都显式声明 |
| `--retries` | 失败重试次数（默认 0） | profile 中 1 |
| `--keep-going` | 一个任务失败其他独立任务继续 | 开 |
| `--rerun-incomplete` | 重跑被标记为不完整的输出 | 开（所以一般不会看到 IncompleteFilesException） |
| `--rerun-triggers` | 默认 mtime、params、input、software-env、code | 想忽略参数变化时 `mtime` |
| `--sdm conda` / `--conda-prefix` | 环境部署方式 / 环境目录 | 共享目录，多次运行复用 |
| `--conda-create-envs-only` | 只建环境 | `wgs_skim install` |
| `--unlock` | 解除工作目录锁 | 上次被中断时 |
| `-n -p` | 干运行并显示命令 | 排错 |

## 5. 已知问题与 issues

| 现象 / 报错 | 原因 | 解决 | 链接 |
|---|---|---|---|
| `LockException: Directory cannot be locked` | 上次运行被杀，`.snakemake/locks` 残留；或另一个 Snakemake 正在同一 `-d` 运行 | 确认没有其他进程后 `wgs_skim run ... -- --unlock`，再运行 | CLI 文档 `--unlock` |
| `IncompleteFilesException` | 中断时正在写的输出被标记为不完整 | `--rerun-incomplete`（profile 已开）；profile 未加载时手动加 | CLI 文档；Ctrl-C 后输出未删见 https://github.com/snakemake/snakemake/issues/1871 |
| `CreateCondaEnvironmentException` | 网络 / 频道 / 求解失败；**Snakemake 8 要求 conda ≥ 24.7.1** | 升级 conda（本机 25.9.1 满足）；重跑 `wgs_skim install`（已建好的保留） | https://github.com/snakemake/snakemake/issues/3163 、/2788 |
| `libmamba Non-conda folder exists at prefix` | 上次建环境失败留下半成品目录 | 删掉 `<conda-prefix>/<hash>_` 半成品目录后重建（先记录在案） | https://github.com/snakemake/snakemake/issues/3249 |
| 改了 `workflow/scripts/*.py` 但规则没重跑 / 无法忽略 code 触发 | `script:` 规则的 code 触发器行为不一致 | 需要时 `--forcerun <rule>`；或 `--rerun-triggers mtime` | https://github.com/snakemake/snakemake/issues/2192 、/3456 |
| rerun-triggers 不传递到下游 | 已知 bug | `--forcerun` 指定上游 | https://github.com/snakemake/snakemake/issues/3095 |
| 命令行 `--default-resources` 整体覆盖 profile 中的设置 | 合并行为 | 命令行给全 | https://github.com/snakemake/snakemake/issues/3429 |
| `--workflow-profile` 指向不存在的目录时不报错 | 已知（open） | 检查日志第一行 "Using workflow specific profile" | https://github.com/snakemake/snakemake/issues/4276 |
| 输出文件已在但 Snakemake 说缺失（NFS） | 文件系统延迟 | `--latency-wait 60` | https://github.com/snakemake/snakemake/issues/39 、/2739 |
| 声明的 mem_mb 大于 `--resources` 总量 | Snakemake 8 把该任务的 mem_mb **静默降为总量**后照常运行（本机 8.18.2 实测：48000 → 20000） | 不会卡住，但工具自己的内存上限（如 SPAdes `-m`）若取自 `{resources.mem_mb}` 也会被降低 | 本项目实测 |

## 6. 本项目中遇到的问题

- **TESTING.md #16**：整机内存紧张时后台运行被停止 → 每个规则声明 `mem_mb`，按核数 + 内存调度；长任务用 `nohup` / `tmux`。
- **TESTING.md #11**：GetOrganelle 中途被杀仍写 `done`，以后不重跑 → 失败时不写标记（与 Snakemake 的 incomplete 机制互补：
  规则内部自建的标记文件 Snakemake 不知道）。
- **TESTING.md #14**：个人配置 `~/.wgs_skim/config.yaml` 的 `db_dir` 被命令行默认值覆盖 → 用户写的值优先（`--config` 优先级问题）。
- 早期 blast 规则失败后 profile 的 `retries: 1` 让它再跑一次（`Trying to restart job 14`，run1.log）：注意 retries 对**所有**规则生效，
  确定性错误（路径错、参数错）也会多跑一次。
- SPAdes 重试：`--restart-from last`；内存 1.5×，依赖 `attempt`（assembly.smk / common.smk）。

## 7. 排错清单

1. 终端第二行 `[wgs_skim] snakemake ...` 是完整命令；加 `-n -p` 看每个任务的 shell。
2. `.snakemake/log/<时间>.snakemake.log`：搜 `Error in rule`、`Exception`、`MissingOutputException`、`Trying to restart job`。
3. 规则日志 `logs/<rule>/<sample>.log`；资源 `benchmarks/<rule>/<sample>.tsv`（max_rss vs 声明的 mem_mb）。
4. 锁：`.snakemake/locks/`；不完整标记：`.snakemake/incomplete/`。
5. 环境：`<conda-prefix>/<hash>_.yaml` 与同名目录；`env_setup_done` 存在才算建好（Nextflow 复用时也看这个）。
6. 版本：`$(dirname $(python -c 'import sys;print(sys.executable)'))/snakemake --version`，不要只看 `which snakemake`。

## 8. 资料来源（查阅日期 2026-10-29）

- https://github.com/snakemake/snakemake/releases （v9.27.1、v9.0.0）
- https://snakemake.readthedocs.io/en/stable/executing/cli.html
- Issues：#39 #1871 #2192 #2739 #2788 #3095 #3163 #3249 #3429 #3456 #4276
- 本机实测（scratchpad 玩具工作流，Snakemake 8.18.2）：资源封顶、workflow profile 查找位置
- 本项目：`src/wgs_skim/cli.py`、`workflow/profiles/default/config.yaml`、`environment.yml`、`pyproject.toml`、
  `docs/TESTING.md` §3.4 §4、README Troubleshooting、`../wgs_skim_test/run1.log`

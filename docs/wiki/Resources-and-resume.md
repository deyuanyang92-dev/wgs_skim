# 资源调度与续跑

所有流程（`wgs_skim run`、`wgs_skim nuclear`、`wgs_skim compare` 批量）都由 Snakemake 管理：按 **核数** 和 **内存** 同时调度，
每个任务一个日志、一条运行时间 / 峰值内存记录，中断后重跑只补缺少的结果。

## 两个总量

```bash
wgs_skim run -r raw/ -o res/ -j 64 --mem-gb 200
```

- `-j`：最多同时使用的核数（默认全部）。
- `--mem-gb`：最多同时占用的内存（默认本机内存的 90%）。传给 Snakemake `--resources mem_mb=…`。

Snakemake 只在 **核数和内存都够** 时启动一个任务。例如 200 GB 时，一个 150 GB 的 SPAdes 运行中，其他 SPAdes 等待，
而 fastp、GetOrganelle、注释这些小任务照常并行 —— 不会因为同时开太多大任务而内存不足、一晚上没有结果。

## 每个任务的用量

| 任务 | 线程 | 内存 |
|---|---|---|
| fastp | 8 | 2 GB |
| meangs | 8 | 12 GB |
| getorganelle_mt（第 0 轮） | 16 | 16 GB |
| mt_iterate（迭代） | 16 | 1.2 × 使用的数据量（≤ `mt_iter_max_data_gb`），至少 8 GB |
| getorganelle_nr / getorganelle_histone | 16 | 1.2 × 最大一轮数据量（`nuclear_max_gb` 10 → 12 GB） |
| spades（genome） | 32 | `spades_gib_per_gb` × clean Gb，至少 `spades_min_mem_gb`；SPAdes `-m` 同值 |
| spades_mode（其他模式） | 16 | `spades_mode_gib_per_gb` × clean Gb，至少 16 GB |
| trinity | 16 | `trinity_gib_per_gb` × clean Gb，至少 16 GB（`--max_memory` 同值） |
| busco / busco_assembly | 16 | 16 GB |
| jellyfish | 16 | 12 GB |
| nuc_markers（SPAdes 上的 MitoM） | 16 | 8 GB |
| mitom_annotate / nr_annotate / histone_annotate | 8 / 4 / 4 | 4 / 2 / 2 GB |
| topology | 8 | 8 GB |
| blast | 8 | 16 GB |

失败的任务自动重试一次（profile `retries: 1`），内存按 1.5 倍（SPAdes、其他模式、Trinity）或 2 倍（GetOrganelle）增加；
SPAdes 从最后的检查点继续。

改线程：`wgs_skim run ... -- --set-threads spades=48 getorganelle_mt=24`，或写进 `workflow/profiles/default/config.yaml`。
改内存估计：`--config spades_gib_per_gb=8`。

## 看进度

```bash
tail -f res/.snakemake/log/*.log            # 总进度："23 of 39 steps (59%) done"
ls res/logs/<任务>/                          # 每个任务的工具输出
cat res/benchmarks/<任务>/<样品>.tsv         # 实际时间（s）、峰值内存（max_rss, MB）
wgs_skim nuclear --status -o nuclear/       # 核基因批量：每个样品每个标记的状态
```

每个样品的结果一完成就写到它的目录里（例如 `mito/<s>/<s>.mt.fasta`），不必等整批结束。
`summary/samples.tsv` 在最后汇总；中途想看，`wgs_skim run -r raw/ -o res/ --steps summary` 汇总已有的结果。

## 续跑与重跑

| 情况 | 做法 |
|---|---|
| 中断（关机、kill、内存不足） | 运行同一条命令；未完成的文件会重做（`rerun-incomplete`） |
| `Directory cannot be locked` | `wgs_skim run ... -- --unlock`，再运行同一条命令 |
| 某个样品失败 | 其他样品继续（`keep-going`）；看 `logs/`，修正后运行同一条命令 |
| 重做某一步 | `wgs_skim run ... -- --forcerun mt_iterate` |
| 只跑到某一步 | `wgs_skim run ... -- --until busco` |
| 改了参数但不想重跑 | `wgs_skim run ... -- --rerun-triggers mtime` |
| 先看会跑什么 | `wgs_skim run ... -n` |
| 加样品 | 把新数据放进同一个目录，运行同一条命令：只跑新样品，汇总表重做 |

长时间运行请放到后台、与终端分离：

```bash
nohup wgs_skim run -r raw/ -o res/ -j 64 --mem-gb 200 > run.log 2>&1 &
# 或者 tmux / screen
```

## 直接用 Snakemake

```bash
snakemake -s workflow/Snakefile --sdm conda --cores 64 --resources mem_mb=200000 \
    --config reads_dir=/path/raw outdir=/path/res
snakemake -s workflow/nuclear_batch.smk --sdm conda --cores 64 --resources mem_mb=200000 \
    --config results=/path/results outdir=/path/nuclear mitom=/path/Mitom
snakemake -s workflow/compare_batch.smk --cores 16 --config results=/path/results outdir=/path/compare
```

集群：加 Snakemake 的执行器插件，如 `--executor slurm`，每个任务的 `threads` / `mem_mb` 会成为作业的资源申请。

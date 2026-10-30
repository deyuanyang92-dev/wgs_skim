# GetOrganelle：线粒体、rDNA、组蛋白

用种子招募 reads、逐步延伸组装，适合高拷贝的串联单元：线粒体、核 rDNA（18S-ITS-28S）、组蛋白基因簇。

## 线粒体（`--target mt`，默认）

```bash
wgs_skim getorganelle -r clean_reads/ -o go_out/                 # GetOrganelle 自带 animal_mt 种子
wgs_skim getorganelle -r clean_reads/ -o go_out/ --meangs        # MEANGS 种子（推荐，= meangs-getorganelle）
wgs_skim getorganelle -r clean_reads/ -o go_out/ --rounds 0      # 不迭代，只跑一轮
```

1. 第 0 轮：`-F animal_mt -R 10 -k 21,45,65,85,105`。
2. 没有成环：以上一轮结果为种子，用更多 reads（≤ `mt_iter_max_data_gb` = 20 Gb）、`-R 20` 再组装，最多 `--rounds`（默认 3）轮；
   成环、或结果不再变化时停止。
- 输出：`mitogenome/<s>/<s>.mt.fasta`（最终线粒体）、`<s>.mt.json`（第几轮、是否成环、停止原因）、`round0..N/`（每轮的 GetOrganelle 结果和组装图）。
- 实测（0.3 Gb）：第 0 轮 2–3 min / 3.6 GB，每轮迭代 2–10 min / 5.8 GB。

## rDNA 与组蛋白（`--target nr`）

```bash
wgs_skim getorganelle -r clean_reads/ -o nr_out/ --target nr
```

`-F anonym`，自带种子 `resources/rdna`、`resources/histone`（可用 `nr_seed`、`histone_seed` 替换）。第 1 轮 5 Gb reads；
单元不完整时第 2 轮 10 Gb、以第 1 轮为种子。之后 MitoM 注释为 GenBank：`rdna/<s>/<s>.rdna.gbf`、`histone/<s>/<s>.histone.gbf`。
详见 [核多拷贝基因](Combo-nuclear.md)。实测（0.3 Gb）：rDNA 3–9 min / 0.9 GB，组蛋白 10–11 min / 2.4 GB。

深入资料：[docs/tools/GetOrganelle.md](../tools/GetOrganelle.md)。

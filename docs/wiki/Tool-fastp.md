# fastp：质控

去接头、去低质量 reads，生成每个样品的质控报告。

```bash
wgs_skim fastp -r raw/ -o qc_out/                 # 目录中全部样品
wgs_skim fastp -s samples.tsv -o qc_out/          # 样品表
```

- 参数：`--detect_adapter_for_pe`（自动识别双端接头），其余为 fastp 默认（质量 Q15、长度 ≥ 15 bp、N ≤ 5）。
- 输出：`qc/<s>/<s>.R1.fq.gz`、`<s>.R2.fq.gz`（clean reads）、`<s>.fastp.html`（报告）、`<s>.fastp.json`（数据量、Q30）。
- 实测：0.3 Gb 样品 30 s、0.4 GB 内存；每个任务 8 线程、2 GB。

其他工具加 `--qc` 时会先运行 fastp；`wgs_skim run` 的组合流程总是先运行 fastp。
已经清洗过的数据，其他工具直接用（不加 `--qc`）：fastp 只统计数据量（所有过滤关闭），reads 原样链接。

深入资料：[docs/tools/fastp.md](../tools/fastp.md)。

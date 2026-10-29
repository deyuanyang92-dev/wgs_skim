# 全流程

```bash
wgs_skim run -r raw/ -o res/                      # 全部（环节动物预设）
wgs_skim run -r raw/ -o res/ --taxon fish -j 64 --mem-gb 300
```

全部 = fastp → MEANGS-GetOrganelle（迭代、SPAdes 兜底）→ MitoM 注释 → reads 验证环 / 线性 → rDNA + 组蛋白（GetOrganelle、MitoM）
→ SPAdes-BUSCO → SPAdes 中的 rDNA / 组蛋白、拷贝数、两种方法比较 → 基因组大小 → BLAST → 汇总。
MitoFinder 默认不运行（加 `--config mitofinder=true` 或 `--steps …,spades-mitofinder`）。

结果：`summary/samples.tsv`（每个样品一行）和 `summary/report.html`，各列见 [输出文件](Output.md)。

双样品测试（各 0.3 Gb）：39 个任务，26 分钟。各步骤实测时间与内存见各工具页面；真实数据（10–50 Gb / 样品）时间随数据量增加，
内存按数据量自动估计（见 [资源调度与续跑](Resources-and-resume.md)）。

中途想看已完成部分的汇总：`wgs_skim run -r raw/ -o res/ --steps summary`。

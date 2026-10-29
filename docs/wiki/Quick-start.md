# 快速开始

```bash
# 1. 看样品名是否识别正确
wgs_skim samples -r /path/to/rawdata

# 2. 全部分析（默认环节动物预设）
wgs_skim run -r /path/to/rawdata -o /path/to/results

# 3. 看结果
#    /path/to/results/summary/samples.tsv     每个样品一行
#    /path/to/results/summary/report.html     浏览器打开
```

长时间运行请放到后台或 tmux，避免终端关闭后中断：

```bash
nohup wgs_skim run -r /path/to/rawdata -o /path/to/results > run.log 2>&1 &
```

中断后**用同一条命令重跑**，已完成的步骤不会重做（见 [[资源调度与续跑|Resources-and-resume]]）。

## 常用变化

```bash
wgs_skim run -r raw/ -o res/ --taxon mollusc             # 软体动物；fish 鱼类；crustacean 甲壳类
wgs_skim run -r raw/ -o res/ -j 64 --mem-gb 300          # 限定 64 核、300 GB 内存
wgs_skim run -r raw/ -o res/ --steps mito,annotate       # 只做线粒体组装 + 注释
wgs_skim run -r raw/ -o res/ --steps nr                  # 只做 rDNA + 组蛋白
wgs_skim run -r rnaseq/ -o res/ --steps rna              # 转录组：rnaSPAdes + Trinity + BUSCO
wgs_skim run -r raw/ -o res/ -n                          # 只看会运行哪些任务（不运行）
```

每个步骤的命令和输出见 [[步骤详解|Steps]]。

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

中断后**用同一条命令重跑**，已完成的步骤不会重做（见 [资源调度与续跑](Resources-and-resume.md)）。

## 一个工具 / 一个组合

```bash
wgs_skim fastp -r raw/ -o qc/                                    # 只做质控
wgs_skim spades -r qc_reads/ -o asm/                             # 只组装（reads 已清洗）
wgs_skim busco -i asm_fastas/ -o busco/                          # 只评估已有组装
wgs_skim mitofinder -i asm_fastas/ -o mf/ --ref related.gb       # 已有组装中找线粒体
wgs_skim run -r raw/ -o res/ --steps meangs-getorganelle         # 线粒体（推荐组合）
wgs_skim run -r raw/ -o res/ --steps spades-busco                # 全基因组
wgs_skim run -r rnaseq/ -o res/ --steps rnaspades-busco          # 转录组
```

## 常用变化

```bash
wgs_skim run -r raw/ -o res/ --taxon mollusc             # 软体动物；fish 鱼类；crustacean 甲壳类
wgs_skim run -r raw/ -o res/ -j 64 --mem-gb 300          # 限定 64 核、300 GB 内存
wgs_skim run -r raw/ -o res/ -n                          # 只看会运行哪些任务（不运行）
```

每个工具见 [单个工具](Tools.md)，推荐组合见 [组合流程](Combinations.md)。

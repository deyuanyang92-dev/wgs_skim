# 与旧脚本（Temp_scripts）的对应

旧脚本的功能都保留，统一成一个由 Snakemake 管理的流程（按核数和内存调度、续跑、每个任务一个日志）。

| 旧脚本 / 仓库 | 本软件 | 变化 |
|---|---|---|
| `use_fastp_cleandata*.py`、`Fastp_Fastqc_mulitqc` | 步骤 `fastp` | 相同参数；fastp.html / json 每个样品一份 |
| `Batch_assembly_of_genomes_using_GetOrganelle*.py`、`batch_getorganelle`、`batch_seeds*.py`（MEANGS → GetOrganelle） | 步骤 `meangs`、`mito`、`getorganelle` | 新增：未成环时以上一轮为种子迭代；SPAdes scaffold 闭合时兜底 |
| `batch_obation_nr_using_getorg.py`（18S-28S）、Hesionidae H3 种子 | 步骤 `nr`（rDNA + 组蛋白）、`wgs_skim nuclear` | 新增：5 / 10 Gb 两轮、MitoM GenBank `.gbf`、与 SPAdes-MitoM 比较、组蛋白类别 |
| `spadesv3*.py`（wgs / rna / meta / bio / sewage / plasmid / metaviral / metaplasmid、--careful、--only-assembler、-k）、`fastp_spades_busco.py`、`busco2summary` | 步骤 `spades` + `busco`，`--spades-modes` | 全部 SPAdes 模式，每个模式一个组装 + BUSCO；内存按数据量估计；断点续跑 |
| `trinity_batch3*.py`、rnaspades snakefile、snakemake_transcriptome_workflow | 步骤 `trinity`、`rna` | rnaSPAdes + Trinity + BUSCO（transcriptome）；TransDecoder / 功能注释未包含 |
| `jeyllfish_genome*.py`、`genome_size.py` | 步骤 `genome_size` | jellyfish + GenomeScope2 |
| `Mitoz-annotate` / `mitoz_get*.py`、`batch_mitos2.py` | 步骤 `annotate` | 用 MitoM 注释（蛋白基因、rRNA、tRNA、NCBI 检查） |
| `mt_circularize.py`、`check_circularity.pl`、MITObim `circules.py` | 步骤 `topology`（reads）、`wgs_skim circularize`（组装图，无需 reads） | 新增：跨连接点 reads 验证、末端诊断 |
| `blast2metadata` | 步骤 `blast`、`wgs_skim build-blastdb` | 本地 nt 或 NCBI 在线；按基因查询（COX1、18S、ITS、28S）；小型标记库 |
| `mitofinder-batch`、`spades_mitofinder.py`、`mitofinder_annotate2summary.py` | `wgs_skim compare` 读取 MitoFinder 结果并与 GetOrganelle 比较 | MitoFinder 本身不再运行（MitoM 注释） |
| `batch_report.py` / `report_pdf.py` / `report_web.py` | `summary/samples.tsv` + `summary/report.html` | 无 PDF |
| `using_bwa_draw_coverage*.py`、`organelle_coverage.py`、`plot_m2_depth.py`（深度图） | 尚未包含 | MitoM `mito_plot.py` 可画带深度环的基因组图 |
| `remove_getorganelle_templefile*.py`、`remove_spades_temp.sh` | 自动 | 每步完成后删除 GetOrganelle 中间文件、SPAdes K* 目录 |
| `Temp_scripts/wgs_skim_pipeline`（旧版一体化流程） | `wgs_skim run`；其结果目录可直接作为 `wgs_skim nuclear --results`、`wgs_skim compare --results` 的输入 | 不改动旧结果 |

## 旧命令 → 新命令

```bash
# spadesv3.5.5.py -m rna ...（批量 rnaSPAdes）
wgs_skim run -r rnaseq/ -o res/ --steps spades,busco --spades-modes rna
# spadesv3.5.5.py -m meta --careful ...
wgs_skim run -r raw/ -o res/ --steps spades --spades-modes meta
# trinity_batch3.1.py
wgs_skim run -r rnaseq/ -o res/ --steps trinity,busco
# fastp_spades_busco.py
wgs_skim run -r raw/ -o res/ --steps spades,busco
# batch_seeds + GetOrganelle
wgs_skim run -r raw/ -o res/ --steps mito
# batch_obation_nr_using_getorg.py
wgs_skim run -r raw/ -o res/ --steps nr
# blast2metadata
wgs_skim run -r raw/ -o res/ --steps blast --blast-db /db/blast_marker/metazoa_marker
```

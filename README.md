# wgs_skim

A Snakemake workflow for **genome skimming of animals**: from raw paired-end reads to mitogenomes, nuclear rDNA
(18S-ITS-28S), a whole-genome survey and a species check, with one command and no manual software setup.

Developed by the YZZ Team (De-Yuan Yang) for annelids, molluscs, fish and crustaceans.

```
raw reads ─ fastp ─┬─ MEANGS seed ─ GetOrganelle mitogenome ─ not circular? reseed with previous result, all reads (≤20 Gb) ─┐
                   ├─ GetOrganelle 18S-ITS-28S rDNA                                                                         ├─ MitoZ annotation ─ read-based circular/linear check
                   ├─ SPAdes ─ BUSCO                                                                                         │
                   ├─ jellyfish ─ GenomeScope2                                                                               ├─ BLAST + GenBank taxonomy
                   └────────────────────────────────────────────────────────────────────────────────────────────────────────┴─ summary table + HTML report
```

## 中文快速上手

```bash
# 1. 安装（一次）
git clone https://github.com/deyuanyang92-dev/wgs_skim.git && cd wgs_skim
conda env create -f environment.yml && conda activate wgs_skim
pip install .
wgs_skim install                 # 自动安装 GetOrganelle、MEANGS、SPAdes、BUSCO、MitoZ 等全部软件（各自独立的 conda 环境）

# 2. 运行：只给原始数据目录，样品名自动识别
wgs_skim run -r /path/to/rawdata                     # 默认环节动物
wgs_skim run -r /path/to/rawdata --taxon mollusc     # 螺类 / 贝类；fish 鱼类；crustacean 甲壳类
```

结果在 `<rawdata 上级目录>/wgs_skim_results/`，汇总表 `summary/samples.tsv`，报告 `summary/report.html`。
中断后用同一条命令重跑即可从断点继续。

## Installation

Requirements: Linux, [conda](https://github.com/conda-forge/miniforge) (Miniforge recommended), internet access for the
first run. All analysis software is installed automatically from conda-forge / bioconda into separate environments.

```bash
git clone https://github.com/deyuanyang92-dev/wgs_skim.git
cd wgs_skim
conda env create -f environment.yml      # python + snakemake
conda activate wgs_skim
pip install .                            # the `wgs_skim` command
wgs_skim install                         # create all tool environments once (20-60 min), shared by every run
```

Environments are kept in `~/.wgs_skim/conda` (change with `--conda-prefix` or `WGS_SKIM_CONDA_PREFIX`), downloaded
databases (GetOrganelle animal_mt, BUSCO lineage) in `~/.wgs_skim/db`.

## Usage

### With the `wgs_skim` command

```bash
wgs_skim run -r rawdata/                         # everything automatic
wgs_skim run -r rawdata/ -o results/ -j 64       # output directory, CPU cores
wgs_skim run -r rawdata/ --taxon fish            # taxon preset (annelid | mollusc | fish | crustacean)
wgs_skim run -r rawdata/ -n                      # dry run: list the jobs without running them
wgs_skim run -s samples.tsv                      # sample sheet instead of a directory
wgs_skim run -r rawdata/ --blast-db /db/nt/nt --email you@example.com
wgs_skim run -r rawdata/ --config mt_iter_max_data_gb=0 assembly=false   # any option, Snakemake style
wgs_skim run -r rawdata/ -- --set-threads spades=48                      # after --: passed to snakemake
wgs_skim samples -r rawdata/                     # show detected samples
wgs_skim config > my.yaml                        # every option with comments; use with --configfile my.yaml
```

| option | meaning | default |
|---|---|---|
| `-r, --reads_dir` | directory of paired FASTQ files (flat, or one sub-directory per sample) | – |
| `-s, --samples` | sample sheet (TSV/CSV: `sample`, `fastq_1`, `fastq_2`) instead of `-r` | – |
| `-o, --outdir` | results directory | `<reads_dir>/../wgs_skim_results` |
| `-t, --taxon` | `annelid`, `mollusc`, `fish`, `crustacean` | `annelid` |
| `-j, --cores` | CPU cores | all |
| `-n, --dry-run` | show jobs only | off |
| `--blast-db` | local nt database prefix | NCBI online BLAST |
| `--email` | e-mail for NCBI Entrez | – |
| `--busco-lineage` | BUSCO lineage name or directory | from taxon |
| `--configfile`, `--config KEY=VALUE` | any workflow option (as in Snakemake) | `config/config.yaml` |

Personal defaults that apply to every run (local nt, e-mail, a downloaded BUSCO lineage ...) can be written once in
`~/.wgs_skim/config.yaml`, e.g.

```yaml
blast_db: /data/db/nt/nt
email: you@example.com
busco_lineage: /data/db/busco/metazoa_odb12
```

### With Snakemake directly

The repository is a standard Snakemake workflow (`workflow/Snakefile`, `workflow/rules/*.smk`, one conda environment
per tool in `workflow/envs/`, defaults in `config/config.yaml`, profile in `workflow/profiles/default`):

```bash
snakemake --cores 32 --sdm conda --config reads_dir=/path/to/rawdata
snakemake --cores 32 --sdm conda --configfile my.yaml -n          # dry run
snakemake --cores 32 --sdm conda --config reads_dir=raw taxon=fish --set-threads spades=48
```

Cluster: add a Snakemake executor profile, e.g. `--executor slurm --jobs 50` (snakemake-executor-plugin-slurm).

### Input

Sample names come from the file names. These are recognised (sample name = part before the read token):
`S_R1_001.fastq.gz`, `S_R1.fq.gz`, `S.R1.raw.fastq.gz` (Novogene / Majorbio), `S_1.fq.gz`, `S_1.clean.fq.gz` (BGI),
`S_S1_L001_R1_001.fastq.gz` (Illumina), `S_R1_trimmed.fq.gz` (Macrogen), `S.1.fastq` ..., in the directory itself or
one sub-directory per sample. Check with `wgs_skim samples -r DIR`. Anything else: use a sample sheet.

### Taxon presets

| `taxon` | MEANGS class | mito genetic code | MitoZ clade | BUSCO lineage |
|---|---|---|---|---|
| `annelid` (default) | A-worms | 5 | Annelida-segmented-worms | metazoa_odb12 |
| `mollusc` | Mollusca | 5 | Mollusca | mollusca_odb12 |
| `fish` | Chordata | 2 | Chordata | actinopterygii_odb12 |
| `crustacean` | Arthropoda | 5 | Arthropoda | arthropoda_odb12 |

Each value can be overridden (`genetic_code`, `meangs_class`, `mitoz_clade`, `busco_lineage`), which also covers other
animal groups.

### Mitogenome strategy

1. MEANGS finds a seed (if it fails, GetOrganelle uses its own animal_mt seeds).
2. GetOrganelle assembles the mitogenome (`-R 10`, GetOrganelle's default read subsample).
3. If the result is not circular, GetOrganelle is rerun up to `mt_iter_rounds` (3) times, each time seeded with the
   previous result, with `-R 20` and all clean reads — capped at `mt_iter_max_data_gb` (20 Gb; `0` = no cap).
   The first circular round is kept, otherwise the last round.
4. `topology` checks the result with the reads: the circle is cut inside *nad5* and reads spanning the original
   ends must exist for a real circle (Ahuja et al. 2024, GBE); `linear_unjoined` = no read joins the ends (true
   linear, or a gap longer than the insert size — confirm with long-range PCR).

## Output

```
wgs_skim_results/
├── summary/samples.tsv        one row per sample: clean data, mitogenome length / circular / round, genes,
│   summary/report.html        rDNA length, topology verdicts, SPAdes N50, BUSCO, genome size, heterozygosity, top BLAST hits
├── qc/<sample>/               clean reads, fastp.html / fastp.json
├── mito/<sample>/<sample>.mt.fasta   final mitogenome (+ .mt.json: circular, chosen round); round0..N GetOrganelle runs
├── rdna/<sample>/<sample>.nr.fasta   18S-ITS-28S
├── annotation/<sample>/<sample>.mt.gb   MitoZ annotation (GenBank)
├── topology/<sample>/         read-based circular / linear check
├── assembly/<sample>/         SPAdes scaffolds, BUSCO short summary
├── genome_size/<sample>/      k-mer histogram, GenomeScope2 plots and summary
├── blast/blast_topn.tsv       top hits with GenBank taxonomy
└── logs/  benchmarks/         one log and one runtime/memory record per job
```

## Resuming and re-running

Snakemake keeps track of every output file: run the same command again after an interruption or an error and only
the missing or out-of-date jobs run. Useful Snakemake options (`wgs_skim run ... -- <option>`):
`--forcerun mt_iterate` (redo a step), `--until busco` (stop after a step), `--unlock` (after a crash),
`--rerun-triggers mtime` (ignore changed parameters).

## Software

fastp, MEANGS, GetOrganelle, SPAdes, BUSCO, Jellyfish, GenomeScope2, MitoZ, BWA, SAMtools, BLAST+ and
[blast2metadata](https://github.com/deyuanyang92-dev/Temp_scripts/tree/main/blast2metadata). Versions are pinned in
`workflow/envs/*.yaml`. Please cite the tools you use, and

> Yang D, Zeng S, Wang Z, Zhang Y, Yang D, Glasby CJ, Hwang J-S, Cai L (2024) Molecular systematics of *Perinereis* …
> *Zoosystematics and Evolution* 100(4): 1297–1314. https://doi.org/10.3897/zse.100.127201

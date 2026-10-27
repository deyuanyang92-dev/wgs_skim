# wgs_skim

A Snakemake workflow for **genome skimming of animals**: from raw paired-end reads to mitogenomes, nuclear rDNA
(18S-ITS-28S), a whole-genome survey and a species check, with one command and no manual software setup.

Developed by the YZZ Team (De-Yuan Yang) for annelids, molluscs, fish and crustaceans.

```
raw reads ─ fastp ─┬─ MEANGS seed ─ GetOrganelle mitogenome ─ not circular? reseed with previous result, all reads (≤20 Gb) ─┐
                   ├─ GetOrganelle 18S-ITS-28S rDNA                                                                         ├─ MitoM annotation ─ read-based circular/linear check
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
wgs_skim install                 # 自动安装 GetOrganelle、MEANGS、SPAdes、BUSCO、MitoM 等全部软件（各自独立的 conda 环境）

# 2. 运行：只给原始数据目录，样品名自动识别
wgs_skim run -r /path/to/rawdata                     # 默认环节动物
wgs_skim run -r /path/to/rawdata --taxon mollusc     # 螺类 / 贝类；fish 鱼类；crustacean 甲壳类
wgs_skim run -r /path/to/rawdata --steps fastp,spades,busco   # 只跑某几步（自动带上依赖的上游步骤）
wgs_skim run -r /path/to/rawdata --steps getorganelle         # 只用 GetOrganelle 组装线粒体（不用 MEANGS 种子）
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
wgs_skim run -r rawdata/ --steps fastp,spades,busco   # only some steps (plus the steps they need)
wgs_skim run -s samples.tsv                      # sample sheet instead of a directory
wgs_skim run -r rawdata/ --blast-db /db/nt/nt --email you@example.com
wgs_skim run -r rawdata/ --blast-mt cox1,16s,full --blast-nr 18S,ITS     # several BLAST queries
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
| `--steps` | run only these steps, comma separated (see below) | all |
| `--blast-db` | local nt database prefix | NCBI online BLAST |
| `--blast-mt` | mitogenome BLAST queries: `full` and/or genes (`cox1` = COI, `16s`, `12s`, `cytb`, `nad1` ...) | `cox1` |
| `--blast-nr` | rDNA BLAST queries: `full` and/or `18S`, `ITS` (ITS1-5.8S-ITS2), `ITS1`, `5.8S`, `ITS2`, `28S` | `18S,ITS,28S` |
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

### Running single steps

`--steps` (or `steps:` in a config file, `--config steps=...` with Snakemake) runs only the named steps; the steps
they depend on run too, and finished ones are reused.

| step | runs | needs |
|---|---|---|
| `fastp` | read QC and trimming | – |
| `meangs` | MEANGS mitochondrial seed | fastp |
| `mito` | MEANGS → GetOrganelle → reseeded rounds if not circular | fastp |
| `getorganelle` | the same without MEANGS (GetOrganelle animal_mt seeds) | fastp |
| `nr` | GetOrganelle 18S-ITS-28S | fastp |
| `spades` | SPAdes whole-genome assembly | fastp |
| `busco` | BUSCO on the SPAdes assembly | spades |
| `genome_size` | jellyfish + GenomeScope2 | fastp |
| `annotate` | MitoM annotation | mito |
| `topology` | read-based circular / linear check | mito, annotate |
| `blast` | BLAST + GenBank taxonomy of mitogenome and rDNA | mito, nr |
| `summary` | summary table and report of what was run | – |

Aliases: `qc`, `meangs-getorganelle`, `rdna`, `assembly`, `jellyfish`, `genomescope`, `mitom`, `report`.

### With Snakemake directly

The repository is a standard Snakemake workflow (`workflow/Snakefile`, `workflow/rules/*.smk`, one conda environment
per tool in `workflow/envs/`, defaults in `config/config.yaml`, profile in `workflow/profiles/default`):

```bash
snakemake --cores 32 --sdm conda --config reads_dir=/path/to/rawdata
snakemake --cores 32 --sdm conda --configfile my.yaml -n          # dry run
snakemake --cores 32 --sdm conda --config reads_dir=raw taxon=fish --set-threads spades=48
snakemake --cores 32 --sdm conda --config reads_dir=raw steps=fastp,spades,busco
```

Cluster: add a Snakemake executor profile, e.g. `--executor slurm --jobs 50` (snakemake-executor-plugin-slurm).

### Input

Sample names come from the file names. These are recognised (sample name = part before the read token):
`S_R1_001.fastq.gz`, `S_R1.fq.gz`, `S.R1.raw.fastq.gz` (Novogene / Majorbio), `S_1.fq.gz`, `S_1.clean.fq.gz` (BGI),
`S_S1_L001_R1_001.fastq.gz` (Illumina), `S_R1_trimmed.fq.gz` (Macrogen), `S.1.fastq` ..., in the directory itself or
one sub-directory per sample. Check with `wgs_skim samples -r DIR`. Anything else: use a sample sheet.

### Taxon presets

| `taxon` | MEANGS class | mito genetic code | MitoM clade | BUSCO lineage |
|---|---|---|---|---|
| `annelid` (default) | A-worms | 5 | Annelida-segmented-worms | metazoa_odb12 |
| `mollusc` | Mollusca | 5 | Mollusca | mollusca_odb12 |
| `fish` | Chordata | 2 | Chordata | actinopterygii_odb12 |
| `crustacean` | Arthropoda | 5 | Arthropoda | arthropoda_odb12 |

Each value can be overridden (`genetic_code`, `meangs_class`, `clade`, `busco_lineage`), which also covers other
animal groups.

### Mitogenome strategy

1. MEANGS finds a seed in the sample's own reads (if it fails, GetOrganelle uses its animal_mt seeds).
2. GetOrganelle assembles the mitogenome (`-R 10`, GetOrganelle's default read subsample).
3. Not circular: GetOrganelle is rerun, seeded with the previous result, `-R 20`, all clean reads up to
   `mt_iter_max_data_gb` (20 Gb). It stops at the first circular round **or as soon as a round returns the same
   sequence as the round before** (more rounds of the same kind cannot help; `mt.json` records why it stopped).
4. Still not circular: the SPAdes whole-genome assembly is searched for the mitogenome (`mt_from_spades`). A scaffold
   that contains the GetOrganelle result and closes into a circle replaces it — whole-genome assembly has no read
   recruitment step and can get through control regions where GetOrganelle stops. Longer *linear* scaffolds are not
   used (their extra ends were collapsed repeats or unsupported in tests), only reported.
5. MitoM annotates the result (cox1 start, NCBI checks).
6. `topology` checks the result with the reads: the circle is cut inside *nad5* and reads spanning the original
   ends must exist for a real circle (Ahuja et al. 2024, GBE). For linear results `end_diagnosis` says why each end
   stops — e.g. "reads continue into a simple repeat ((AT)n)", "coverage drops", "homopolymer" — which tells you
   whether long-range PCR / long reads are needed.
7. BLAST (species check) uses the queries you choose, several at once: whole sequences (`full`) and/or genes cut
   from the annotations — mitogenome genes from MitoM (default `cox1`, the COI barcode; `16s`, `12s`, `cytb` ...),
   rDNA regions from `mitom nuc` (default `18S,ITS,28S`). Hits come with GenBank organism and taxonomy.

## Output

```
wgs_skim_results/
├── summary/samples.tsv        one row per sample: clean data, mitogenome length / circular / round, genes,
│   summary/report.html        rDNA length, topology verdicts, SPAdes N50, BUSCO, genome size, heterozygosity, top BLAST hits
├── qc/<sample>/               clean reads, fastp.html / fastp.json
├── mito/<sample>/<sample>.mt.fasta   final mitogenome (+ .mt.json: source, circular, why iteration stopped); round0..N
├── rdna/<sample>/<sample>.nr.fasta   18S-ITS-28S (+ .nr.gb / .nr_features.tsv: 18S ITS1 5.8S ITS2 28S, mitom nuc)
├── annotation/<sample>/<sample>.mt.gb   MitoM annotation (GenBank; + .mitom.tsv gene counts)
├── topology/<sample>/         read-based circular / linear check, end diagnosis
├── assembly/<sample>/         SPAdes scaffolds, BUSCO short summary
├── genome_size/<sample>/      k-mer histogram, GenomeScope2 plots and summary
├── blast/blast_topn.tsv       COX1 / rDNA top hits with GenBank organism and taxonomy
└── logs/  benchmarks/         one log and one runtime/memory record per job
```

## Resuming and re-running

Snakemake keeps track of every output file: run the same command again after an interruption or an error and only
the missing or out-of-date jobs run. Useful Snakemake options (`wgs_skim run ... -- <option>`):
`--forcerun mt_iterate` (redo a step), `--until busco` (stop after a step), `--unlock` (after a crash),
`--rerun-triggers mtime` (ignore changed parameters).

## Software

fastp, MEANGS, GetOrganelle, SPAdes, BUSCO, Jellyfish, GenomeScope2, [MitoM](https://github.com/deyuanyang92-dev/Mitom), BWA, SAMtools, BLAST+ and
Biopython (GenBank taxonomy of the BLAST hits). Versions are pinned in
`workflow/envs/*.yaml`. Please cite the tools you use, and

> Yang D, Zeng S, Wang Z, Zhang Y, Yang D, Glasby CJ, Hwang J-S, Cai L (2024) Molecular systematics of *Perinereis* …
> *Zoosystematics and Evolution* 100(4): 1297–1314. https://doi.org/10.3897/zse.100.127201

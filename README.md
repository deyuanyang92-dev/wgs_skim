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

**详细中文使用手册（每个功能的命令、实测时间与内存、输出、报错排查）：[docs/wiki](docs/wiki/README.md)**


```bash
# 1. 安装（一次）
git clone https://github.com/deyuanyang92-dev/wgs_skim.git && cd wgs_skim
conda env create -f environment.yml && conda activate wgs_skim
pip install .
wgs_skim install                 # 自动安装 GetOrganelle、MEANGS、SPAdes、BUSCO、MitoM 等全部软件（各自独立的 conda 环境）及数据库

# 2. 运行：只给原始数据目录，样品名自动识别
wgs_skim run -r /path/to/rawdata                     # 默认环节动物
wgs_skim run -r /path/to/rawdata --taxon mollusc     # 螺类 / 贝类；fish 鱼类；crustacean 甲壳类
wgs_skim run -r /path/to/rawdata --steps fastp,spades,busco   # 只跑某几步（自动带上依赖的上游步骤）
wgs_skim run -r /path/to/rawdata --steps getorganelle         # 只用 GetOrganelle 组装线粒体（不用 MEANGS 种子）
wgs_skim run -r /path/to/rawdata --steps nr                   # 核 rDNA（18S-ITS-28S）+ 组蛋白基因簇，MitoM 注释成 .gbf
wgs_skim run -r /path/to/rawdata --mem-gb 200 -j 64           # 限定内存 / 核数，Snakemake 按两者调度并行

# 3. 对已有结果（含旧流程 Temp_scripts 的结果目录）
wgs_skim nuclear --results results/ -o nuclear/               # 所有样本的 rDNA + 组蛋白：组装、注释、与 SPAdes 比较
wgs_skim nuclear --status -o nuclear/                         # 查看批量进度（运行中也可以）
wgs_skim compare --results results/ -o compare/               # GetOrganelle 与 SPAdes-MitoFinder 线粒体是否一致
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
| `-s, --samples` | sample sheet (TSV/CSV: `sample`, `fastq_1`, `fastq_2`; optional `scaffolds` = an existing SPAdes assembly, linked instead of assembled) instead of `-r` | – |
| `-o, --outdir` | results directory | `<reads_dir>/../wgs_skim_results` |
| `-t, --taxon` | `annelid`, `mollusc`, `fish`, `crustacean` | `annelid` |
| `-j, --cores` | CPU cores | all |
| `--mem-gb` | memory all jobs may use together; each rule declares its need (SPAdes from the data size), so Snakemake runs as many jobs side by side as cores *and* memory allow | 90 % of the machine |
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

#### A fast BLAST database for species checks

BLAST reads the whole database for every search and is disk-bound: full nt (3.4 Tbases, ~1.2 TB) took 38-52 min
for one COX1 + rDNA search on MJ20260820041, with 8 or 32 threads and with or without `-taxids` Metazoa. Barcode
and marker queries are matched by short records (mitogenomes, rDNA, gene sequences), so build a Metazoa database of
records up to 50 kb once from the local nt:

```bash
wgs_skim build-blastdb --nt /data/db/nt/nt --out /data/db/metazoa_marker/metazoa_marker   # --taxid, --max-length
```

| database | size | search (same queries) | COX1 top 5 | rDNA top hit |
|---|---|---|---|---|
| nt | ~1.2 TB | 2273-2545 s | NC_063122 79.9% ... | KY753844 *Pisione* sp. 93.6% |
| Metazoa ≤ 50 kb (71.3 M records) | 68 GB | 121 s first, 22 s cached | identical | identical |

Left out: chromosome-scale genome assemblies (e.g. Darwin Tree of Life OZ chromosomes that carry an rDNA copy). They
ranked 2nd-5th for rDNA here, behind the same top hit; point `blast_db` at full nt when a species is known only from
such an assembly. Building needs a temporary FASTA of ~210 GB (deleted at the end) and several hours; makeblastdb
alone took 1.7 h here.

### One tool on its own, for a whole batch

Every tool has its own command. It runs on every sample of a directory (scheduled by Snakemake, one log per job,
resumable); reads tools use the reads as they are (already clean; `--qc` cleans them with fastp first), FASTA
tools take files or directories with `-i` (sample = file name; `--split-by-id` splits one merged FASTA by header).

```bash
wgs_skim fastp        -r raw/        -o qc/                        # fastp cleaning + report
wgs_skim meangs       -r clean/      -o meangs/                    # MEANGS mitochondrial seed
wgs_skim getorganelle -r clean/      -o go/ [--meangs] [--rounds 3] [--target mt|nr]
wgs_skim spades       -r clean/      -o asm/ [--mode genome,rna,meta,plasmid,...]
wgs_skim trinity      -r clean/      -o trinity/
wgs_skim genome-size  -r clean/      -o gs/                        # jellyfish + GenomeScope2
wgs_skim busco        -i assemblies/ -o busco/ [-m genome|transcriptome|proteins] [-l lineage]
wgs_skim mitofinder   -i assemblies/ -o mf/ --ref related.gb      # or -r reads/ [--assembler megahit|metaspades|spades|rnaspades|trinity]
wgs_skim mitom        -i mitogenomes/ -o mitom/                    # MitoM annotation
wgs_skim blast        -i seqs/       -o blast/ [--genes full|cox1,...] [--db local_db]
```

Common options: `-j` cores and `--mem-gb` memory for all jobs together, `-t` taxon preset, `-n` dry run,
`--config KEY=VALUE`, anything after `--` to Snakemake. Chinese manual per tool: [docs/wiki/Tools.md](docs/wiki/Tools.md).

### Combinations (recommended)

| goal | combination | `--steps` |
|---|---|---|
| mitogenome | MEANGS-GetOrganelle, reseeded rounds until circular | `meangs-getorganelle` |
| mitogenome (second method) | SPAdes-MitoFinder (reference: `--mitofinder-ref`, or the sample's own MitoM annotation with `meangs-getorganelle`) | `spades-mitofinder` |
| whole genome | SPAdes-BUSCO (fastp first) | `spades-busco` / `fastp-spades-busco` |
| transcriptome | rnaSPAdes-BUSCO (recommended), Trinity-BUSCO; both: `rna` | `rnaspades-busco` / `trinity-busco` |
| mt genes from RNA-seq | rnaSPAdes-MitoFinder, Trinity-MitoFinder; all samples gathered per gene (`mitofinder/<source>/genes/`) | `rnaspades-mitofinder` / `trinity-mitofinder` |

```bash
wgs_skim run -r raw/ -o res/ --steps meangs-getorganelle,spades-mitofinder,spades-busco
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
| `blast` | BLAST + GenBank taxonomy of mitogenome and rDNA | mito, nr (+ annotate for gene queries) |
| `nuclear` | rDNA + histones in the SPAdes assembly, copy numbers, comparison with the GetOrganelle units | spades, busco, nr |
| `rnaspades` | rnaSPAdes (SPAdes `--rna`) | fastp |
| `mitofinder` | MitoFinder on the SPAdes assembly (or its own MEGAHIT assembly) | spades, (annotate as reference) |
| `trinity` | Trinity transcriptome assembly | fastp |
| `rna` | transcriptome: rnaSPAdes + Trinity + BUSCO `-m transcriptome` on both | fastp |
| `summary` | summary table and report; alone: of everything the results directory already has | – |

`nr` also assembles and annotates the histone cluster. Aliases: `qc`, `mitogenome` (mito),
`rdna` / `histone` / `histones` (nr), `assembly` (spades), `jellyfish` / `genomescope` (genome_size), `mitom` / `mitoz`
(annotate), `nuc` (nuclear), `report` (summary).

### Recipes: one step for a whole batch

Every command below works on a directory of raw reads of any number of samples (sample names from the file names),
reuses what an earlier command in the same results directory already made, and writes one log per sample under
`logs/<step>/`. Replace `-r rawdata/` by `-s samples.tsv` for a sample sheet.

```bash
# QC only (fastp): clean reads + fastp.html per sample
wgs_skim run -r rawdata/ --steps fastp                   # -> qc/<sample>/<sample>.R1.fq.gz, .fastp.html/.json

# SPAdes in any mode (as spadesv3.5.5.py): genome (default) isolate sc meta rna rnaviral plasmid metaplasmid
# metaviral bio sewage corona; several at once, each with BUSCO (-m transcriptome for rna / rnaviral)
wgs_skim run -r rawdata/ --steps spades,busco --spades-modes genome,meta,plasmid   # -> assembly/<s>/<s>.meta.fasta ...
wgs_skim run -r rawdata/ --steps spades --spades-modes isolate --config spades_careful=true spades_kmers=21,33,55

# transcriptome (RNA-seq reads): rnaSPAdes and Trinity, BUSCO on each
wgs_skim run -r rnaseq/ --steps rna                     # -> assembly/<s>/<s>.rna.fasta, <s>.trinity.fasta, *.busco.txt
wgs_skim run -r rnaseq/ --steps trinity,busco --spades-modes rna --config trinity_extra="--min_contig_length 300"

# SPAdes whole-genome assembly, then BUSCO on it
wgs_skim run -r rawdata/ --steps spades                  # -> assembly/<sample>/<sample>.scaffolds.fasta
wgs_skim run -r rawdata/ --steps busco                   # -> assembly/<sample>/<sample>.busco.txt (runs spades first)
wgs_skim run -r rawdata/ --steps spades -- --set-threads spades=48 --set-resources spades:mem_mb=250000

# mitogenome: MEANGS seed -> GetOrganelle -> reseeded rounds until circular (-> SPAdes scaffold if that closes)
wgs_skim run -r rawdata/ --steps mito                    # -> mito/<sample>/<sample>.mt.fasta + .mt.json
wgs_skim run -r rawdata/ --steps getorganelle            # the same with GetOrganelle's own animal_mt seeds, no MEANGS
wgs_skim run -r rawdata/ --steps annotate                # + MitoM annotation -> annotation/<sample>/<sample>.mt.gb
wgs_skim run -r rawdata/ --steps topology                # + read-based circular / linear check -> topology/

# nuclear multi-copy markers: 18S-ITS-28S and the histone cluster (GetOrganelle, MitoM .gbf)
wgs_skim run -r rawdata/ --steps nr                      # -> rdna/<sample>/<sample>.rdna.gbf, histone/<sample>/<sample>.histone.gbf
wgs_skim run -r rawdata/ --steps nuclear                 # + the same genes in the SPAdes assembly, copy numbers,
                                                         #   GetOrganelle vs SPAdes per gene -> nuclear/<sample>/

# species check: COX1 and rDNA genes against nt, top hits with GenBank taxonomy
wgs_skim run -r rawdata/ --steps blast --blast-db /db/blast_marker/metazoa_marker   # local (fast, see build-blastdb)
wgs_skim run -r rawdata/ --steps blast --email you@example.com                      # NCBI online BLAST (slow)
wgs_skim run -r rawdata/ --steps blast --blast-mt cox1,16s,cytb --blast-nr 18S,ITS,28S

# genome size and heterozygosity (jellyfish + GenomeScope2)
wgs_skim run -r rawdata/ --steps genome_size             # -> genome_size/<sample>/<sample>.genomescope.txt + plots

# the report of everything made so far (one row per sample + HTML; with other steps: of those steps)
wgs_skim run -r rawdata/ --steps summary                 # -> summary/samples.tsv, summary/report.html
```

Several steps at once: `--steps fastp,spades,busco`. Everything: no `--steps`. Only some samples: a sample sheet
with those rows. Rerun one step after a change: `-- --forcerun busco`.

Outside the workflow, on results made by other pipelines:

```bash
wgs_skim nuclear --results old_results/ -o nuclear/            # rDNA + histones of every sample (see above)
wgs_skim compare --results old_results/ -o compare/            # GetOrganelle vs SPAdes-MitoFinder mitogenomes
wgs_skim circularize -i <sample>_all_mt.fasta --spades spades/ # close MitoFinder / SPAdes contigs without reads
wgs_skim build-blastdb --nt /db/blast/nt --out /db/blast_marker/metazoa_marker   # once: fast marker BLAST database
```

### With Snakemake directly

The repository is a standard Snakemake workflow (`workflow/Snakefile`, `workflow/rules/*.smk`, one conda environment
per tool in `workflow/envs/`, defaults in `config/config.yaml`, profile in `workflow/profiles/default`):

```bash
snakemake --cores 32 --sdm conda --config reads_dir=/path/to/rawdata
snakemake --cores 32 --sdm conda --configfile my.yaml -n          # dry run
snakemake --cores 32 --sdm conda --config reads_dir=raw taxon=fish --set-threads spades=48
snakemake --cores 32 --sdm conda --config reads_dir=raw steps=fastp,spades,busco
```

Scheduling: Snakemake runs jobs in parallel within `--cores` and `--mem-gb`. Every rule declares `threads` and
`resources: mem_mb` (measured peaks are in `workflow/rules/common.smk`; GetOrganelle and SPAdes scale with the data),
so on a 128-core / 880 GB machine many samples run at once, and on a 16-core / 64 GB laptop they queue. Change a
rule's share with `-- --set-threads spades=48` or `-- --set-resources spades:mem_mb=200000`.

Cluster: add a Snakemake executor profile, e.g. `--executor slurm --jobs 50` (snakemake-executor-plugin-slurm).

### With Nextflow

`nextflow/main.nf` (DSL2, Nextflow >= 24.04) runs the same steps with the same scripts (`workflow/scripts`), the same
conda environments (`workflow/envs`) and the same output layout. The easiest way is the same command with
`--engine nextflow`:

```bash
wgs_skim run --engine nextflow -r /path/to/rawdata -o results -j 32 --mem-gb 120
```

It passes every option of `wgs_skim run` as a Nextflow parameter (`-j` -> `--max_cpus`, `--mem-gb` -> `--max_memory`,
`--config KEY=VALUE` -> `--KEY VALUE`; `~/.wgs_skim/config.yaml` and `--configfile` -> `-params-file`), always with
`-resume`, the task directories in `results/work` and Nextflow's own files (`.nextflow/`, `.nextflow.log`) in
`results/`. Directly:

```bash
nextflow run nextflow/main.nf -profile conda --reads_dir /path/to/rawdata --outdir results \
    --max_cpus 32 --max_memory '120 GB' -resume
nextflow run nextflow/main.nf -profile conda -params-file my.yaml -resume     # my.yaml: `wgs_skim config > my.yaml`
nextflow run nextflow/main.nf -profile conda --samples sheet.tsv --steps fastp,spades,busco --taxon fish -resume
```

- **Parameters** are the keys of `config/config.yaml` (its values are the defaults, read by `nextflow.config`), plus
  `--max_cpus`, `--max_memory`, `--publish_mode` (`link`: hard links from `work/`, the default; `copy` when `work/`
  is on another disk), `--python` (runs `samples.py` at start-up) and `--conda_prefix`. `--steps` takes the same
  names and aliases as with Snakemake and runs only what they need; `--steps summary` alone makes the report of an
  existing results directory.
- **Resume**: `-resume` reuses every finished task (same inputs and parameters) from `work/`; run the same command
  again after an error or an interruption, from the same directory (`wgs_skim` runs Nextflow in the results
  directory). Keep `work/` until the batch is done: results are hard links into it, so deleting `work/` afterwards
  leaves them intact. Databases (GetOrganelle `animal_mt`, MitoM at `mitom_commit` / `mitom_nuc_commit`, MITOS2
  refseq89m, BUSCO lineage) are made once under `db_dir` (`storeDir`) and shared by all runs, Snakemake or Nextflow.
- **Software**: profile `conda` (`conda.enabled`). Environments that `wgs_skim install` (Snakemake) made in the
  conda prefix (`WGS_SKIM_CONDA_PREFIX` / `--conda_prefix`) from an identical `workflow/envs/*.yaml` are used as
  they are; otherwise Nextflow builds them in `conda.cacheDir` (= `WGS_SKIM_CONDA_PREFIX` when set).
- **Resources**: every process has the Snakemake rule's `threads` as `cpus` and its `mem_mb` as `memory`;
  SPAdes and GetOrganelle get theirs from the clean data in `fastp.json` (as in `common.smk`). The local executor
  runs tasks side by side within `--max_cpus` / `--max_memory` (default: all cores, 90 % of the memory); larger
  requests are capped (`resourceLimits`). A failed task is tried once more (SPAdes with 1.5x memory, GetOrganelle
  with 2x), then running tasks finish and the run stops. Per-process threads, as `--set-threads`: a config file
  `-c my.config` with `process { withName: SPADES { cpus = 48 } }`.
- **Slurm**: `-profile conda,slurm` submits every task as a job (`--slurm_queue`, `--slurm_options '--account=x'`,
  up to 50 jobs queued); `--max_cpus` / `--max_memory` then cap a single job.
- **Run reports**: `results/pipeline_info/` (trace with the peak memory and time of every task, HTML report,
  timeline), in place of Snakemake's `benchmarks/`.
- **Scripts**: the scripts written for Snakemake's `script:` directive (`mt_iterate.py`, `mt_select.py`,
  `topology.py`, `blast_taxonomy.py`, `summary.py`, `nuc_copies.py`) run unchanged through `nextflow/bin/smk_shim.py`,
  which gives them the `snakemake` object (input / output / params / log / threads / wildcards) from a JSON file.
- **Test**: `nextflow run nextflow/main.nf -stub-run --reads_dir DIR --outdir x --db_dir x/db` builds the whole DAG
  in seconds without running a tool; `-profile test` sets the short test settings of `tests/e2e.sh`.

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
   used (their extra ends were collapsed repeats or unsupported in tests), only reported. Circularization needs no
   reads (see below).
5. MitoM annotates the result (cox1 start, NCBI checks).
6. `topology` (default `auto`: only linear mitogenomes; `yes`: everything, also the rDNA repeat unit; `no`) checks
   with the reads: the circle is cut inside *nad5* and reads spanning the original ends must exist for a real
   circle (Ahuja et al. 2024, GBE). For linear results `end_diagnosis` says why each end stops — e.g. "reads
   continue into a dinucleotide repeat, in different directions: repeat shared with other loci", "coverage drops",
   "homopolymer" — which tells you whether long-range PCR / long reads are needed.
7. BLAST (species check) uses the queries you choose, several at once: whole sequences (`full`) and/or genes cut
   from the annotations — mitogenome genes from MitoM (default `cox1`, the COI barcode; `16s`, `12s`, `cytb` ...),
   rDNA regions from `mitom nuc` (default `18S,ITS,28S`). Hits come with GenBank organism and taxonomy.

### Circular or linear?

| evidence | needs | shows |
|---|---|---|
| ends overlap | sequence | the assembler *treated* it as a circle |
| path closes in the assembly graph | SPAdes graph | stronger: the graph itself closes |
| reads span the junction | short reads | a real circle (`topology`) |
| no reads span the junction | short reads | only "not joined": linear **or** a gap longer than the insert size |
| long reads ending at the same place, terminal inverted repeats / telomere-like ends, PCR across the ends fails while internal PCR works | long reads / lab | a truly linear molecule |

Short reads can never prove a linear mitogenome. All annelid mitogenomes known so far are circular; a linear annelid
assembly almost always means an unfinished control region (see `end_diagnosis`).

**Circularization without reads** (`wgs_skim circularize`, also used inside the workflow for SPAdes): MitoFinder and
MitoHiFi accept a circle only for a self-BLAST overlap >= 220 bp within 40 bp of both ends, but short-read assemblers
leave ~k-1 bp, often with ragged ends — circular contigs are then reported as "not found". `circularize` checks, in
order:

1. **graph** — the contig's path in the SPAdes assembly graph links its end back to its start;
2. **walk** — otherwise a short walk through the graph (≤ 12 segments, each at mt-level depth) closes the gap
   SPAdes left, usually in the control region. Tandem repeat units get the copy number their depth implies
   (segment depth / mt depth); copy numbers the path and depth disagree on are reported as uncertain;
3. **overlap** — no graph: an end overlap ≥ 25 bp tolerant of ragged ends (after the graph, because an overlap
   inside a tandem repeat collapses copies);
4. **join** — several pieces joined end to end.

```bash
wgs_skim circularize -i mitofinder/<sample>_all_mt.fasta --spades spades_dir/ -o circ/ --prefix <sample>
```

MJ20260820041, three individuals of one species, MitoFinder contigs on full-data SPAdes (~2 min, ~7 GB RAM each):

| sample | MitoFinder | circularize | vs. GetOrganelle (reads) |
|---|---|---|---|
| CHY006 | not found | walk, 9 segments + 2 turns of a 94-bp repeat (2500x vs 858x): 17,052 bp | identical (100%) |
| CHY007 | trimmed 17,040 of a 17,305 bp scaffold | scaffold closes: 17,228 bp; repeat unit used 4-5x at ~1x depth: copies uncertain | GetOrganelle only linear 16,392 bp |
| CHY008 | circular with a gap | walk, 3 segments (404-601x): 17,048 bp | 99.99% (one 2-bp indel) |

Without the graph the end overlap gave CHY006 16,770 bp (58-bp overlap inside the repeat: 3 copies collapsed).
Exact tandem copy numbers need long reads.

**Comparing assemblies** (`wgs_skim compare`): are GetOrganelle and SPAdes-MitoFinder (or any other assembly) the
same mitogenome? GetOrganelle is the reference because its graph shows where short reads cannot settle the
sequence (segments used several times in the path or with several times the depth: tandem repeats, multi-copy
regions). The other assembly is circularized first if needed (as above), turned to the same strand and start,
aligned (MAFFT), and every difference is placed on a gene and a GetOrganelle-flagged region.

```bash
wgs_skim compare -a getorganelle_out/ -b mitofinder=<sample>_all_mt.fasta --spades spades_dir/ \
                 --gb <sample>.gb -o compare/ --prefix <sample>      # -b can be repeated (label=fasta)
```

Verdict: **identical**; **consistent** — every difference is in a flagged repeat / gap, at the end of a linear
reference, or a ≤ 2-bp homopolymer length difference; **differ** otherwise. Outputs `<sample>.compare.tsv` (one row
per comparison: lengths, topology, identity, blocks, verdict, GetOrganelle problem regions), `<sample>.diffs.tsv`
(each difference: position, type, gene, graph region) and the alignment. MJ20260820041 (~2 min per sample):

| sample | GetOrganelle | MitoFinder after circularize | verdict | differences |
|---|---|---|---|---|
| CHY006 | 17,052 circular | 17,052 circular | identical | — |
| CHY007 | 16,318 linear | 17,228 circular | consistent | 910 bp beyond the linear end (control region, after trnH) |
| CHY008 | 17,050 circular | 17,048 circular | consistent | two 1-bp deletions in poly-T runs of 14 and 11 bp (control region) |

Structure and genes (QUAST-style, without QUAST): every contig of every assembly is also aligned locally to the
reference with minimap2 and QUAST 5.2's own options (`-r 200`: any jump over 200 bp ends an alignment), so
fragmented assemblies, moved blocks and inversions are seen. Taken from QUAST's source:
misassembly classes between consecutive alignments of a contig (other strand = inversion; reference/contig distance
mismatch > 1000 bp relocation, > 200 bp local), genome fraction, duplication ratio, unaligned length, complete /
partial genes. Changed for mitogenomes: circular reference handled whole (doubled), genes judged on the bases present
(complete >= 90 %, partial >= 10 %, missing below; QUAST's 100-bp partial threshold would hide every tRNA), each
junction placed against GetOrganelle's repeat regions, and each gene gets its base differences and, for protein
genes, whether the indels break the reading frame (`<sample>.genes.tsv`; junctions, long indels, reference parts not
in the assembly and unaligned contig parts in `<sample>.structure.tsv`). Verdict **incomplete**: no wrong sequence,
but genes or > 1 % of the reference missing. Annotation: `--gb`, or `--gff` on the reference sequence (as QUAST
`-g`; Geneious / MitoFinder GFF cleaned as the user's Quast/modify_gff2.py did: gene rows dropped, `Name=COX1 gene`
-> COX1, `trnN%28guu%29` -> trnN). Tools: mafft + minimap2 in `<conda_prefix>/wgs_skim_compare`, created on first
use (69 MB; minimap2's k8 dependency, only for paftools.js, would add a 470 MB sysroot).

Batch: every sample of a results directory at once (every assembly found is compared; samples without any are
listed as skipped), or a sample sheet for any other layout:

```bash
wgs_skim compare --results results/ -o compare/ [--samples A,B] [-j 2]   # wgs_skim_pipeline layout: GetOrganelle
                                     # 02.organelle/02_mt vs MitoFinder 10.mitofinder (+ SPAdes graph)
wgs_skim compare --sheet samples.tsv -o compare/    # TSV: sample, ref, other ([label=]fasta,...), [spades, gb, gff, gfa]
```

-> `compare/all.compare.tsv`, `all.diffs.tsv`, `all.structure.tsv`, `all.genes.tsv` and `all.compare.html` (table +
every figure). Each comparison gets `<sample>.<label>.compare.svg`: genes of the reference annotation on both
assemblies (carried over through the alignments), contigs in reference order, alignment links (orange = inverted),
differences and misassembly junctions red, GetOrganelle repeat regions purple, reference parts missing from the
assembly light red; hover for details. Plain SVG, no extra software (pyGenomeViz draws similar figures but needs a
389 MB matplotlib / numpy / Biopython environment). `--no-plot` skips the figures.

MJ20260820041, 3 at a time: 105 s, 7.2 GB peak memory (SPAdes graphs). MEANGS rows are from a test of
fragmented input only: MEANGS just gives GetOrganelle its seed and is not compared in batch mode.

| sample | assembly | verdict | genome fraction | genes complete / partial / missing | notes |
|---|---|---|---|---|---|
| CHY006 | MitoFinder | identical | 100 % | 36 / 0 / 0 | |
| CHY006 | MEANGS (4 contigs) | incomplete | 79.8 % | 32 / 2 / 2 | rrnL, trnV missing; no wrong base |
| CHY007 | MitoFinder | consistent | 100 % | 36 / 0 / 0 | 910 bp beyond the linear GetOrganelle end |
| CHY007 | MEANGS (11 contigs) | incomplete | 60.1 % | 7 / 14 / 15 | no wrong base |
| CHY008 | MitoFinder | consistent | 100 % | 36 / 0 / 0 | two poly-T lengths (control region) |
| CHY008 | MEANGS (4 contigs) | differ | 76.9 % | 27 / 4 / 5 | one mismatch in ND2 (8741) |

Synthetic checks (CHY006 reference): rotated / reverse-complemented circle - no event; 1-bp deletion in COX1 -
COX1 frameshift; ND4 inverted - 2 inversion junctions, ND4 reverse; ATP8 deleted - ATP8 missing, 182-bp deletion;
2-kb block moved - 3 relocations (+2000, -6000, +4000); 3 linear fragments - COX1 split, rrnL partial, two
reference gaps.

### Multi-copy nuclear markers

`nuc_markers` (default on, needs the SPAdes assembly) annotates rDNA (18S, ITS1, 5.8S, ITS2, 28S) and histones
(H1, H2A, H2B, H3, H4) in the whole-genome assembly with `mitom nuc`: every copy, one consensus per gene, and the
copy number = coverage of the gene's contigs / coverage of BUSCO single-copy contigs.

Two ways for each marker, compared (rules in `workflow/rules/nuclear.smk`, code in `workflow/scripts/nuclear_markers.py`):

- **GetOrganelle** assembles the tandem unit from the reads (`-F anonym`): rDNA with `resources/rdna`, histones with
  `resources/histone` (H1-H4 cluster copies of annelids + Hesionidae H3; an H3-only label kept only the H3 piece).
  MitoM annotates it: `rdna/<sample>/<sample>.rdna.gbf`, `histone/<sample>/<sample>.histone.gbf` (GenBank).
- **SPAdes-MitoM**: the copies scattered over the SPAdes contigs, one consensus per gene (`nuclear/<sample>/`).
- `nuclear/<sample>/<sample>.nuc_compare.tsv`: per gene both lengths, partial ends, histone class, mismatches and
  indels on the overlap, verdict identical / contained (the partial one inside the other) / differ / only one side.

Histone copies are classed by MitoM: replication-dependent (the cluster; H3 residues 31/87/89/90 = A S V M),
replacement (H3.3 = S A I G, or spliced) or single-copy (coverage < 2x the genome); each gene's consensus takes a
cluster copy first (MJ20260820041: SPAdes splits the cluster into 84-134 nt contigs, and the old rule took H3.3).

For an existing results directory (this workflow's or Temp_scripts/wgs_skim_pipeline's), the same steps in one
command, resumable (finished assemblies reused; annotations redone when MitoM changes):

```bash
wgs_skim nuclear --results results/ -o nuclear_markers/ [--samples A,B] [--markers rdna,histone] [-j N] [--max-gb 5,10]
wgs_skim nuclear --status -o nuclear_markers/        # what a running / stopped batch has done, per sample
```

-> `nuclear_markers/{rdna,histone}/{getorganelle,annotation}/<sample>/`, `nuclear_markers/spades_mitom/<sample>/`
(MitoM on the SPAdes scaffolds, both markers in one run), `summary.tsv` (every sample x
gene; rewritten after every finished sample, so a stopped batch keeps its results and a second run reuses them),
`README.txt`. Each step prints a time-stamped line as it starts and ends. `-j` (samples at a time) comes from free
memory and cores (1-4; ~12 GB and 16 threads each) unless given. This command is for results made outside the
workflow (e.g. Temp_scripts/wgs_skim_pipeline); within `wgs_skim run` the same steps are Snakemake rules
(`workflow/rules/nuclear.smk`) and scheduled with everything else. Reads: multi-copy units need little data - round 1 at most 5 Gb; only when
the unit is not complete, round 2 with 10 Gb seeded with round 1 (`--max-gb`, config `nuclear_max_gb`;
`<marker>/getorganelle/<sample>/rounds.tsv` lists the rounds with the label coverage of each gene). Complete = the
path is circular, or every core gene (18S, 5.8S, 28S; H2A, H2B, H3, H4) is covered end to end by a label sequence:
a linear rDNA unit with all genes whole is complete and gets no second round.
An rDNA assembly the results directory already has counts as round 1.

## Output

```
wgs_skim_results/
├── summary/samples.tsv        one row per sample: clean data, mitogenome length / circular / round, genes,
│   summary/report.html        rDNA length, topology verdicts, SPAdes N50, BUSCO, genome size, heterozygosity, top BLAST hits
├── qc/<sample>/               clean reads, fastp.html / fastp.json
├── mito/<sample>/<sample>.mt.fasta   final mitogenome (+ .mt.json: source, circular, why iteration stopped); round0..N
├── rdna/<sample>/<sample>.nr.fasta   18S-ITS-28S (+ .rdna.gbf / .nr_features.tsv: 18S ITS1 5.8S ITS2 28S, mitom nuc)
├── histone/<sample>/<sample>.histone.fasta   histone cluster (+ .histone.gbf: H1 H2A H2B H3 H4, mitom nuc)
├── annotation/<sample>/<sample>.mt.gb   MitoM annotation (GenBank; + .mitom.tsv gene counts)
├── topology/<sample>/         read-based circular / linear check, end diagnosis
├── nuclear/<sample>/          rDNA + histone copies in SPAdes, consensus, copy_number.tsv, nuc_compare.tsv
├── assembly/<sample>/         SPAdes scaffolds, BUSCO short summary; other SPAdes modes <sample>.<mode>.fasta,
│                              Trinity <sample>.trinity.fasta, each with <sample>.<mode>.busco.txt
├── genome_size/<sample>/      k-mer histogram, GenomeScope2 plots and summary
├── blast/blast_topn.tsv       COX1 / rDNA top hits with GenBank organism and taxonomy
└── logs/  benchmarks/         one log and one runtime/memory record per job, in directories named by rule
                               (logs/mitom, logs/getorganelle_mt + mt_iterate + mt_select, logs/getorganelle_nr ...)
```

### From the earlier scripts (Temp_scripts)

| earlier script / repository | here |
|---|---|
| `use_fastp_cleandata*.py`, `Fastp_Fastqc_mulitqc` | step `fastp` |
| `Batch_assembly_of_genomes_using_GetOrganelle*.py`, `batch_getorganelle`, `batch_seeds*.py` (MEANGS → GetOrganelle) | steps `meangs`, `mito`, `getorganelle` (+ reseeded rounds, SPAdes scaffold fallback) |
| `batch_obation_nr_using_getorg.py` (18S-28S), Hesionidae H3 seeds | step `nr` (rDNA + histone cluster, MitoM `.gbf`) |
| `spadesv3*.py` (wgs / rna / meta / bio / sewage / plasmid / metaviral / metaplasmid, --careful, --only-assembler, -k), `fastp_spades_busco.py`, `busco2summary` | steps `spades` + `busco` with `--spades-modes` (every SPAdes mode, one assembly + BUSCO each); BUSCO line in the report |
| `trinity_batch3*.py`, rnaspades snakefile, snakemake_transcriptome_workflow | steps `trinity`, `rna` (rnaSPAdes + Trinity + BUSCO transcriptome); TransDecoder / annotation not included |
| `jeyllfish_genome*.py`, `genome_size.py` | step `genome_size` |
| `Mitoz-annotate` / `mitoz_get*.py`, `batch_mitos2.py` | step `annotate` (MitoM: protein genes, rRNA, tRNA, NCBI checks) |
| `mt_circularize.py`, `check_circularity.pl`, MITObim `circules.py` | step `topology` (reads) and `wgs_skim circularize` (graph, no reads) |
| `blast2metadata` | step `blast` (local nt or NCBI online, GenBank taxonomy); `wgs_skim build-blastdb` |
| `mitofinder-batch`, `spades_mitofinder.py`, `mitofinder_annotate2summary.py` | `wgs_skim compare` reads MitoFinder results and compares them with GetOrganelle; MitoFinder itself is not run (MitoM annotates) |
| `batch_report.py` / `report_pdf.py` / `report_web.py` | `summary/samples.tsv` + `summary/report.html` (no PDF) |
| `using_bwa_draw_coverage*.py`, `organelle_coverage.py`, `plot_m2_depth.py` (depth plots) | not included yet (MitoM `mito_plot.py` draws genome maps with a depth ring) |
| `remove_getorganelle_templefile*.py`, `remove_spades_temp.sh` | automatic: GetOrganelle intermediates and the SPAdes K* directories are removed when a step finishes |

## Resuming and re-running

Snakemake keeps track of every output file: run the same command again after an interruption or an error and only
the missing or out-of-date jobs run. Useful Snakemake options (`wgs_skim run ... -- <option>`):
`--forcerun mt_iterate` (redo a step), `--until busco` (stop after a step), `--unlock` (after a crash),
`--rerun-triggers mtime` (ignore changed parameters).
With `--engine nextflow`, `-resume` is always on: the same command reuses every finished task (see "With Nextflow").

## Troubleshooting

Where to look: `logs/<step>/<sample>.log` is the tool's own output for that job; `benchmarks/<step>/<sample>.tsv`
its runtime and peak memory; `.snakemake/log/` the Snakemake log of each run. A failed job keeps `keep-going`: the
other samples continue, the summary waits. `wgs_skim run ... -n` shows what would run without running it.
With Nextflow: the error message names the task's work directory (`.command.sh` = the command, `.command.err`,
`logs/<step>/<sample>.log` inside it); `results/.nextflow.log` is the run log. Fix the cause, run the same command:
finished tasks are reused.

| message / symptom | cause | what to do |
|---|---|---|
| `no paired FASTQ files ... found` | file names not recognised | `wgs_skim samples -r DIR` shows what is detected; else a sample sheet (`-s`) |
| `sample 'X' found twice` | two lanes / two name patterns of one sample | merge the lanes, or a sample sheet |
| `Directory cannot be locked` / `IncompleteFilesException` | a previous run was interrupted | `wgs_skim run ... -- --unlock`, then the same command again (incomplete files are redone) |
| `CreateCondaEnvironmentException`, `wgs_skim install` stops | network / channel problem while building an environment | run `wgs_skim install` again (finished environments are kept); `--conda-prefix` on a disk with > 20 GB free |
| GetOrganelle: `Too few seed reads found` | the seed is too far from the sample (other family), or too few reads | histones map seeds locally already; give your own seed (`histone_seed`, `nr_seed`), or `nuclear_max_gb: 10,20` |
| GetOrganelle: `No space left on device` | it unzips the reads it uses into its output directory (up to 1.5x the data) | results directory on a large disk; the data caps (`nuclear_max_gb`, `mt_iter_max_data_gb`) bound it |
| mitogenome not circular (`mt_circular: no`) | reseeded rounds did not close it: repeat in the control region, low coverage | `topology/` says whether the ends can be joined (end diagnosis); `mt_iter_max_data_gb: 0` uses all reads; `wgs_skim compare` with a MitoFinder result |
| SPAdes: `<jemalloc>: ... out of memory` or killed | the memory estimate (`spades_gib_per_gb` x clean data) was too low | the retry (profile `retries: 1`) gives 1.5x and resumes from the last checkpoint; raise `spades_min_mem_gb` |
| SPAdes not run (`clean data ... < spades_min_clean_gb`) | too little data for an assembly | lower `spades_min_clean_gb` (0 = always) |
| BUSCO: `lineage ... not found` / download fails | no network, or wrong name | `--busco-lineage /path/to/metazoa_odb12` (a downloaded lineage directory) |
| BLAST takes hours | full nt is read for every search | `wgs_skim build-blastdb` once (Metazoa, <= 50 kb records: 68 GB, 1-2 min per search), then `--blast-db` |
| `organism: (efetch failed)` in `blast/blast_topn.tsv` | NCBI Entrez did not answer for those hits | the hits and identities are in the table; `-- --forcerun blast` later repeats the search and the metadata |
| `mitom_db`: `git clone` fails | no network for the pinned MitoM checkout / MITOS2 references | `mitom_dir: ~/mybiosofts/Mitom` and `mitos_refdir:` in the config (no clone needed then) |
| a background run disappears | the machine ran out of memory for other jobs (the OOM killer / the session manager stops it) | start long runs with `nohup ... &` or in `tmux`; `--mem-gb` keeps this workflow within its share |
| `rdna_copies` / `histones` empty in the report | no SPAdes assembly or no BUSCO single-copy genes to normalise against | needs steps `spades` + `busco` (`nuclear` runs them) |

Still stuck: the full command Snakemake ran is the second line `wgs_skim run` prints; run it with `-n -p` to see
each job's shell command, and open an issue with `logs/<step>/<sample>.log` attached.

## Testing

```bash
python -m pytest -q tests/                  # unit tests of the scripts (seconds, no tools needed)
tests/e2e.sh small_reads/ e2e_out/          # whole workflow on 2+ small samples, every output checked per sample
```

What was tested, how, with which data, the results and the bugs found: [docs/TESTING.md](docs/TESTING.md).

The end-to-end check used before every release is one 0.6 Gb polychaete sample (HNP1sub, 2 M read pairs): every
step runs, `summary/samples.tsv` has one row, and the nuclear markers are compared (rDNA 5/5 genes identical
between GetOrganelle and SPAdes-MitoM). Make such a sample from any run with
`seqtk sample -s1 R1.fq.gz 2000000 | gzip > sub_R1.fq.gz` (and R2 with the same seed).

## Software

fastp, MEANGS, GetOrganelle, SPAdes, BUSCO, Jellyfish, GenomeScope2, [MitoM](https://github.com/deyuanyang92-dev/Mitom), BWA, SAMtools, BLAST+ and
Biopython (GenBank taxonomy of the BLAST hits). Versions are pinned in
`workflow/envs/*.yaml`. Please cite the tools you use, and

> Yang D, Zeng S, Wang Z, Zhang Y, Yang D, Glasby CJ, Hwang J-S, Cai L (2024) Molecular systematics of *Perinereis* …
> *Zoosystematics and Evolution* 100(4): 1297–1314. https://doi.org/10.3897/zse.100.127201

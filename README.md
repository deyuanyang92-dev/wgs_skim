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
wgs_skim nuclear --results results/ -o nuclear_markers/ [--samples A,B] [--markers rdna,histone] [-j 3] [--all-reads]
```

-> `nuclear_markers/{rdna,histone}/{getorganelle,annotation,spades_mitom}/<sample>/`, `summary.tsv` (every sample x
gene), `README.txt`. `--all-reads`: the histone GetOrganelle with all reads, -R 20, k 21-127 (the user's 2024 setting).

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

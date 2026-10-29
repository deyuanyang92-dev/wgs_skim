# 配置

所有参数都在 `config/config.yaml`（下表由它生成：`python docs/wiki/tools/config_page.py`）。三种方式修改，后者优先：

1. `~/.wgs_skim/config.yaml`：个人长期设置（数据库目录、`blast_db`、`email`、`mitom_dir`），每次运行自动读取；
2. 命令行选项：`--taxon`、`--steps`、`--spades-modes`、`--blast-db`、`--mem-gb` ……（`wgs_skim run -h`）；
3. `--config 键=值 键=值`：任何一个键，例如 `--config mt_iter_max_data_gb=0 nuclear_max_gb=5,10`。

线程数和内存不在这里设置：`-j`（总核数）、`--mem-gb`（总内存），每个任务的用量见 [资源与续跑](Resources-and-resume)。

## Input / output

| 键 | 默认 | 说明 |
|---|---|---|
| `reads_dir` | `""` | directory of paired FASTQ files (flat, or one sub-directory per sample); names detected automatically |
| `samples` | `""` | optional sample sheet instead of reads_dir: TSV/CSV with columns sample, fastq_1, fastq_2 |
| `outdir` | `results` | results directory (relative to the working directory) |
| `db_dir` | `""` | where GetOrganelle / BUSCO databases are downloaded once (empty = resources/ of the workflow) |

## Taxon preset

| 键 | 默认 | 说明 |
|---|---|---|
| `taxon` | `annelid` | annelid / mollusc / fish / crustacean：同时设定 MEANGS 类群、线粒体遗传密码、MitoM 类群、BUSCO 谱系（下面 4 个键可单独覆盖） |
| `genetic_code` | `""` | NCBI mitochondrial code (5 invertebrate, 2 vertebrate) |
| `meangs_class` | `""` | A-worms Arthropoda Bryozoa Chordata Echinodermata Mollusca Nematoda N-worms Porifera-sponges |
| `clade` | `""` | MitoM protein set: Chordata Arthropoda Echinodermata Annelida-segmented-worms Bryozoa Mollusca Nematoda Nemertea-ribbon-worms Porifera-sponges Animal (all animals: other phyla) |
| `mitos_refdir` | `""` | existing MITOS2 reference directory (contains refseq89m); empty = downloaded once |
| `mitom_dir` | `""` | your own MitoM checkout (e.g. ~/mybiosofts/Mitom) instead of the pinned versions below |
| `mitom_commit` | `"27264d2cf672c718deff6915317a7ee01852830b"` | MitoM version for mitogenome annotation |
| `mitom_nuc_commit` | `"6f75f73373ad5198c200b6da222a5975d7ee73ea"` | MitoM version with `mitom nuc` (rDNA + histones; histone classes, 2026-10-28) |
| `busco_lineage` | `""` | lineage name (downloaded once, e.g. metazoa_odb12) or path to a downloaded lineage directory |

## Steps

| 键 | 默认 | 说明 |
|---|---|---|
| `steps` | `""` | 只运行部分步骤（逗号分隔），所需的上游步骤自动加入；可用值见 [步骤](Steps)；空 = 下面打开的全部 |
| `meangs` | `true` | MEANGS seed for GetOrganelle (false: GetOrganelle's own animal_mt seeds) |
| `mito` | `true` | MEANGS seed -> GetOrganelle mitogenome (+ reseeded iterations if not circular) |
| `nr` | `true` | GetOrganelle nuclear rDNA (18S-ITS-28S) |
| `assembly` | `true` | SPAdes whole-genome assembly + BUSCO |
| `nuc_markers` | `true` | rDNA (18S ITS1 5.8S ITS2 28S) + histones (H1-H4) in the assembly, copy numbers (mitom nuc) |
| `genome_size` | `true` | jellyfish + GenomeScope2 |
| `annotate` | `true` | MitoM annotation of the final mitogenome |
| `topology` | `auto` | read-based check of circular / linear (junction-spanning reads): auto = only linear mitogenomes; yes = every sequence, also the rDNA (tandem repeat unit); no = never |
| `blast` | `true` | BLAST mitogenome and rDNA, top hits with GenBank taxonomy (species check) |

## Mitogenome

| 键 | 默认 | 说明 |
|---|---|---|
| `getorganelle_kmers` | `"21,45,65,85,105"` |  |
| `getorganelle_rounds` | `10` | GetOrganelle -R for the first assembly |
| `meangs_reads` | `2000000` | reads MEANGS uses to find a seed |
| `insert_size` | `350` |  |
| `mt_iter_rounds` | `3` | reseeded GetOrganelle rounds when the first assembly is not circular |
| `mt_iter_max_data_gb` | `20` | reseeded rounds use all clean reads up to this many Gb (0 = no limit) |
| `mt_from_spades` | `true` | use the SPAdes scaffold only when it closes the circle and GetOrganelle did not |
| `topology_pairs` | `10000000` | read pairs used for the circular / linear check (0 = all) |

## 18S-28S rDNA

| 键 | 默认 | 说明 |
|---|---|---|
| `nr_seed` | `""` | seed FASTA; empty = bundled resources/rdna (polychaete 18S/5.8S/28S, works across animals) |
| `nr_label` | `""` | label FASTA (--genes); empty = bundled |
| `histone_seed` | `""` | histone cluster seed FASTA; empty = bundled resources/histone (annelid H1-H4 cluster copies) |
| `histone_label` | `""` | its label FASTA (H1-H4 genes); empty = bundled |
| `nuclear_max_gb` | `"5,10"` | Gb of reads per GetOrganelle round for rDNA / histones: round 2 (10 Gb, seeded with round 1) only when round 1 does not close the unit |

## Whole-genome assembly

| 键 | 默认 | 说明 |
|---|---|---|
| `spades_min_clean_gb` | `0` | skip SPAdes for samples with less clean data than this (0 = assemble all) |
| `spades_gib_per_gb` | `6.5` | SPAdes memory estimate: GiB per Gb of clean data (annelids 5.3–6.1) |
| `spades_min_mem_gb` | `32` |  |

## SPAdes modes and Trinity (the user's spadesv3.5.5.py / trinity_batch3.1.py)

| 键 | 默认 | 说明 |
|---|---|---|
| `spades_modes` | `"genome"` | comma separated: genome (whole-genome skim: BUSCO, nuclear markers, mt fallback) and any SPAdes mode: isolate sc meta rna rnaviral plasmid metaplasmid metaviral bio sewage corona -> assembly/<s>/<s>.<mode>.fasta + .<mode>.busco.txt |
| `spades_kmers` | `""` | -k for the non-genome modes ("" = SPAdes chooses) |
| `rnaspades_kmers` | `"21,33,55"` | -k of rna / rnaviral (the user's rnaspades snakefile); "auto" = SPAdes chooses |
| `spades_careful` | `false` | --careful (isolate, sc and plasmid modes only) |
| `spades_only_assembler` | `false` | --only-assembler for the non-genome modes (the genome mode always uses it) |
| `spades_mode_gib_per_gb` | `5` | memory estimate of the non-genome modes: GiB per Gb of clean data (>= 16 GB) |
| `trinity_gib_per_gb` | `6` | Trinity --max_memory from the clean data (>= 16 GB) |
| `trinity_extra` | `""` | more Trinity options, e.g. "--min_contig_length 300" |

## Genome size

| 键 | 默认 | 说明 |
|---|---|---|
| `kmer_size` | `21` |  |
| `ploidy` | `2` |  |

## BLAST

| 键 | 默认 | 说明 |
|---|---|---|
| `blast_db` | `""` | local nt database prefix (fast); empty = NCBI online BLAST (slow, no download) |
| `blast_mt` | `"cox1"` | mitogenome queries, comma separated: full (whole sequence) and/or genes from the MitoM annotation: cox1 (COI) cox2 cox3 cytb nad1-6 nad4l atp6 atp8 16s (rrnL) 12s (rrnS) |
| `blast_nr` | `"18S,ITS,28S"` | rDNA queries: full and/or 18S ITS (ITS1-5.8S-ITS2) ITS1 5.8S ITS2 28S (mitom nuc annotation) |
| `blast_topn` | `5` |  |
| `email` | `""` | required by NCBI Entrez for GenBank metadata, e.g. you@example.com |

# Usage: commands, configuration, presets

wgs_skim runs the same workflow in three ways. All three read the same options (the keys of
`config/config.yaml`; `wgs_skim config` prints them with comments).

| Way | When | Example |
|---|---|---|
| `wgs_skim` command | everyday use: builds the Snakemake / Nextflow command for you | `wgs_skim run -r rawdata/ -o results/ --preset full` |
| Snakemake directly | you want every Snakemake option, a cluster profile, your own config file | `snakemake -s workflow/Snakefile --configfile my.yaml ...` |
| Nextflow directly | you prefer Nextflow (`-resume`, Tower, cloud executors) | `nextflow run nextflow/main.nf -params-file my.yaml ...` |

## 1. One run, everything (raw data from the sequencing company)

```bash
wgs_skim run -r /data/MJ2026_rawdata -o /data/MJ2026_results --preset full -j 64
```
`full` = QC (once) -> mitogenome (MEANGS seed -> GetOrganelle, reseeded with the previous result while not
circular) -> MitoM annotation -> read-based circular / linear check -> nuclear rDNA (18S-ITS1-5.8S-ITS2-28S) and
histone cluster (GetOrganelle) -> SPAdes + BUSCO -> rDNA / histone copies in SPAdes -> genome size (GenomeScope) ->
SPAdes-MitoFinder -> BLAST species check -> `summary/samples.tsv` and `summary/report.html`.

Add `-n` first to see what would run. The run is resumable: after an interruption, run the same command again.
Results: one directory per tool, one sub-directory per sample, the tool's summary table next to them
(see [Output](Output.md)).

## 2. Presets: your own habits

A preset is a small YAML file of options. Built-in presets (`wgs_skim presets` lists them):

| Preset | Runs |
|---|---|
| `full` | everything above, including SPAdes-MitoFinder |
| `mito` | mitogenome only: MEANGS + GetOrganelle (reseeded), MitoM, topology check, BLAST cox1 |
| `nuclear` | rDNA unit and histone cluster (GetOrganelle), annotation, BLAST 18S / ITS / 28S |
| `genome` | QC, SPAdes, BUSCO, genome size |
| `transcriptome` | RNA-seq: rnaSPAdes + BUSCO, mitochondrial genes with MitoFinder, contamination check |

Your own preset: `~/.wgs_skim/presets/<name>.yaml`, e.g. `~/.wgs_skim/presets/my_skim.yaml`:
```yaml
# my_skim: mitogenome with 5 reseeded rounds, rDNA + histones, BLAST on the local database
steps: "mito,annotate,nr,blast,summary"
mt_iter_rounds: 5
blast_db: /db/blast/metazoa_marker
```
then `wgs_skim run -r rawdata/ -o results/ --preset my_skim`.

**Without writing YAML**: `wgs_skim gui` opens a page where you tick the functions and fill in the options; it
saves the preset and shows the command (and can preview the jobs with a dry run). On a server:
`ssh -L 8765:localhost:8765 server`, run `wgs_skim gui --no-browser` there, open http://localhost:8765.

**Priority** (later wins): `config/config.yaml` (defaults) < `~/.wgs_skim/config.yaml` (your machine: database
paths, e-mail, MitoFinder image ...) < `--preset` < `--configfile` < command-line options (`--steps`,
`--mt-rounds`, `--config KEY=VALUE`).

### Mitogenome iteration

MEANGS gives the seed, GetOrganelle assembles (round 0); when the result is not circular, GetOrganelle runs again
with the previous result as the seed and all clean reads (round 1, 2, ...). It stops at the first circular result,
when a round gives the same sequence as the one before, or after `mt_iter_rounds` rounds (default 3).

| You want | Command |
|---|---|
| MEANGS + GetOrganelle, reseeded while not circular | `wgs_skim run -r raw/ --steps meangs-getorganelle` |
| more rounds | `... --steps meangs-getorganelle --mt-rounds 6` |
| one GetOrganelle run, no reseeding | `... --mt-rounds 0` |
| GetOrganelle's own seeds instead of MEANGS | `wgs_skim getorganelle -r clean/` (or `--config meangs=false`) |

## 3. Every function: command, Snakemake, Nextflow

`wgs_skim <tool>` runs one tool for every sample (reads tools: `-r READS_DIR`, no QC unless `--qc`; FASTA tools:
`-i FILES_OR_DIRS`). `--steps` takes the same names in a full run. Snakemake / Nextflow: `steps=` / `--steps`.

| Function | Tool command | Step (`--steps`) | Main options |
|---|---|---|---|
| QC | `wgs_skim fastp -r raw/` | `fastp` | `qc` |
| Mitochondrial seed | `wgs_skim meangs -r clean/` | `meangs` | `meangs_class`, `meangs_reads` |
| Mitogenome | `wgs_skim getorganelle -r clean/ [--meangs] [--rounds N]` | `mito` / `meangs-getorganelle`, `getorganelle` | `mt_iter_rounds`, `mt_iter_max_data_gb`, `getorganelle_kmers` |
| rDNA + histones | `wgs_skim getorganelle -r clean/ --target nr` | `nr` | `nuclear_max_gb`, `nr_seed`, `histone_seed` |
| SPAdes (any mode) | `wgs_skim spades -r clean/ [--mode genome,rna,meta]` | `spades`, `rnaspades` | `spades_modes`, `spades_kmers` |
| Trinity | `wgs_skim trinity -r clean/` | `trinity` | `trinity_extra` |
| BUSCO | `wgs_skim busco -i asm/ [-m transcriptome] [-l lineage]` | `busco` | `busco_lineage`, `busco_mode` |
| Genome size | `wgs_skim genome-size -r clean/` | `genome_size` | `kmer_size`, `ploidy` |
| MitoM annotation | `wgs_skim mitom -i mitogenomes/` | `annotate` | `clade`, `genetic_code` |
| MitoFinder | `wgs_skim mitofinder -i asm/ --ref refs/` or `-r clean/ --assembler megahit` | `mitofinder` (`spades-mitofinder`, `rnaspades-mitofinder`, `trinity-mitofinder`) | `mitofinder_ref`, `mitofinder_assembly`, `mitofinder_blast` |
| BLAST species check | `wgs_skim blast -i seqs/ [--genes cox1,16s]` | `blast` | `blast_db`, `blast_mt`, `blast_nr`, `email` |
| Nuclear-marker copies (SPAdes) | `wgs_skim nuclear --results results/` | `nuclear` | `nuclear_max_gb` |
| Circular / linear check | - | `topology` | `topology`, `topology_pairs` |
| Summary report | - | `summary` | - |

Combinations: `meangs-getorganelle`, `spades-mitofinder`, `spades-busco`, `fastp-spades-busco`,
`rnaspades-busco`, `trinity-busco`, `rnaspades-mitofinder`, `trinity-mitofinder` (see [Combinations](Combinations.md)).

## 4. Snakemake directly

From the wgs_skim directory (or give `-s /path/to/wgs_skim/workflow/Snakefile`):
```bash
snakemake -s workflow/Snakefile --sdm conda --cores 64 --resources mem_mb=400000 \
    --configfile ~/.wgs_skim/config.yaml config/presets/full.yaml \
    --config reads_dir=/data/raw outdir=. -d /data/results
```
- `--configfile A B`: several files in one option (later ones win). A repeated `--configfile` keeps only the last.
- `--config KEY=VALUE ...`: any option of `config/config.yaml`.
- `-d DIR`: the results directory (`outdir=.` writes there); `-n` dry run; `--rerun-incomplete` after a crash.
- Cluster: add a profile (`--profile slurm`) or `--executor slurm --default-resources slurm_account=...`.

Without a config file at all:
```bash
snakemake -s workflow/Snakefile --sdm conda --cores 32 -d results \
    --config reads_dir=/data/raw steps=meangs-getorganelle,annotate mt_iter_rounds=5
```

## 5. Nextflow directly

```bash
NXF_SYNTAX_PARSER=v1 nextflow run nextflow/main.nf -profile conda -resume \
    -params-file config/presets/genome.yaml --reads_dir /data/raw --outdir /data/results --max_cpus 64
```
- `-params-file`: one YAML file (a preset, or your own); options on the command line as `--key value` win.
- `-resume` continues after an interruption; `work/` holds the cache.
- Cluster: `-profile conda,slurm --slurm_queue normal`.
- Not in the Nextflow version yet (wgs_skim refuses instead of running half): MitoFinder, Trinity, SPAdes modes
  other than genome, transcriptome steps. Use Snakemake for those.

`wgs_skim run --engine nextflow --preset genome` builds this command for you (and merges
`~/.wgs_skim/config.yaml` with the preset into one params file).

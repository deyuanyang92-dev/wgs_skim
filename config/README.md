# Configuration

`config/config.yaml` holds every option with its default and a comment. Only `reads_dir` (or `samples`) is required:

```bash
snakemake --cores 32 --sdm conda --config reads_dir=/path/to/rawdata
```

- `reads_dir`: directory of paired FASTQ files; sample names are taken from the file names (see README, *Input*).
- `samples`: alternatively a sample sheet, TSV or CSV with columns `sample`, `fastq_1`, `fastq_2`.
- `taxon`: `annelid` (default), `mollusc`, `fish` or `crustacean` — sets MEANGS class, genetic code, MitoM clade and
  BUSCO lineage; override any of them with `meangs_class`, `genetic_code`, `clade`, `busco_lineage`.
- step switches: `mito`, `nr`, `assembly`, `genome_size`, `annotate`, `topology`, `blast` (true/false).
- `blast_db`: local nt prefix; empty uses NCBI online BLAST. `email` is needed for GenBank metadata.

Threads and memory are runtime settings: use `--cores`, `--set-threads RULE=N`, `--set-resources RULE:mem_mb=N` or a
profile, as for any Snakemake workflow. The options are validated against `workflow/schemas/config.schema.yaml`.

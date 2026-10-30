rule blast:
    """Species check of the batch: chosen mitogenome / rDNA queries (whole sequences or genes: blast_mt, blast_nr)
    against nt, TopN hits with GenBank organism and taxonomy. Local nt if blast_db is set, else NCBI online BLAST."""
    input:
        unpack(blast_inputs),
    output:
        tsv=f"{OUT}blast/blast_topn.tsv",
    params:
        db=config["blast_db"],
        topn=config["blast_topn"],
        email=config["email"],
        samples=list(SAMPLES),
        mt=BLAST_MT,
        nr=BLAST_NR,
        cache=str(DB / "genbank_metadata.tsv"),       # GenBank organism / taxonomy of hits, fetched once
    resources:
        mem_mb=48000,                     # a local nt / marker database is memory-mapped: 47.5 GB peak (metazoa_marker)
    log:
        f"{OUT}pipeline_info/logs/blast.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/blast.tsv"
    threads: 8
    conda:
        "../envs/blast.yaml"
    script:
        "../scripts/blast_taxonomy.py"

rule blast:
    """Species check of the batch: chosen mitogenome / rDNA queries (whole sequences or genes: blast_mt, blast_nr)
    against nt, TopN hits with GenBank organism and taxonomy. Local nt if blast_db is set, else NCBI online BLAST."""
    input:
        unpack(blast_inputs),
    output:
        tsv=f"{OUT}/blast/blast_topn.tsv",
    params:
        db=config["blast_db"],
        topn=config["blast_topn"],
        email=config["email"],
        samples=list(SAMPLES),
        mt=BLAST_MT,
        nr=BLAST_NR,
    log:
        f"{OUT}/logs/blast.log",
    benchmark:
        f"{OUT}/benchmarks/blast.tsv"
    threads: 8
    conda:
        "../envs/blast.yaml"
    script:
        "../scripts/blast_taxonomy.py"

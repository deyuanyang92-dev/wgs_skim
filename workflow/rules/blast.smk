rule blast:
    """All mitogenomes and rDNA of the batch -> blast2metadata (TopN hits + GenBank taxonomy). Local nt if
    blast_db is set, otherwise NCBI online BLAST."""
    input:
        blast_inputs,
    output:
        tsv=f"{OUT}/blast/blast_topn.tsv",
    params:
        out=lambda w, output: str(Path(output.tsv).parent / "blast2metadata"),
        topn=config["blast_topn"],
        mode=(f"-bm local --blast_db_nt {config['blast_db']}" if config["blast_db"]
              else "-bm online --online_engine biopython"),
        email=f"--email {config['email']}" if config["email"] else "",
    log:
        f"{OUT}/logs/blast.log",
    benchmark:
        f"{OUT}/benchmarks/blast.tsv"
    threads: 8
    conda:
        "../envs/blast.yaml"
    shell:
        """
        mkdir -p {params.out}
        cat {input} | sed -e 's/ .*//' > {params.out}/query.fasta
        blast2metadata -i {params.out}/query.fasta -o {params.out} --topN {params.topn} --threads {threads} \
            {params.mode} {params.email} --resume --no_smk > {log} 2>&1
        cp {params.out}/summary_results/query_summary.tsv {output.tsv}
        """

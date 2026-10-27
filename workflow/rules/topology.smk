rule topology:
    """Are circular / linear results real? Junction-spanning reads after splitting inside nad5 (Ahuja et al. 2024)."""
    input:
        mt=rules.mt_iterate.output.fasta,
        nr=lambda wc: rules.getorganelle_nr.output.fasta if config["nr"] else [],
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
        annotation=lambda wc: rules.mitoz_annotate.output.gb if config["annotate"] else [],
    output:
        tsv=f"{OUT}/topology/{{sample}}/{{sample}}.topology.tsv",
    params:
        out=lambda w, output: str(Path(output.tsv).parent),
        script=str(Path(workflow.basedir) / "scripts" / "mt_circularize.py"),
        annotation=lambda wc, input: f"--annotation {input.annotation}" if input.annotation else "",
    log:
        f"{OUT}/logs/topology/{{sample}}.log",
    benchmark:
        f"{OUT}/benchmarks/topology/{{sample}}.tsv"
    threads: 8
    conda:
        "../envs/topology.yaml"
    shell:
        """
        cat {input.mt} {input.nr} > {params.out}/input.fasta
        if [ ! -s {params.out}/input.fasta ]; then echo "nothing assembled" > {log}; : > {output.tsv}; exit 0; fi
        python {params.script} -i {params.out}/input.fasta -o {params.out} \
            --prefix {wildcards.sample} -1 {input.r1} -2 {input.r2} -t {threads} {params.annotation} > {log} 2>&1
        cp {params.out}/circularize.tsv {output.tsv}
        """

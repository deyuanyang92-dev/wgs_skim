rule topology:
    """Are circular / linear results real? Junction-spanning reads after cutting inside nad5 (Ahuja et al. 2024)."""
    input:
        fasta=lambda wc: [rules.mt_select.output.fasta.format(sample=wc.sample)]
        + ([rules.getorganelle_nr.output.fasta.format(sample=wc.sample)] if config["nr"] else []),
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
        annotation=lambda wc: rules.mitom_annotate.output.gb if config["annotate"] else [],
    output:
        tsv=f"{OUT}/topology/{{sample}}/{{sample}}.topology.tsv",
    params:
        pairs=config["topology_pairs"],
        mode=TOPOLOGY,
    log:
        f"{OUT}/logs/topology/{{sample}}.log",
    benchmark:
        f"{OUT}/benchmarks/topology/{{sample}}.tsv"
    threads: 8
    conda:
        "../envs/topology.yaml"
    script:
        "../scripts/topology.py"

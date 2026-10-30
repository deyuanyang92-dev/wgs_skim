rule topology:
    """Are circular / linear results real? Junction-spanning reads after cutting inside nad5 (Ahuja et al. 2024)."""
    input:
        unpack(clean_reads),
        fasta=lambda wc: [mt_of(wc.sample)]
        + ([rules.getorganelle_nr.output.fasta.format(sample=wc.sample)] if config["nr"] else []),
        annotation=lambda wc: rules.mitom_annotate.output.gb if config["annotate"] else [],
        status=lambda wc: mt_of(wc.sample, "status"),
    output:
        tsv=f"{OUT}topology/{{sample}}/{{sample}}.topology.tsv",
    params:
        pairs=config["topology_pairs"],
        mode=TOPOLOGY,
    resources:
        mem_mb=8000,
    log:
        f"{OUT}pipeline_info/logs/topology/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/topology/{{sample}}.tsv"
    threads: 8
    conda:
        "../envs/topology.yaml"
    script:
        "../scripts/topology.py"

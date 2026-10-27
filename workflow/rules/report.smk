rule summary:
    """One row per sample (samples.tsv) and a self-contained HTML report."""
    input:
        unpack(summary_inputs),
    output:
        tsv=f"{OUT}/summary/samples.tsv",
        html=f"{OUT}/summary/report.html",
    params:
        samples=list(SAMPLES),
        taxon=config["taxon"],
    log:
        f"{OUT}/logs/summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summary.py"

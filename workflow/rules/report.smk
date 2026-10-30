rule summary:
    """One row per sample (samples.tsv) and a self-contained HTML report."""
    input:
        unpack(summary_inputs),
    output:
        tsv=f"{OUT}summary/samples.tsv",
        html=f"{OUT}summary/report.html",
    params:
        samples=list(SAMPLES),
        taxon=config["taxon"],
    resources:
        mem_mb=1000,
    log:
        f"{OUT}pipeline_info/logs/summary.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/summary.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summary.py"


rule tool_summary:
    """one tool, one row per sample (and per assembly): <tool>/<tool>_summary.tsv (scripts/summary.py tool_table)"""
    input:
        lambda w: tool_summaries()[w.tool],
    output:
        tsv=f"{OUT}{{tool}}/{{tool}}_summary.tsv",
    resources:
        mem_mb=1000,
    log:
        f"{OUT}pipeline_info/logs/tool_summary/{{tool}}.log",
    wildcard_constraints:
        tool="fastp|meangs|getorganelle|mitogenome|spades|trinity|busco|genomescope|mitom|rdna|histone",
    localrule: True
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summary.py"

rule spades:
    """Whole-genome assembly. -m is a hard limit: estimated from clean data, x1.5 on each retry (--retries),
    resuming from the last SPAdes checkpoint instead of starting over."""
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
        json=rules.fastp.output.json,
    output:
        scaffolds=f"{OUT}/assembly/{{sample}}/{{sample}}.scaffolds.fasta",
    params:
        out=lambda w, output: str(Path(output.scaffolds).parent / "spades"),
        min_gb=config["spades_min_clean_gb"],
    resources:
        mem_mb=spades_mem_mb,
    log:
        f"{OUT}/logs/spades/{{sample}}.log",
    benchmark:
        f"{OUT}/benchmarks/spades/{{sample}}.tsv"
    threads: 32
    conda:
        "../envs/spades.yaml"
    shell:
        """
        gb=$(python -c "import json;print(json.load(open('{input.json}'))['summary']['after_filtering']['total_bases']/1e9)")
        if python -c "import sys;sys.exit(not ({params.min_gb} > 0 and $gb < {params.min_gb}))"; then
            echo "[wgs_skim] clean data $gb Gb < spades_min_clean_gb {params.min_gb}: not assembled" > {log}
            : > {output.scaffolds}; exit 0
        fi
        mem=$(( {resources.mem_mb} / 1024 ))
        if [ -f {params.out}/params.txt ]; then
            spades.py -o {params.out} --restart-from last -t {threads} -m $mem > {log} 2>&1
        else
            spades.py -1 {input.r1} -2 {input.r2} -o {params.out} -t {threads} -m $mem --only-assembler > {log} 2>&1
        fi
        cp {params.out}/scaffolds.fasta {output.scaffolds}
        """


rule busco_lineage:
    """Download a BUSCO lineage once (skipped when busco_lineage is a local directory)."""
    output:
        directory(str(DB / "busco_downloads" / "lineages" / "{lineage}")),
    params:
        path=lambda w, output: str(Path(output[0]).parent.parent),
    log:
        f"{OUT}/logs/busco_download/{{lineage}}.log",
    conda:
        "../envs/busco.yaml"
    shell:
        "busco --download {wildcards.lineage} --download_path {params.path} > {log} 2>&1"


rule busco:
    input:
        scaffolds=rules.spades.output.scaffolds,
        lineage=busco_lineage_dir,
    output:
        summary=f"{OUT}/assembly/{{sample}}/{{sample}}.busco.txt",
    params:
        out=lambda w, output: str(Path(output.summary).parent),
    log:
        f"{OUT}/logs/busco/{{sample}}.log",
    benchmark:
        f"{OUT}/benchmarks/busco/{{sample}}.tsv"
    threads: 16
    conda:
        "../envs/busco.yaml"
    shell:
        """
        if [ ! -s {input.scaffolds} ]; then echo "no assembly" > {output.summary}; exit 0; fi
        busco -i {input.scaffolds} -l {input.lineage} -m genome -o busco --out_path {params.out} \
            -c {threads} --offline -f > {log} 2>&1
        cp {params.out}/busco/short_summary.*.txt {output.summary}
        """

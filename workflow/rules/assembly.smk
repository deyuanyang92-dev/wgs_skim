rule spades:
    """Whole-genome assembly. -m is a hard limit: estimated from clean data, x1.5 on each retry (--retries),
    resuming from the last SPAdes checkpoint instead of starting over. Samples with an assembly made elsewhere
    (sheet column `scaffolds`) are linked by rule assembly_given instead (pipeline_info/inputs/assembly/)."""
    input:
        unpack(clean_reads),
    output:
        scaffolds=f"{OUT}spades/{{sample}}/{{sample}}.scaffolds.fasta",
    params:
        out=lambda w, output: str(Path(output.scaffolds).parent / "spades"),
        min_gb=config["spades_min_clean_gb"],
    resources:
        mem_mb=spades_mem_mb,
    log:
        f"{OUT}pipeline_info/logs/spades/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/spades/{{sample}}.tsv"
    threads: 32
    wildcard_constraints:
        sample=one_of([x for x in SAMPLES if x not in GIVEN_ASM]),
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
        # continue an interrupted run (killed, out of memory) only with the same reads and options; otherwise start
        # clean, so changed reads / options are never mixed into an old run (--restart-from keeps the old ones)
        sig="$(stat -L -c '%n %s %Y' {input.r1} {input.r2} | tr '\n' ' ') --only-assembler"
        if [ -f {params.out}/params.txt ] && [ "$(cat {params.out}/wgs_skim.sig 2>/dev/null)" = "$sig" ]; then
            spades.py -o {params.out} --restart-from last -t {threads} -m $mem > {log} 2>&1
        else
            rm -rf {params.out}; mkdir -p {params.out}; echo "$sig" > {params.out}/wgs_skim.sig
            spades.py -1 {input.r1} -2 {input.r2} -o {params.out} -t {threads} -m $mem --only-assembler > {log} 2>&1
        fi
        cp {params.out}/scaffolds.fasta {output.scaffolds}
        rm -rf {params.out}/K* {params.out}/tmp {params.out}/misc {params.out}/before_rr.fasta   # intermediates (~2 GB/Gb)
        rm -f {params.out}/wgs_skim.sig        # finished: a forced rerun starts clean
        """


rule assembly_given:
    """an assembly made elsewhere (sample sheet column scaffolds, wgs_skim busco / mitofinder -i): linked"""
    input:
        lambda w: GIVEN[w.sample]["scaffolds"],
    output:
        scaffolds=f"{OUT}pipeline_info/inputs/assembly/{{sample}}/{{sample}}.scaffolds.fasta",
    wildcard_constraints:
        sample=one_of(GIVEN_ASM),
    localrule: True
    shell:
        "ln -sf $(realpath {input}) {output.scaffolds}"


rule busco_lineage:
    """Download a BUSCO lineage once (skipped when busco_lineage is a local directory)."""
    output:
        directory(str(DB / "busco_downloads" / "lineages" / "{lineage}")),
    params:
        path=lambda w, output: str(Path(output[0]).parent.parent),
    resources:
        mem_mb=1000,
    log:
        f"{OUT}pipeline_info/logs/busco_download/{{lineage}}.log",
    conda:
        "../envs/busco.yaml"
    shell:
        "busco --download {wildcards.lineage} --download_path {params.path} > {log} 2>&1"


rule busco:
    input:
        scaffolds=lambda w: scaffolds_of(w.sample),
        lineage=busco_lineage_dir,
    output:
        summary=f"{OUT}busco/{{sample}}/{{sample}}.busco.txt",
    params:
        out=lambda w, output: str(Path(output.summary).parent),
        mode=config["busco_mode"],
    resources:
        mem_mb=16000,
    log:
        f"{OUT}pipeline_info/logs/busco/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/busco/{{sample}}.tsv"
    threads: 16
    conda:
        "../envs/busco.yaml"
    shell:
        """
        if [ ! -s {input.scaffolds} ]; then echo "no assembly" > {output.summary}; exit 0; fi
        rm -rf {params.out}/busco
        busco -i {input.scaffolds} -l {input.lineage} -m {params.mode} -o busco --out_path {params.out} \
            -c {threads} --offline -f > {log} 2>&1
        cp {params.out}/busco/short_summary.*.txt {output.summary}
        """

rule jellyfish:
    input:
        unpack(clean_reads),
    output:
        histo=f"{OUT}genomescope/{{sample}}/{{sample}}.k{config['kmer_size']}.histo",
    params:
        k=config["kmer_size"],
        jf=lambda w, output: str(Path(output.histo).parent / "reads.jf"),
    resources:
        mem_mb=12000,
    log:
        f"{OUT}pipeline_info/logs/jellyfish/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/jellyfish/{{sample}}.tsv"
    threads: 16
    conda:
        "../envs/genome_size.yaml"
    shell:
        """
        jellyfish count -C -m {params.k} -s 2G -t {threads} -o {params.jf} \
            <(zcat {input.r1}) <(zcat {input.r2}) > {log} 2>&1
        jellyfish histo -t {threads} {params.jf} > {output.histo} 2>> {log}
        rm -f {params.jf}
        """


rule genomescope:
    input:
        histo=rules.jellyfish.output.histo,
    output:
        summary=f"{OUT}genomescope/{{sample}}/{{sample}}.genomescope.txt",
    params:
        out=lambda w, output: str(Path(output.summary).parent / "genomescope"),
        k=config["kmer_size"],
        ploidy=config["ploidy"],
    resources:
        mem_mb=2000,
    log:
        f"{OUT}pipeline_info/logs/genomescope/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/genomescope/{{sample}}.tsv"
    conda:
        "../envs/genome_size.yaml"
    shell:
        """
        genomescope2 -i {input.histo} -o {params.out} -k {params.k} -p {params.ploidy} -n {wildcards.sample} > {log} 2>&1 \
            || echo "[wgs_skim] GenomeScope did not converge (low coverage?)" >> {log}
        if [ -s {params.out}/{wildcards.sample}_summary.txt ]; then cp {params.out}/{wildcards.sample}_summary.txt {output.summary}
        else echo "GenomeScope failed, see {log}" > {output.summary}; fi
        """

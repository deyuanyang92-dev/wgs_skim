rule jellyfish:
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
    output:
        histo=f"{OUT}genome_size/{{sample}}/{{sample}}.k{config['kmer_size']}.histo",
    params:
        k=config["kmer_size"],
        jf=lambda w, output: str(Path(output.histo).parent / "reads.jf"),
    resources:
        mem_mb=64000,
    log:
        f"{OUT}logs/jellyfish/{{sample}}.log",
    benchmark:
        f"{OUT}benchmarks/jellyfish/{{sample}}.tsv"
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
        summary=f"{OUT}genome_size/{{sample}}/{{sample}}.genomescope.txt",
    params:
        out=lambda w, output: str(Path(output.summary).parent / "genomescope"),
        k=config["kmer_size"],
        ploidy=config["ploidy"],
    log:
        f"{OUT}logs/genomescope/{{sample}}.log",
    conda:
        "../envs/genome_size.yaml"
    shell:
        """
        genomescope2 -i {input.histo} -o {params.out} -k {params.k} -p {params.ploidy} -n {wildcards.sample} > {log} 2>&1 \
            || echo "[wgs_skim] GenomeScope did not converge (low coverage?)" >> {log}
        if [ -s {params.out}/{wildcards.sample}_summary.txt ]; then cp {params.out}/{wildcards.sample}_summary.txt {output.summary}
        else echo "GenomeScope failed, see {log}" > {output.summary}; fi
        """

rule mitoz_annotate:
    """MitoZ annotation (PCGs, tRNAs, rRNAs) of the final mitogenome; GenBank file for topology and the report."""
    input:
        fasta=rules.mt_iterate.output.fasta,
    output:
        gb=f"{OUT}/annotation/{{sample}}/{{sample}}.mt.gb",
    params:
        wd=lambda w, output: str(Path(output.gb).parent),
        clade=config["mitoz_clade"],
        code=config["genetic_code"],
    log:
        f"{OUT}/logs/mitoz_annotate/{{sample}}.log",
    benchmark:
        f"{OUT}/benchmarks/mitoz_annotate/{{sample}}.tsv"
    threads: 8
    conda:
        "../envs/mitoz.yaml"
    shell:
        """
        if [ ! -s {input.fasta} ]; then echo "no mitogenome" > {log}; : > {output.gb}; exit 0; fi
        fa=$(readlink -f {input.fasta}); log=$(readlink -f {log})
        cd {params.wd}
        mitoz annotate --fastafiles $fa --outprefix {wildcards.sample} --thread_number {threads} \
            --clade {params.clade} --genetic_code {params.code} --species_name {wildcards.sample} > $log 2>&1
        gb=$(find {wildcards.sample}.result -maxdepth 1 -name '*.gbf' | head -n 1)
        cp "$gb" {wildcards.sample}.mt.gb
        """

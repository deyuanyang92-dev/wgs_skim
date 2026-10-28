MITOS_REF_URL = "https://zenodo.org/records/4284483/files/refseq89m.tar.bz2"
MITOS_REF_MD5 = "bb6325e27612e61a6995e88e8ffcecf8"


rule mitom_db:
    """MitoM code (pinned commit) and the MITOS2 reference data refseq89m it uses (skipped if mitos_refdir is set)."""
    output:
        directory(MITOM_DB),
    params:
        refdir=config["mitos_refdir"],
        commit=config["mitom_commit"],
        url=MITOS_REF_URL,
        md5=MITOS_REF_MD5,
    log:
        f"{OUT}logs/mitom_db.log",
    conda:
        "../envs/mitom.yaml"
    shell:
        """
        git clone -q https://github.com/deyuanyang92-dev/Mitom.git {output}/Mitom > {log} 2>&1
        git -C {output}/Mitom checkout -q {params.commit} >> {log} 2>&1
        if [ -z "{params.refdir}" ]; then
            mkdir -p {output}/mitos_refs
            wget -q -O {output}/mitos_refs/refseq89m.tar.bz2 {params.url} >> {log} 2>&1
            echo "{params.md5}  {output}/mitos_refs/refseq89m.tar.bz2" | md5sum -c >> {log} 2>&1
            tar -xjf {output}/mitos_refs/refseq89m.tar.bz2 -C {output}/mitos_refs && rm {output}/mitos_refs/refseq89m.tar.bz2
        fi
        """


rule mitom_nuc_db:
    """MitoM version with `mitom nuc` (nuclear rDNA 18S-ITS1-5.8S-ITS2-28S annotation; reference library included)."""
    output:
        directory(MITOM_NUC_DB),
    params:
        commit=config["mitom_nuc_commit"],
    log:
        f"{OUT}logs/mitom_nuc_db.log",
    conda:
        "../envs/mitom.yaml"
    shell:
        """
        git clone -q https://github.com/deyuanyang92-dev/Mitom.git {output}/Mitom > {log} 2>&1
        git -C {output}/Mitom checkout -q {params.commit} >> {log} 2>&1
        """


rule mitom_annotate:
    """MitoM annotation (protein genes, rRNAs, tRNAs; cox1 start; NCBI checks) of the final mitogenome."""
    input:
        fasta=rules.mt_select.output.fasta,
        db=MITOM_DB,
    output:
        gb=f"{OUT}annotation/{{sample}}/{{sample}}.mt.gb",
        tsv=f"{OUT}annotation/{{sample}}/{{sample}}.mitom.tsv",
    params:
        wd=lambda w, output: str(Path(output.gb).parent / "mitom"),
        refdir=lambda w, input: config["mitos_refdir"] or f"{input.db}/mitos_refs",
        mitom=lambda w, input: config["mitom_dir"] or f"{input.db}/Mitom",
        clade=config["clade"],
        code=config["genetic_code"],
    log:
        f"{OUT}logs/mitom/{{sample}}.log",
    benchmark:
        f"{OUT}benchmarks/mitom/{{sample}}.tsv"
    threads: 8
    conda:
        "../envs/mitom.yaml"
    shell:
        """
        if [ ! -s {input.fasta} ]; then echo "no mitogenome" > {log}; : > {output.gb}; : > {output.tsv}; exit 0; fi
        rm -rf {params.wd}; mkdir -p {params.wd}
        # record IDs "S|mt|1" -> "S_mt_1" (topology matches the annotation to the sequence by content)
        awk '{{if (substr($0, 1, 1) == ">") {{sub(" .*", ""); gsub("[^A-Za-z0-9_.>-]", "_")}} print}}' {input.fasta} > {params.wd}/input.fasta
        python {params.mitom}/mitom run -i {params.wd}/input.fasta -o {params.wd}/out --clade {params.clade} \
            --genetic_code {params.code} --threads {threads} --max_tasks 1 --id_maxlen 100 \
            --mitos_refdir {params.refdir} > {log} 2>&1
        cp {params.wd}/out/final_annotation_all.gb {output.gb}
        cp {params.wd}/out/summary.tsv {output.tsv}
        """

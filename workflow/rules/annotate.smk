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
        f"{OUT}/logs/mitom_db.log",
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
        f"{OUT}/logs/mitom_nuc_db.log",
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
        gb=f"{OUT}/annotation/{{sample}}/{{sample}}.mt.gb",
        tsv=f"{OUT}/annotation/{{sample}}/{{sample}}.mitom.tsv",
    params:
        wd=lambda w, output: str(Path(output.gb).parent / "mitom"),
        refdir=lambda w, input: config["mitos_refdir"] or f"{input.db}/mitos_refs",
        mitom=lambda w, input: config["mitom_dir"] or f"{input.db}/Mitom",
        clade=config["clade"],
        code=config["genetic_code"],
    log:
        f"{OUT}/logs/mitom/{{sample}}.log",
    benchmark:
        f"{OUT}/benchmarks/mitom/{{sample}}.tsv"
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


rule nr_annotate:
    """18S, ITS1, 5.8S, ITS2 and 28S in the rDNA with `mitom nuc annotate` (reference-based, boundaries projected
    from NCBI records; per Mitom's benchmark 92.6 % of boundaries within 5 nt)."""
    input:
        fasta=rules.getorganelle_nr.output.fasta,
        db=MITOM_NUC_DB,
    output:
        gb=f"{OUT}/rdna/{{sample}}/{{sample}}.nr.gb",
        tsv=f"{OUT}/rdna/{{sample}}/{{sample}}.nr_features.tsv",
    params:
        wd=lambda w, output: str(Path(output.gb).parent / "mitom_nuc"),
        mitom=lambda w, input: config["mitom_dir"] or f"{input.db}/Mitom",
    log:
        f"{OUT}/logs/nr_annotate/{{sample}}.log",
    threads: 4
    conda:
        "../envs/mitom.yaml"
    shell:
        """
        if [ ! -s {input.fasta} ]; then echo "no rDNA" > {log}; : > {output.gb}; : > {output.tsv}; exit 0; fi
        rm -rf {params.wd}; mkdir -p {params.wd}
        awk '{{if (substr($0, 1, 1) == ">") {{sub(" .*", ""); gsub("[^A-Za-z0-9_.>-]", "_")}} print}}' {input.fasta} > {params.wd}/input.fasta
        python {params.mitom}/mitom nuc annotate -i {params.wd}/input.fasta -o {params.wd}/out --genes rdna \
            --threads {threads} > {log} 2>&1
        cp {params.wd}/out/nuc_annotation.gb {output.gb}
        cp {params.wd}/out/nuc_features.tsv {output.tsv}
        """


rule nuc_markers:
    """Multi-copy nuclear markers in the whole-genome assembly with `mitom nuc`: rDNA (18S, ITS1, 5.8S, ITS2, 28S)
    and histones (H1, H2A, H2B, H3, H4); every copy plus one consensus per gene (contigs used / dropped and why)."""
    input:
        scaffolds=rules.spades.output.scaffolds,
        db=MITOM_NUC_DB,
    output:
        features=f"{OUT}/nuclear/{{sample}}/{{sample}}.nuc_features.tsv",
        consensus=f"{OUT}/nuclear/{{sample}}/{{sample}}.nuc_consensus.tsv",
    params:
        wd=lambda w, output: str(Path(output.features).parent / "mitom_nuc"),
        mitom=lambda w, input: config["mitom_dir"] or f"{input.db}/Mitom",
    log:
        f"{OUT}/logs/nuc_markers/{{sample}}.log",
    benchmark:
        f"{OUT}/benchmarks/nuc_markers/{{sample}}.tsv"
    threads: 16
    conda:
        "../envs/mitom.yaml"
    shell:
        """
        if [ ! -s {input.scaffolds} ]; then echo "no assembly" > {log}; : > {output.features}; : > {output.consensus}; exit 0; fi
        rm -rf {params.wd}
        python {params.mitom}/mitom nuc annotate -i {input.scaffolds} -o {params.wd} --genes rdna,histone \
            --threads {threads} > {log} 2>&1
        cp {params.wd}/nuc_features.tsv {output.features}
        cp {params.wd}/consensus.tsv {output.consensus}
        """


rule nuc_copies:
    """Copy number of each marker: k-mer coverage of its contigs / median coverage of the contigs carrying BUSCO
    complete single-copy genes (both from the SPAdes headers, so read length and k cancel out)."""
    input:
        features=rules.nuc_markers.output.features,
        consensus=rules.nuc_markers.output.consensus,
        busco=rules.busco.output.summary,
    output:
        tsv=f"{OUT}/nuclear/{{sample}}/{{sample}}.copy_number.tsv",
    params:
        busco_dir=lambda w, input: str(Path(input.busco).parent / "busco"),
    log:
        f"{OUT}/logs/nuc_copies/{{sample}}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/nuc_copies.py"

# Multi-copy nuclear markers: rDNA (18S-ITS1-5.8S-ITS2-28S) and histones (H1-H4). GetOrganelle assembles the
# tandem unit (rdna/, histone/), MitoM annotates it (.gbf) and the SPAdes assembly (nuclear/), nuc_compare compares
# the two. The same steps for an existing results directory: `wgs_skim nuclear` (scripts/nuclear_markers.py).


rule getorganelle_nr:
    """Nuclear rDNA (18S-ITS-28S) with GetOrganelle -F anonym and the bundled seed / label."""
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
        seed=NR_SEED,
        label=NR_LABEL,
    output:
        fasta=f"{OUT}rdna/{{sample}}/{{sample}}.nr.fasta",
    params:
        out=lambda w, output: str(Path(output.fasta).parent / "getorganelle"),
        collect=str(Path(workflow.basedir) / "scripts" / "collect_fasta.py"),
        k=config["getorganelle_kmers"],
        R=config["getorganelle_rounds"],
    log:
        f"{OUT}logs/getorganelle_nr/{{sample}}.log",
    benchmark:
        f"{OUT}benchmarks/getorganelle_nr/{{sample}}.tsv"
    threads: 16
    conda:
        "../envs/getorganelle.yaml"
    shell:
        """
        get_organelle_from_reads.py -1 {input.r1} -2 {input.r2} -o {params.out} -F anonym \
            -s {input.seed} --genes {input.label} -P 0 -R {params.R} -k {params.k} -t {threads} --overwrite > {log} 2>&1 \
            || echo "[wgs_skim] GetOrganelle rDNA failed (rc=$?)" >> {log}
        python {params.collect} {params.out} {wildcards.sample} nr > {output.fasta}
        """


rule nr_annotate:
    """18S, ITS1, 5.8S, ITS2 and 28S in the rDNA with `mitom nuc annotate` (reference-based, boundaries projected
    from NCBI records; per Mitom's benchmark 92.6 % of boundaries within 5 nt)."""
    input:
        fasta=rules.getorganelle_nr.output.fasta,
        db=MITOM_NUC_DB,
    output:
        gb=f"{OUT}rdna/{{sample}}/{{sample}}.rdna.gbf",
        tsv=f"{OUT}rdna/{{sample}}/{{sample}}.nr_features.tsv",
    params:
        wd=lambda w, output: str(Path(output.gb).parent / "mitom_nuc"),
        mitom=lambda w, input: config["mitom_dir"] or f"{input.db}/Mitom",
    log:
        f"{OUT}logs/nr_annotate/{{sample}}.log",
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
        features=f"{OUT}nuclear/{{sample}}/{{sample}}.nuc_features.tsv",
        consensus=f"{OUT}nuclear/{{sample}}/{{sample}}.nuc_consensus.tsv",
    params:
        wd=lambda w, output: str(Path(output.features).parent / "mitom_nuc"),
        mitom=lambda w, input: config["mitom_dir"] or f"{input.db}/Mitom",
    log:
        f"{OUT}logs/nuc_markers/{{sample}}.log",
    benchmark:
        f"{OUT}benchmarks/nuc_markers/{{sample}}.tsv"
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
        tsv=f"{OUT}nuclear/{{sample}}/{{sample}}.copy_number.tsv",
    params:
        busco_dir=lambda w, input: str(Path(input.busco).parent / "busco"),
    log:
        f"{OUT}logs/nuc_copies/{{sample}}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/nuc_copies.py"


rule getorganelle_histone:
    """Histone cluster (H1, H2A, H2B, H3, H4) with GetOrganelle -F anonym and the bundled H1-H4 seed / label
    (resources/histone; an H3-only label kept only the H3 piece). Same read settings as the rDNA."""
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
        seed=HIST_SEED,
        label=HIST_LABEL,
    output:
        fasta=f"{OUT}histone/{{sample}}/{{sample}}.histone.fasta",
    params:
        out=lambda w, output: str(Path(output.fasta).parent / "getorganelle"),
        nm=str(Path(workflow.basedir) / "scripts" / "nuclear_markers.py"),
        collect=str(Path(workflow.basedir) / "scripts" / "collect_fasta.py"),
        k=config["getorganelle_kmers"],
        R=config["getorganelle_rounds"],
    log:
        f"{OUT}logs/getorganelle_histone/{{sample}}.log",
    benchmark:
        f"{OUT}benchmarks/getorganelle_histone/{{sample}}.tsv"
    threads: 16
    conda:
        "../envs/getorganelle.yaml"
    shell:
        """
        rm -f {params.out}/done
        python {params.nm} histone -1 {input.r1} -2 {input.r2} -o {params.out} -t {threads} --kmers {params.k} \
            --rounds {params.R} --seed {input.seed} --label {input.label} > {log} 2>&1 \
            || echo "[wgs_skim] GetOrganelle histone failed (rc=$?)" >> {log}
        python {params.collect} {params.out} {wildcards.sample} histone > {output.fasta}
        """


rule histone_annotate:
    """H1-H4 in the histone cluster with `mitom nuc annotate`: gene models from 677 NCBI histone CDS, each copy
    classed replication-dependent (cluster) / replacement (H3.3, introns) / single-copy."""
    input:
        fasta=rules.getorganelle_histone.output.fasta,
        db=MITOM_NUC_DB,
    output:
        gb=f"{OUT}histone/{{sample}}/{{sample}}.histone.gbf",
    params:
        wd=lambda w, output: str(Path(output.gb).parent / "mitom_nuc"),
        nm=str(Path(workflow.basedir) / "scripts" / "nuclear_markers.py"),
        mitom=lambda w, input: config["mitom_dir"] or f"{input.db}/Mitom",
    log:
        f"{OUT}logs/histone_annotate/{{sample}}.log",
    threads: 4
    conda:
        "../envs/mitom.yaml"
    shell:
        """
        if [ ! -s {input.fasta} ]; then echo "no histone cluster" > {log}; : > {output.gb}; exit 0; fi
        python {params.nm} annotate -i {input.fasta} -o {params.wd} --kind histone --sample {wildcards.sample} \
            --mitom {params.mitom} -t {threads} > {log} 2>&1
        cp {params.wd}/{wildcards.sample}.histone.gbf {output.gb}
        """


rule nuc_compare:
    """GetOrganelle (rDNA unit, histone cluster) vs SPAdes-MitoM (copies over the whole assembly), per gene."""
    input:
        rdna=rules.nr_annotate.output.gb,
        histone=rules.histone_annotate.output.gb,
        spades=rules.nuc_markers.output.consensus,
    output:
        tsv=f"{OUT}nuclear/{{sample}}/{{sample}}.nuc_compare.tsv",
    params:
        nm=str(Path(workflow.basedir) / "scripts" / "nuclear_markers.py"),
        rdna=lambda w, input: str(Path(input.rdna).parent / "mitom_nuc" / "out"),
        histone=lambda w, input: str(Path(input.histone).parent / "mitom_nuc"),
        spades=lambda w, input: str(Path(input.spades).parent / "mitom_nuc"),
    log:
        f"{OUT}logs/nuc_compare/{{sample}}.log",
    conda:
        "../envs/mitom.yaml"
    shell:
        """
        python {params.nm} compare --sample {wildcards.sample} --getorganelle-rdna {params.rdna} \
            --getorganelle-histone {params.histone} --spades {params.spades} -o {output.tsv} > {log} 2>&1
        """

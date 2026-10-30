# Multi-copy nuclear markers: rDNA (18S-ITS1-5.8S-ITS2-28S) and histones (H1-H4). GetOrganelle assembles the
# tandem unit (rdna/, histone/), MitoM annotates it (.gbf) and the SPAdes assembly (nuclear/), nuc_compare compares
# the two. The same steps for an existing results directory: `wgs_skim nuclear` (scripts/nuclear_markers.py).


rule getorganelle_nr:
    """Nuclear rDNA (18S-ITS-28S) with GetOrganelle -F anonym and the bundled seed / label; rounds of at most
    nuclear_max_gb Gb of reads (5, then 10 seeded with round 1 when the unit does not close)."""
    input:
        unpack(clean_reads),
        seed=NR_SEED,
        label=NR_LABEL,
    output:
        fasta=f"{OUT}rdna/{{sample}}/{{sample}}.nr.fasta",
    params:
        out=lambda w, output: str(Path(output.fasta).parent / "getorganelle"),
        nm=str(Path(workflow.basedir) / "scripts" / "nuclear_markers.py"),
        collect=str(Path(workflow.basedir) / "scripts" / "collect_fasta.py"),
        k=config["getorganelle_kmers"],
        R=config["getorganelle_rounds"],
        gb=config["nuclear_max_gb"],
    resources:
        mem_mb=NUC_MEM,
    log:
        f"{OUT}pipeline_info/logs/getorganelle_nr/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/getorganelle_nr/{{sample}}.tsv"
    threads: 16
    conda:
        "../envs/getorganelle.yaml"
    shell:
        """
        rm -rf {params.out}
        python {params.nm} assemble --kind rdna -1 {input.r1} -2 {input.r2} -o {params.out} -t {threads} \
            --kmers {params.k} --rounds {params.R} --max-gb {params.gb} --seed {input.seed} --label {input.label} \
            > {log} 2>&1 || echo "[wgs_skim] GetOrganelle rDNA failed (rc=$?)" >> {log}
        python {params.collect} {params.out} {wildcards.sample} nr > {output.fasta}
        """


rule nr_annotate:
    """18S, ITS1, 5.8S, ITS2 and 28S in the rDNA with `mitom nuc annotate` (reference-based, boundaries projected
    from NCBI records; per Mitom's benchmark 92.6 % of boundaries within 5 nt)."""
    input:
        fasta=rules.getorganelle_nr.output.fasta,
        db=[] if config["mitom_dir"] else MITOM_NUC_DB,
    output:
        gb=f"{OUT}rdna/{{sample}}/{{sample}}.rdna.gbf",
        tsv=f"{OUT}rdna/{{sample}}/{{sample}}.nr_features.tsv",
    params:
        wd=lambda w, output: str(Path(output.gb).parent / "mitom_nuc"),
        mitom=lambda w, input: config["mitom_dir"] or f"{MITOM_NUC_DB}/Mitom",
    resources:
        mem_mb=2000,
    log:
        f"{OUT}pipeline_info/logs/nr_annotate/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/nr_annotate/{{sample}}.tsv"
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
        scaffolds=lambda w: scaffolds_of(w.sample),
        db=[] if config["mitom_dir"] else MITOM_NUC_DB,
    output:
        features=f"{OUT}nuclear_markers/{{sample}}/{{sample}}.nuc_features.tsv",
        consensus=f"{OUT}nuclear_markers/{{sample}}/{{sample}}.nuc_consensus.tsv",
    params:
        wd=lambda w, output: str(Path(output.features).parent / "mitom_nuc"),
        mitom=lambda w, input: config["mitom_dir"] or f"{MITOM_NUC_DB}/Mitom",
    resources:
        mem_mb=8000,
    log:
        f"{OUT}pipeline_info/logs/nuc_markers/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/nuc_markers/{{sample}}.tsv"
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
        tsv=f"{OUT}nuclear_markers/{{sample}}/{{sample}}.copy_number.tsv",
    params:
        busco_dir=lambda w, input: str(Path(input.busco).parent / "busco"),
    resources:
        mem_mb=1000,
    log:
        f"{OUT}pipeline_info/logs/nuc_copies/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/nuc_copies/{{sample}}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/nuc_copies.py"


rule getorganelle_histone:
    """Histone cluster (H1, H2A, H2B, H3, H4) with GetOrganelle -F anonym and the bundled H1-H4 seed / label
    (resources/histone; an H3-only label kept only the H3 piece). Rounds of reads as for the rDNA."""
    input:
        unpack(clean_reads),
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
        gb=config["nuclear_max_gb"],
    resources:
        mem_mb=NUC_MEM,
    log:
        f"{OUT}pipeline_info/logs/getorganelle_histone/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/getorganelle_histone/{{sample}}.tsv"
    threads: 16
    conda:
        "../envs/getorganelle.yaml"
    shell:
        """
        rm -rf {params.out}
        python {params.nm} assemble --kind histone -1 {input.r1} -2 {input.r2} -o {params.out} -t {threads} \
            --kmers {params.k} --rounds {params.R} --max-gb {params.gb} --seed {input.seed} --label {input.label} \
            > {log} 2>&1 || echo "[wgs_skim] GetOrganelle histone failed (rc=$?)" >> {log}
        python {params.collect} {params.out} {wildcards.sample} histone > {output.fasta}
        """


rule histone_annotate:
    """H1-H4 in the histone cluster with `mitom nuc annotate`: gene models from 677 NCBI histone CDS, each copy
    classed replication-dependent (cluster) / replacement (H3.3, introns) / single-copy."""
    input:
        fasta=rules.getorganelle_histone.output.fasta,
        db=[] if config["mitom_dir"] else MITOM_NUC_DB,
    output:
        gb=f"{OUT}histone/{{sample}}/{{sample}}.histone.gbf",
    params:
        wd=lambda w, output: str(Path(output.gb).parent / "mitom_nuc"),
        nm=str(Path(workflow.basedir) / "scripts" / "nuclear_markers.py"),
        mitom=lambda w, input: config["mitom_dir"] or f"{MITOM_NUC_DB}/Mitom",
    resources:
        mem_mb=2000,
    log:
        f"{OUT}pipeline_info/logs/histone_annotate/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/histone_annotate/{{sample}}.tsv"
    threads: 4
    conda:
        "../envs/mitom.yaml"
    shell:
        """
        if [ ! -s {input.fasta} ]; then echo "no histone cluster" > {log}; : > {output.gb}; exit 0; fi
        rm -rf {params.wd}
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
        tsv=f"{OUT}nuclear_markers/{{sample}}/{{sample}}.nuc_compare.tsv",
    params:
        nm=str(Path(workflow.basedir) / "scripts" / "nuclear_markers.py"),
        rdna=lambda w, input: str(Path(input.rdna).parent / "mitom_nuc" / "out"),
        histone=lambda w, input: str(Path(input.histone).parent / "mitom_nuc"),
        spades=lambda w, input: str(Path(input.spades).parent / "mitom_nuc"),
    resources:
        mem_mb=1000,
    log:
        f"{OUT}pipeline_info/logs/nuc_compare/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/nuc_compare/{{sample}}.tsv"
    conda:
        "../envs/mitom.yaml"
    shell:
        """
        python {params.nm} compare --sample {wildcards.sample} --getorganelle-rdna {params.rdna} \
            --getorganelle-histone {params.histone} --spades {params.spades} -o {output.tsv} > {log} 2>&1
        """

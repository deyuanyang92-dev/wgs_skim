GETORG_DB = str(DB / "getorganelle_db")


rule getorganelle_db:
    """GetOrganelle label/seed databases (animal_mt for the mitogenome; rDNA uses the bundled seed/label)."""
    output:
        directory(GETORG_DB),
    log:
        f"{OUT}/logs/getorganelle_db.log",
    conda:
        "../envs/getorganelle.yaml"
    shell:
        "get_organelle_config.py -a animal_mt --config-dir {output} > {log} 2>&1"


rule meangs:
    """MEANGS finds a mitochondrial seed; an empty seed file means GetOrganelle falls back to its animal_mt seeds."""
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
    output:
        seed=f"{OUT}/mito/{{sample}}/meangs/{{sample}}.seed.fasta",
    params:
        wd=lambda w, output: str(Path(output.seed).parent),
        cls=config["meangs_class"],
        n=config["meangs_reads"],
        insert_size=config["insert_size"],
    log:
        f"{OUT}/logs/meangs/{{sample}}.log",
    benchmark:
        f"{OUT}/benchmarks/meangs/{{sample}}.tsv"
    threads: 8
    conda:
        "../envs/meangs.yaml"
    shell:
        """
        r1=$(readlink -f {input.r1}); r2=$(readlink -f {input.r2}); log=$(readlink -f {log})
        cd {params.wd}
        meangs.py --silence -1 $r1 -2 $r2 -o meangs -t {threads} -i {params.insert_size} \
            --species_class {params.cls} -n {params.n} --deepin > $log 2>&1 \
            || echo "[wgs_skim] MEANGS failed (rc=$?): GetOrganelle will use its own animal_mt seeds" >> $log
        seed=$(find . -name 'scaffold_seeds.fa*' | head -n 1)
        if [ -n "$seed" ] && [ -s "$seed" ]; then cp "$seed" {wildcards.sample}.seed.fasta; else : > {wildcards.sample}.seed.fasta; fi
        rm -rf meangs
        """


rule getorganelle_mt:
    """First GetOrganelle mitogenome assembly from the MEANGS seed (GetOrganelle's default read subsample)."""
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
        seed=mt_seed,
        db=GETORG_DB,
    output:
        done=f"{OUT}/mito/{{sample}}/round0/.done",
    params:
        out=lambda w, output: str(Path(output.done).parent),
        k=config["getorganelle_kmers"],
        R=config["getorganelle_rounds"],
    log:
        f"{OUT}/logs/getorganelle_mt/{{sample}}.log",
    benchmark:
        f"{OUT}/benchmarks/getorganelle_mt/{{sample}}.tsv"
    threads: 16
    conda:
        "../envs/getorganelle.yaml"
    shell:
        """
        seed=""; [ -s {input.seed} ] && seed="-s {input.seed}"
        get_organelle_from_reads.py -1 {input.r1} -2 {input.r2} -o {params.out} -F animal_mt \
            -R {params.R} -k {params.k} -t {threads} --config-dir {input.db} $seed --overwrite > {log} 2>&1 \
            || echo "[wgs_skim] GetOrganelle round 0 failed (rc=$?): mt_iterate reseeds it" >> {log}
        touch {output.done}
        """


rule mt_iterate:
    """Keep round 0 if circular; otherwise reseed GetOrganelle with the previous round, all reads (capped)."""
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
        json=rules.fastp.output.json,
        round0=rules.getorganelle_mt.output.done,
        db=GETORG_DB,
    output:
        fasta=f"{OUT}/mito/{{sample}}/{{sample}}.getorganelle.fasta",
        status=f"{OUT}/mito/{{sample}}/{{sample}}.getorganelle.json",
    params:
        rounds=config["mt_iter_rounds"],
        R=config["mt_iter_R"],
        k=config["getorganelle_kmers"],
        max_gb=config["mt_iter_max_data_gb"],
    log:
        f"{OUT}/logs/mt_iterate/{{sample}}.log",
    benchmark:
        f"{OUT}/benchmarks/mt_iterate/{{sample}}.tsv"
    threads: 16
    conda:
        "../envs/getorganelle.yaml"
    script:
        "../scripts/mt_iterate.py"


rule mt_select:
    """Final mitogenome: GetOrganelle, or the SPAdes scaffold when it closes the circle / is longer (see script)."""
    input:
        unpack(mt_select_inputs),
    output:
        fasta=f"{OUT}/mito/{{sample}}/{{sample}}.mt.fasta",
        status=f"{OUT}/mito/{{sample}}/{{sample}}.mt.json",
    log:
        f"{OUT}/logs/mt_select/{{sample}}.log",
    threads: 8
    conda:
        "../envs/blast.yaml"
    script:
        "../scripts/mt_select.py"


rule getorganelle_nr:
    """Nuclear rDNA (18S-ITS-28S) with GetOrganelle -F anonym and the bundled seed / label."""
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
        seed=NR_SEED,
        label=NR_LABEL,
    output:
        fasta=f"{OUT}/rdna/{{sample}}/{{sample}}.nr.fasta",
    params:
        out=lambda w, output: str(Path(output.fasta).parent / "getorganelle"),
        collect=str(Path(workflow.basedir) / "scripts" / "collect_fasta.py"),
        k=config["getorganelle_kmers"],
        R=config["getorganelle_rounds"],
    log:
        f"{OUT}/logs/getorganelle_nr/{{sample}}.log",
    benchmark:
        f"{OUT}/benchmarks/getorganelle_nr/{{sample}}.tsv"
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

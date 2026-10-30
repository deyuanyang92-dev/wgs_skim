// Mitogenome: MEANGS seed -> GetOrganelle -> reseeded rounds -> final choice (workflow/rules/mito.smk)
include { cfg; condaEnv; getorganelleMem; smkJson; scriptsDir } from './common'

process GETORGANELLE_DB {
    // made once under db_dir and reused by every run (Snakemake: directory(GETORG_DB))
    storeDir "${cfg().db}"
    memory '1000 MB'
    conda { condaEnv('getorganelle') }

    output:
    path 'getorganelle_db'

    script:
    """
    get_organelle_config.py -a animal_mt --config-dir getorganelle_db > getorganelle_db.log 2>&1
    """

    stub:
    """
    mkdir -p getorganelle_db
    """
}

process MEANGS {
    tag "$s"
    cpus 8
    memory '12000 MB'
    conda { condaEnv('meangs') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(r1), path(r2)

    output:
    tuple val(s), path("meangs/${s}/${s}.seed.fasta"), emit: seed
    path "meangs/${s}/scaffold_seeds.fa*", optional: true
    path "pipeline_info/logs/meangs/${s}.log"

    script:
    def c = cfg()
    """
    mkdir -p meangs/${s} pipeline_info/logs/meangs
    r1=\$(readlink -f ${r1}); r2=\$(readlink -f ${r2}); log=\$PWD/pipeline_info/logs/meangs/${s}.log
    cd meangs/${s}
    meangs.py --silence -1 \$r1 -2 \$r2 -o meangs -t ${task.cpus} -i ${c.insert_size} \\
        --species_class ${c.meangs_class} -n ${c.meangs_reads} --deepin > \$log 2>&1 \\
        || echo "[wgs_skim] MEANGS failed (rc=\$?): GetOrganelle will use its own animal_mt seeds" >> \$log
    seed=\$(find . -name 'scaffold_seeds.fa*' | head -n 1)
    if [ -n "\$seed" ] && [ -s "\$seed" ]; then cp "\$seed" ${s}.seed.fasta; else : > ${s}.seed.fasta; fi
    rm -rf meangs
    """

    stub:
    """
    mkdir -p meangs/${s} pipeline_info/logs/meangs
    echo -e '>seed\\nACGT' > meangs/${s}/${s}.seed.fasta; touch pipeline_info/logs/meangs/${s}.log
    """
}

process GETORGANELLE_MT {
    tag "$s"
    cpus 16
    memory { "${16000 * task.attempt} MB" }
    conda { condaEnv('getorganelle') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(r1), path(r2), path(seed, stageAs: 'seed/*')
    path db

    output:
    tuple val(s), path("getorganelle/${s}/round0"), emit: round0
    path "pipeline_info/logs/getorganelle_mt/${s}.log"

    script:
    def c = cfg()
    """
    mkdir -p getorganelle/${s} pipeline_info/logs/getorganelle_mt
    log=pipeline_info/logs/getorganelle_mt/${s}.log
    sd=""; [ -s ${seed} ] && sd="-s ${seed}"
    get_organelle_from_reads.py -1 ${r1} -2 ${r2} -o getorganelle/${s}/round0 -F animal_mt \\
        -R ${c.getorganelle_rounds} -k ${c.getorganelle_kmers} -t ${task.cpus} --config-dir ${db} \$sd --overwrite > \$log 2>&1 \\
        || echo "[wgs_skim] GetOrganelle round 0 failed (rc=\$?): mt_iterate reseeds it" >> \$log
    mkdir -p getorganelle/${s}/round0; touch getorganelle/${s}/round0/.done
    # mt_iterate removes these intermediates anyway; doing it here keeps this task's output unchanged afterwards,
    # which -resume needs (a modified input directory would rerun mt_iterate)
    python ${scriptsDir()}/collect_fasta.py --clean getorganelle/${s}/round0 2>> \$log
    """

    stub:
    """
    mkdir -p getorganelle/${s}/round0 pipeline_info/logs/getorganelle_mt
    touch getorganelle/${s}/round0/.done pipeline_info/logs/getorganelle_mt/${s}.log
    """
}

process MT_ITERATE {
    tag "$s"
    cpus 16
    memory { getorganelleMem(gb, cfg().mt_iter_max_data_gb, task.attempt) }
    conda { condaEnv('getorganelle') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(r1), path(r2), path(json), val(gb), path(round0)
    path db

    output:
    tuple val(s), path("getorganelle/${s}/${s}.getorganelle.fasta"), path("getorganelle/${s}/${s}.getorganelle.json"), emit: result
    path "getorganelle/${s}/round[1-9]*", optional: true
    path "pipeline_info/logs/mt_iterate/${s}.log"

    script:
    def c = cfg()
    def spec = smkJson(
        input: [r1: "${r1}", r2: "${r2}", json: "${json}", round0: "getorganelle/${s}/round0/.done", db: "${db}"],
        output: [fasta: "getorganelle/${s}/${s}.getorganelle.fasta", status: "getorganelle/${s}/${s}.getorganelle.json"],
        params: [rounds: c.mt_iter_rounds, R: c.mt_iter_R, k: c.getorganelle_kmers, max_gb: c.mt_iter_max_data_gb],
        log: ["pipeline_info/logs/mt_iterate/${s}.log"], threads: task.cpus, wildcards: [sample: s], rule: 'mt_iterate')
    """
    mkdir -p getorganelle/${s} pipeline_info/logs/mt_iterate
    ln -s ../../${round0} getorganelle/${s}/round0     # the script reads round0 next to its outputs
    ${spec}
    smk_shim.py snakemake.json ${scriptsDir()}/mt_iterate.py
    """

    stub:
    """
    mkdir -p getorganelle/${s} pipeline_info/logs/mt_iterate
    echo -e '>${s}|mt|1 linear length=4\\nACGT' > getorganelle/${s}/${s}.getorganelle.fasta
    echo '{"sample": "${s}", "circular": false, "round": "round0", "length": 4}' > getorganelle/${s}/${s}.getorganelle.json
    touch pipeline_info/logs/mt_iterate/${s}.log
    """
}

process MT_SELECT {
    tag "$s"
    cpus 8
    memory '2000 MB'
    conda { condaEnv('blast') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    // asm: [scaffolds, spades directory] (mt_from_spades) or []: staged side by side, as in spades/<sample>/
    tuple val(s), path(fasta), path(status), path(seed, stageAs: 'seed/*'), path(asm, stageAs: 'asm/*')

    output:
    tuple val(s), path("mitogenome/${s}/${s}.mt.fasta"), path("mitogenome/${s}/${s}.mt.json"), emit: result
    path "pipeline_info/logs/mt_select/${s}.log"

    script:
    def inp = [fasta: "${fasta}", status: "${status}", seed: "${seed}"]
    if (asm) inp.scaffolds = "asm/${s}.scaffolds.fasta".toString()
    def spec = smkJson(input: inp, output: [fasta: "mitogenome/${s}/${s}.mt.fasta", status: "mitogenome/${s}/${s}.mt.json"],
        log: ["pipeline_info/logs/mt_select/${s}.log"], threads: task.cpus, wildcards: [sample: s], rule: 'mt_select')
    """
    mkdir -p mitogenome/${s} pipeline_info/logs/mt_select
    ${spec}
    smk_shim.py snakemake.json ${scriptsDir()}/mt_select.py
    """

    stub:
    """
    mkdir -p mitogenome/${s} pipeline_info/logs/mt_select
    cp ${fasta} mitogenome/${s}/${s}.mt.fasta
    echo '{"sample": "${s}", "circular": false, "round": "round0", "length": 4, "source": "getorganelle", "note": ""}' > mitogenome/${s}/${s}.mt.json
    touch pipeline_info/logs/mt_select/${s}.log
    """
}

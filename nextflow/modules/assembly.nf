// Whole-genome assembly + BUSCO; genome size (workflow/rules/assembly.smk, genome_size.smk)
include { cfg; condaEnv; spadesMem } from './common'

process SPADES {
    tag "$s"
    cpus 32
    memory { spadesMem(gb, given, task.attempt) }
    // -m is a hard limit: one retry with 1.5x memory (Snakemake: retries: 1); a retry starts in a new work dir,
    // so SPAdes starts over instead of --restart-from last
    errorStrategy { task.attempt <= 1 ? 'retry' : 'finish' }
    maxRetries 1
    conda { condaEnv('spades') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(r1), path(r2), path(json), val(gb), val(given)

    output:
    // the spades directory goes with the scaffolds: mt_select walks its assembly graph (scaffolds.paths)
    tuple val(s), path("assembly/${s}/${s}.scaffolds.fasta"), path("assembly/${s}/spades"), emit: asm
    path "logs/spades/${s}.log"

    script:
    def out = "assembly/${s}/spades"
    def log = "logs/spades/${s}.log"
    def min_gb = cfg().spades_min_clean_gb
    """
    mkdir -p assembly/${s} logs/spades
    if [ -n "${given}" ]; then
        echo "[wgs_skim] assembly given in the sample sheet: ${given}" > ${log}
        mkdir -p ${out}; ln -sf "${given}" assembly/${s}/${s}.scaffolds.fasta; exit 0
    fi
    gb=\$(python -c "import json;print(json.load(open('${json}'))['summary']['after_filtering']['total_bases']/1e9)")
    if python -c "import sys;sys.exit(not (${min_gb} > 0 and \$gb < ${min_gb}))"; then
        echo "[wgs_skim] clean data \$gb Gb < spades_min_clean_gb ${min_gb}: not assembled" > ${log}
        mkdir -p ${out}; : > assembly/${s}/${s}.scaffolds.fasta; exit 0
    fi
    mem=${task.memory.toMega().intdiv(1024)}
    if [ -f ${out}/params.txt ]; then
        spades.py -o ${out} --restart-from last -t ${task.cpus} -m \$mem > ${log} 2>&1
    else
        spades.py -1 ${r1} -2 ${r2} -o ${out} -t ${task.cpus} -m \$mem --only-assembler > ${log} 2>&1
    fi
    cp ${out}/scaffolds.fasta assembly/${s}/${s}.scaffolds.fasta
    rm -rf ${out}/K* ${out}/tmp ${out}/misc ${out}/before_rr.fasta   # intermediates (~2 GB/Gb)
    """

    stub:
    """
    mkdir -p assembly/${s}/spades logs/spades
    echo -e '>NODE_1_length_4_cov_10.0\\nACGT' > assembly/${s}/${s}.scaffolds.fasta; touch logs/spades/${s}.log
    """
}

process BUSCO_LINEAGE {
    // downloaded once into db_dir/busco_downloads/lineages (skipped when busco_lineage is a local directory)
    tag "$lineage"
    storeDir "${cfg().db}/busco_downloads/lineages"
    memory '1000 MB'
    conda { condaEnv('busco') }

    input:
    val lineage

    output:
    path "${lineage}"

    script:
    """
    busco --download ${lineage} --download_path dl > busco_download.log 2>&1
    mv dl/lineages/${lineage} .
    """

    stub:
    """
    mkdir -p ${lineage}
    """
}

process BUSCO {
    tag "$s"
    cpus 16
    memory { "${16000 * task.attempt} MB" }
    conda { condaEnv('busco') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(scaffolds)
    path lineage

    output:
    // the busco directory: nuc_copies reads the coverage of single-copy genes from its full_table.tsv
    tuple val(s), path("assembly/${s}/${s}.busco.txt"), path("assembly/${s}/busco"), emit: summary
    path "logs/busco/${s}.log", optional: true

    script:
    def out = "assembly/${s}"
    """
    mkdir -p ${out} logs/busco
    if [ ! -s ${scaffolds} ]; then echo "no assembly" > ${out}/${s}.busco.txt; mkdir -p ${out}/busco; exit 0; fi
    rm -rf ${out}/busco
    # a full path: BUSCO takes a bare name (the staged link) for a lineage to look up in ./busco_downloads
    busco -i ${scaffolds} -l \$(readlink -f ${lineage}) -m genome -o busco --out_path ${out} \\
        -c ${task.cpus} --offline -f > logs/busco/${s}.log 2>&1
    cp ${out}/busco/short_summary.*.txt ${out}/${s}.busco.txt
    """

    stub:
    """
    mkdir -p assembly/${s}/busco logs/busco
    echo 'C:1.0%[S:1.0%,D:0.0%],F:0.0%,M:99.0%,n:954' > assembly/${s}/${s}.busco.txt; touch logs/busco/${s}.log
    """
}

process JELLYFISH {
    tag "$s"
    cpus 16
    memory { "${12000 * task.attempt} MB" }
    conda { condaEnv('genome_size') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(r1), path(r2)

    output:
    tuple val(s), path("genome_size/${s}/${s}.k${params.kmer_size}.histo"), emit: histo
    path "logs/jellyfish/${s}.log"

    script:
    def k = params.kmer_size
    """
    mkdir -p genome_size/${s} logs/jellyfish
    jellyfish count -C -m ${k} -s 2G -t ${task.cpus} -o genome_size/${s}/reads.jf \\
        <(zcat ${r1}) <(zcat ${r2}) > logs/jellyfish/${s}.log 2>&1
    jellyfish histo -t ${task.cpus} genome_size/${s}/reads.jf > genome_size/${s}/${s}.k${k}.histo 2>> logs/jellyfish/${s}.log
    rm -f genome_size/${s}/reads.jf
    """

    stub:
    """
    mkdir -p genome_size/${s} logs/jellyfish
    echo '1 100' > genome_size/${s}/${s}.k${params.kmer_size}.histo; touch logs/jellyfish/${s}.log
    """
}

process GENOMESCOPE {
    tag "$s"
    memory '2000 MB'
    conda { condaEnv('genome_size') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(histo)

    output:
    tuple val(s), path("genome_size/${s}/${s}.genomescope.txt"), emit: summary
    path "genome_size/${s}/genomescope", optional: true
    path "logs/genomescope/${s}.log"

    script:
    def out = "genome_size/${s}/genomescope"
    def log = "logs/genomescope/${s}.log"
    """
    mkdir -p genome_size/${s} logs/genomescope
    genomescope2 -i ${histo} -o ${out} -k ${params.kmer_size} -p ${params.ploidy} -n ${s} > ${log} 2>&1 \\
        || echo "[wgs_skim] GenomeScope did not converge (low coverage?)" >> ${log}
    if [ -s ${out}/${s}_summary.txt ]; then cp ${out}/${s}_summary.txt genome_size/${s}/${s}.genomescope.txt
    else echo "GenomeScope failed, see ${log}" > genome_size/${s}/${s}.genomescope.txt; fi
    """

    stub:
    """
    mkdir -p genome_size/${s} logs/genomescope
    echo 'GenomeScope stub' > genome_size/${s}/${s}.genomescope.txt; touch logs/genomescope/${s}.log
    """
}

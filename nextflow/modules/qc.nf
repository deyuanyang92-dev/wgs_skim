include { condaEnv } from './common'

process FASTP {
    tag "$s"
    cpus 8
    memory '2000 MB'
    conda { condaEnv('fastp') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(r1, stageAs: 'raw/*'), path(r2, stageAs: 'raw/*')

    output:
    tuple val(s), path("qc/${s}/${s}.R1.fq.gz"), path("qc/${s}/${s}.R2.fq.gz"), path("qc/${s}/${s}.fastp.json"), emit: reads
    path "qc/${s}/${s}.fastp.html"
    path "logs/fastp/${s}.log"

    script:
    """
    mkdir -p qc/${s} logs/fastp
    fastp -i ${r1} -I ${r2} -o qc/${s}/${s}.R1.fq.gz -O qc/${s}/${s}.R2.fq.gz --detect_adapter_for_pe \\
        -w ${task.cpus} -h qc/${s}/${s}.fastp.html -j qc/${s}/${s}.fastp.json > logs/fastp/${s}.log 2>&1
    """

    stub:
    """
    mkdir -p qc/${s} logs/fastp
    echo | gzip > qc/${s}/${s}.R1.fq.gz; echo | gzip > qc/${s}/${s}.R2.fq.gz
    echo '{"summary": {"after_filtering": {"total_bases": 300000000, "q30_rate": 0.97, "read1_mean_length": 150}}}' > qc/${s}/${s}.fastp.json
    touch qc/${s}/${s}.fastp.html logs/fastp/${s}.log
    """
}

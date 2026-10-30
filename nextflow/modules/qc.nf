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
    tuple val(s), path("fastp/${s}/${s}.R1.fq.gz"), path("fastp/${s}/${s}.R2.fq.gz"), path("fastp/${s}/${s}.fastp.json"), emit: reads
    path "fastp/${s}/${s}.fastp.html"
    path "pipeline_info/logs/fastp/${s}.log"

    script:
    """
    mkdir -p fastp/${s} pipeline_info/logs/fastp
    fastp -i ${r1} -I ${r2} -o fastp/${s}/${s}.R1.fq.gz -O fastp/${s}/${s}.R2.fq.gz --detect_adapter_for_pe \\
        -w ${task.cpus} -h fastp/${s}/${s}.fastp.html -j fastp/${s}/${s}.fastp.json > pipeline_info/logs/fastp/${s}.log 2>&1
    """

    stub:
    """
    mkdir -p fastp/${s} pipeline_info/logs/fastp
    echo | gzip > fastp/${s}/${s}.R1.fq.gz; echo | gzip > fastp/${s}/${s}.R2.fq.gz
    echo '{"summary": {"after_filtering": {"total_bases": 300000000, "q30_rate": 0.97, "read1_mean_length": 150}}}' > fastp/${s}/${s}.fastp.json
    touch fastp/${s}/${s}.fastp.html pipeline_info/logs/fastp/${s}.log
    """
}

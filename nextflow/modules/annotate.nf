// MitoM annotation of the final mitogenome (workflow/rules/annotate.smk)
include { cfg; condaEnv } from './common'

process MITOM_DB {
    // MitoM code (pinned commit) and the MITOS2 reference data refseq89m (skipped if mitos_refdir is set);
    // made once under db_dir
    storeDir "${cfg().db}"
    memory '1000 MB'
    conda { condaEnv('mitom') }

    output:
    path 'mitom'

    script:
    def c = cfg()
    """
    git clone -q https://github.com/deyuanyang92-dev/Mitom.git mitom/Mitom > mitom_db.log 2>&1
    git -C mitom/Mitom checkout -q ${c.mitom_commit} >> mitom_db.log 2>&1
    if [ -z "${c.mitos_refdir}" ]; then
        mkdir -p mitom/mitos_refs
        wget -q -O mitom/mitos_refs/refseq89m.tar.bz2 https://zenodo.org/records/4284483/files/refseq89m.tar.bz2 >> mitom_db.log 2>&1
        echo "bb6325e27612e61a6995e88e8ffcecf8  mitom/mitos_refs/refseq89m.tar.bz2" | md5sum -c >> mitom_db.log 2>&1
        tar -xjf mitom/mitos_refs/refseq89m.tar.bz2 -C mitom/mitos_refs && rm mitom/mitos_refs/refseq89m.tar.bz2
    fi
    """

    stub:
    """
    mkdir -p mitom/Mitom mitom/mitos_refs
    """
}

process MITOM_NUC_DB {
    // MitoM version with `mitom nuc` (rDNA + histone annotation; reference library included)
    storeDir "${cfg().db}"
    memory '1000 MB'
    conda { condaEnv('mitom') }

    output:
    path 'mitom_nuc'

    script:
    """
    git clone -q https://github.com/deyuanyang92-dev/Mitom.git mitom_nuc/Mitom > mitom_nuc_db.log 2>&1
    git -C mitom_nuc/Mitom checkout -q ${cfg().mitom_nuc_commit} >> mitom_nuc_db.log 2>&1
    """

    stub:
    """
    mkdir -p mitom_nuc/Mitom
    """
}

process MITOM_ANNOTATE {
    tag "$s"
    cpus 8
    memory '4000 MB'
    conda { condaEnv('mitom') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(fasta)
    path db, stageAs: "staged_db"             // MITOM_DB output, or [] when both mitom_dir and mitos_refdir are given

    output:
    tuple val(s), path("mitom/${s}/${s}.mt.gb"), emit: gb
    path "mitom/${s}/${s}.mitom.tsv"
    path "mitom/${s}/mitom", optional: true
    path "pipeline_info/logs/mitom/${s}.log"

    script:
    def c = cfg()
    def wd = "mitom/${s}/mitom"
    def refdir = c.mitos_refdir ?: "\$PWD/${db}/mitos_refs"
    def mitom = c.mitom_dir ?: "\$PWD/${db}/Mitom"
    def log = "pipeline_info/logs/mitom/${s}.log"
    """
    mkdir -p mitom/${s} pipeline_info/logs/mitom
    if [ ! -s ${fasta} ]; then echo "no mitogenome" > ${log}; : > mitom/${s}/${s}.mt.gb; : > mitom/${s}/${s}.mitom.tsv; exit 0; fi
    rm -rf ${wd}; mkdir -p ${wd}
    # record IDs "S|mt|1" -> "S_mt_1" (topology matches the annotation to the sequence by content)
    awk '{if (substr(\$0, 1, 1) == ">") {sub(" .*", ""); gsub("[^A-Za-z0-9_.>-]", "_")} print}' ${fasta} > ${wd}/input.fasta
    python ${mitom}/mitom run -i ${wd}/input.fasta -o ${wd}/out --clade ${c.clade} \\
        --genetic_code ${c.genetic_code} --threads ${task.cpus} --max_tasks 1 --id_maxlen 100 \\
        --mitos_refdir ${refdir} > ${log} 2>&1
    cp ${wd}/out/final_annotation_all.gb mitom/${s}/${s}.mt.gb
    cp ${wd}/out/summary.tsv mitom/${s}/${s}.mitom.tsv
    """

    stub:
    """
    mkdir -p mitom/${s}/mitom pipeline_info/logs/mitom
    touch mitom/${s}/${s}.mt.gb mitom/${s}/${s}.mitom.tsv pipeline_info/logs/mitom/${s}.log
    """
}

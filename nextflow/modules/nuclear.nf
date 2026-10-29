// Multi-copy nuclear markers: rDNA (18S-ITS1-5.8S-ITS2-28S) and histones (H1-H4) (workflow/rules/nuclear.smk).
// GetOrganelle assembles the tandem unit (rdna/, histone/), MitoM annotates it (.gbf) and the SPAdes assembly
// (nuclear/), NUC_COMPARE compares the two.
include { cfg; condaEnv; getorganelleMem; nucMaxGb; smkJson; scriptsDir } from './common'

process GETORGANELLE_NUC {
    // one definition for both markers (included as GETORGANELLE_NR and GETORGANELLE_HISTONE): GetOrganelle -F anonym
    // with the bundled seed / label, rounds of at most nuclear_max_gb Gb of reads (5, then 10 seeded with round 1)
    tag "$s"
    cpus 16
    memory { getorganelleMem(gb, nucMaxGb(), task.attempt) }
    conda { condaEnv('getorganelle') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(r1), path(r2), path(json), val(gb)
    tuple val(kind), path(seed, stageAs: 'ref/seed.fasta'), path(label, stageAs: 'ref/label.fasta')

    output:
    tuple val(s), path("${nucDir(kind)}/${s}/${s}.${nucLocus(kind)}.fasta"), emit: fasta
    path "${nucDir(kind)}/${s}/getorganelle", optional: true
    path "logs/getorganelle_${kind == 'rdna' ? 'nr' : 'histone'}/${s}.log"

    script:
    def c = cfg()
    def out = "${nucDir(kind)}/${s}/getorganelle"
    def log = "logs/getorganelle_${kind == 'rdna' ? 'nr' : 'histone'}/${s}.log"
    """
    mkdir -p ${nucDir(kind)}/${s} \$(dirname ${log})
    rm -rf ${out}
    python ${scriptsDir()}/nuclear_markers.py assemble --kind ${kind} -1 ${r1} -2 ${r2} -o ${out} -t ${task.cpus} \\
        --kmers ${c.getorganelle_kmers} --rounds ${c.getorganelle_rounds} --max-gb ${c.nuclear_max_gb} \\
        --seed ${seed} --label ${label} > ${log} 2>&1 || echo "[wgs_skim] GetOrganelle ${kind} failed (rc=\$?)" >> ${log}
    python ${scriptsDir()}/collect_fasta.py ${out} ${s} ${nucLocus(kind)} > ${nucDir(kind)}/${s}/${s}.${nucLocus(kind)}.fasta
    """

    stub:
    def log = "logs/getorganelle_${kind == 'rdna' ? 'nr' : 'histone'}/${s}.log"
    """
    mkdir -p ${nucDir(kind)}/${s}/getorganelle \$(dirname ${log}); touch ${log}
    echo -e '>${s}|${nucLocus(kind)}|1 linear length=4\\nACGT' > ${nucDir(kind)}/${s}/${s}.${nucLocus(kind)}.fasta
    """
}

def nucDir(kind) { kind == 'rdna' ? 'rdna' : 'histone' }

def nucLocus(kind) { kind == 'rdna' ? 'nr' : 'histone' }

process NR_ANNOTATE {
    // 18S, ITS1, 5.8S, ITS2 and 28S with `mitom nuc annotate` (reference-based)
    tag "$s"
    cpus 4
    memory '2000 MB'
    conda { condaEnv('mitom') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(fasta)
    path db             // MITOM_NUC_DB output, or [] with mitom_dir

    output:
    // mitom_nuc/out is read by NUC_COMPARE; created empty without rDNA (a missing directory means the same there)
    tuple val(s), path("rdna/${s}/${s}.rdna.gbf"), path("rdna/${s}/mitom_nuc"), emit: gb
    path "rdna/${s}/${s}.nr_features.tsv"
    path "logs/nr_annotate/${s}.log"

    script:
    def wd = "rdna/${s}/mitom_nuc"
    def mitom = cfg().mitom_dir ?: "\$PWD/${db}/Mitom"
    def log = "logs/nr_annotate/${s}.log"
    """
    mkdir -p rdna/${s} logs/nr_annotate
    if [ ! -s ${fasta} ]; then echo "no rDNA" > ${log}; : > rdna/${s}/${s}.rdna.gbf; : > rdna/${s}/${s}.nr_features.tsv; mkdir -p ${wd}; exit 0; fi
    rm -rf ${wd}; mkdir -p ${wd}
    awk '{if (substr(\$0, 1, 1) == ">") {sub(" .*", ""); gsub("[^A-Za-z0-9_.>-]", "_")} print}' ${fasta} > ${wd}/input.fasta
    python ${mitom}/mitom nuc annotate -i ${wd}/input.fasta -o ${wd}/out --genes rdna \\
        --threads ${task.cpus} > ${log} 2>&1
    cp ${wd}/out/nuc_annotation.gb rdna/${s}/${s}.rdna.gbf
    cp ${wd}/out/nuc_features.tsv rdna/${s}/${s}.nr_features.tsv
    """

    stub:
    """
    mkdir -p rdna/${s}/mitom_nuc logs/nr_annotate
    touch rdna/${s}/${s}.rdna.gbf rdna/${s}/${s}.nr_features.tsv logs/nr_annotate/${s}.log
    """
}

process HISTONE_ANNOTATE {
    // H1-H4 in the histone cluster with `mitom nuc annotate` (nuclear_markers.py annotate)
    tag "$s"
    cpus 4
    memory '2000 MB'
    conda { condaEnv('mitom') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(fasta)
    path db

    output:
    tuple val(s), path("histone/${s}/${s}.histone.gbf"), path("histone/${s}/mitom_nuc"), emit: gb
    path "logs/histone_annotate/${s}.log"

    script:
    def wd = "histone/${s}/mitom_nuc"
    def mitom = cfg().mitom_dir ?: "\$PWD/${db}/Mitom"
    def log = "logs/histone_annotate/${s}.log"
    """
    mkdir -p histone/${s} logs/histone_annotate
    if [ ! -s ${fasta} ]; then echo "no histone cluster" > ${log}; : > histone/${s}/${s}.histone.gbf; mkdir -p ${wd}; exit 0; fi
    rm -rf ${wd}
    python ${scriptsDir()}/nuclear_markers.py annotate -i ${fasta} -o ${wd} --kind histone --sample ${s} \\
        --mitom ${mitom} -t ${task.cpus} > ${log} 2>&1
    cp ${wd}/${s}.histone.gbf histone/${s}/${s}.histone.gbf
    """

    stub:
    """
    mkdir -p histone/${s}/mitom_nuc logs/histone_annotate
    touch histone/${s}/${s}.histone.gbf logs/histone_annotate/${s}.log
    """
}

process NUC_MARKERS {
    // rDNA and histones in the whole-genome assembly with `mitom nuc`: every copy plus one consensus per gene
    tag "$s"
    cpus 16
    memory '8000 MB'
    conda { condaEnv('mitom') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(scaffolds)
    path db

    output:
    tuple val(s), path("nuclear/${s}/${s}.nuc_features.tsv"), path("nuclear/${s}/${s}.nuc_consensus.tsv"),
        path("nuclear/${s}/mitom_nuc"), emit: markers
    path "logs/nuc_markers/${s}.log"

    script:
    def wd = "nuclear/${s}/mitom_nuc"
    def mitom = cfg().mitom_dir ?: "\$PWD/${db}/Mitom"
    def log = "logs/nuc_markers/${s}.log"
    """
    mkdir -p nuclear/${s} logs/nuc_markers
    if [ ! -s ${scaffolds} ]; then echo "no assembly" > ${log}; : > nuclear/${s}/${s}.nuc_features.tsv; : > nuclear/${s}/${s}.nuc_consensus.tsv; mkdir -p ${wd}; exit 0; fi
    rm -rf ${wd}
    python ${mitom}/mitom nuc annotate -i ${scaffolds} -o ${wd} --genes rdna,histone \\
        --threads ${task.cpus} > ${log} 2>&1
    cp ${wd}/nuc_features.tsv nuclear/${s}/${s}.nuc_features.tsv
    cp ${wd}/consensus.tsv nuclear/${s}/${s}.nuc_consensus.tsv
    """

    stub:
    """
    mkdir -p nuclear/${s}/mitom_nuc logs/nuc_markers
    touch nuclear/${s}/${s}.nuc_features.tsv nuclear/${s}/${s}.nuc_consensus.tsv logs/nuc_markers/${s}.log
    """
}

process NUC_COPIES {
    // copy number = k-mer coverage of the marker contigs / median coverage of BUSCO single-copy contigs
    tag "$s"
    memory '1000 MB'
    conda { condaEnv('python') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(features), path(consensus), path(busco), path(busco_dir, stageAs: 'busco_run')

    output:
    tuple val(s), path("nuclear/${s}/${s}.copy_number.tsv"), emit: tsv
    path "logs/nuc_copies/${s}.log"

    script:
    def spec = smkJson(input: [features: "${features}", consensus: "${consensus}", busco: "${busco}"],
        output: [tsv: "nuclear/${s}/${s}.copy_number.tsv"], params: [busco_dir: 'busco_run'],
        log: ["logs/nuc_copies/${s}.log"], wildcards: [sample: s], rule: 'nuc_copies')
    """
    mkdir -p nuclear/${s} logs/nuc_copies
    ${spec}
    smk_shim.py snakemake.json ${scriptsDir()}/nuc_copies.py
    """

    stub:
    """
    mkdir -p nuclear/${s} logs/nuc_copies
    printf 'gene\\tcopy_number\\n' > nuclear/${s}/${s}.copy_number.tsv; touch logs/nuc_copies/${s}.log
    """
}

process NUC_COMPARE {
    // GetOrganelle (rDNA unit, histone cluster) vs SPAdes-MitoM (copies over the whole assembly), per gene
    tag "$s"
    memory '1000 MB'
    conda { condaEnv('mitom') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(s), path(rdna, stageAs: 'rdna_nuc'), path(histone, stageAs: 'histone_nuc'), path(spades, stageAs: 'spades_nuc')

    output:
    tuple val(s), path("nuclear/${s}/${s}.nuc_compare.tsv"), emit: tsv
    path "logs/nuc_compare/${s}.log"

    script:
    """
    mkdir -p nuclear/${s} logs/nuc_compare
    python ${scriptsDir()}/nuclear_markers.py compare --sample ${s} --getorganelle-rdna rdna_nuc/out \\
        --getorganelle-histone histone_nuc --spades spades_nuc -o nuclear/${s}/${s}.nuc_compare.tsv > logs/nuc_compare/${s}.log 2>&1
    """

    stub:
    """
    mkdir -p nuclear/${s} logs/nuc_compare
    touch nuclear/${s}/${s}.nuc_compare.tsv logs/nuc_compare/${s}.log
    """
}

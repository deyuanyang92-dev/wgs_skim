// Topology check, BLAST species check, summary (workflow/rules/topology.smk, blast.smk, report.smk)
include { cfg; condaEnv; smkJson; scriptsDir; nestLinks; nestInputs; names } from './common'

process TOPOLOGY {
    // are circular / linear results real? junction-spanning reads after cutting inside nad5
    tag "$s"
    cpus 8
    memory '8000 MB'
    conda { condaEnv('topology') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    // fastas: the mitogenome (+ the rDNA when nr is on); annotation: the MitoM .gb or []
    tuple val(s), path(fastas, stageAs: 'fasta/*'), path(r1), path(r2), path(annotation, stageAs: 'ann/*'), path(status)

    output:
    tuple val(s), path("topology/${s}/${s}.topology.tsv"), emit: tsv
    path "topology/${s}/*.circular.fasta", optional: true
    path "logs/topology/${s}.log"

    script:
    def c = cfg()
    def spec = smkJson(
        input: [fasta: names(fastas), r1: "${r1}", r2: "${r2}",
                annotation: annotation ? names(annotation)[0] : [], status: "${status}"],
        output: [tsv: "topology/${s}/${s}.topology.tsv"], params: [pairs: c.topology_pairs, mode: c.topology_mode],
        log: ["logs/topology/${s}.log"], threads: task.cpus, wildcards: [sample: s], rule: 'topology')
    """
    mkdir -p topology/${s} logs/topology
    ${spec}
    smk_shim.py snakemake.json ${scriptsDir()}/topology.py
    """

    stub:
    """
    mkdir -p topology/${s} logs/topology
    printf 'sample\\tcontig\\tverdict\\n' > topology/${s}/${s}.topology.tsv; touch logs/topology/${s}.log
    """
}

process BLAST {
    // species check of the batch: mitogenome / rDNA queries against nt, TopN hits with GenBank taxonomy
    cpus 8
    memory '48000 MB'
    conda { condaEnv('blast') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    // meta: [[key (mt, gb, nr, nr_gb), sample, file name], ...] for files
    tuple val(meta), path(files)

    output:
    path "blast/blast_topn.tsv", emit: tsv
    path "blast/*"
    path "logs/blast.log"

    script:
    def c = cfg()
    def spec = smkJson(input: nestInputs(meta), output: [tsv: 'blast/blast_topn.tsv'],
        params: [db: c.blast_db, topn: c.blast_topn, email: c.email, samples: meta.collect { it[1] }.unique().sort(),
                 mt: c.blast_mt_list, nr: c.blast_nr_list],
        log: ['logs/blast.log'], threads: task.cpus, rule: 'blast')
    """
    mkdir -p blast logs
    ${nestLinks(meta)}
    ${spec}
    smk_shim.py snakemake.json ${scriptsDir()}/blast_taxonomy.py
    """

    stub:
    """
    mkdir -p blast logs
    printf 'qseqid\\tsacc\\tpident\\n' > blast/blast_topn.tsv; touch logs/blast.log
    """
}

process SUMMARY {
    // one row per sample (samples.tsv) and a self-contained HTML report
    memory '1000 MB'
    conda { condaEnv('python') }
    publishDir params.outdir, mode: params.publish_mode

    input:
    tuple val(meta), path(files)
    val samples

    output:
    path "summary/samples.tsv", emit: tsv
    path "summary/report.html"
    path "logs/summary.log"

    script:
    def spec = smkJson(input: nestInputs(meta, ['blast']),
        output: [tsv: 'summary/samples.tsv', html: 'summary/report.html'],
        params: [samples: samples, taxon: cfg().taxon], log: ['logs/summary.log'], rule: 'summary')
    """
    mkdir -p summary logs
    ${nestLinks(meta)}
    ${spec}
    smk_shim.py snakemake.json ${scriptsDir()}/summary.py
    """

    stub:
    """
    mkdir -p summary logs
    printf 'sample\\n${samples.join('\\n')}\\n' > summary/samples.tsv; touch summary/report.html logs/summary.log
    """
}

#!/usr/bin/env nextflow
// wgs_skim with Nextflow: the same steps, scripts, software environments and output layout as the Snakemake
// workflow (workflow/Snakefile). Options are the keys of config/config.yaml (see nextflow.config).
//   nextflow run nextflow/main.nf -profile conda --reads_dir /path/to/rawdata --outdir results -resume
// or: wgs_skim run --engine nextflow -r /path/to/rawdata
nextflow.enable.dsl = 2

include { cfg; steps; loadSamples; cleanGb } from './modules/common'
include { FASTP } from './modules/qc'
include { GETORGANELLE_DB; MEANGS; GETORGANELLE_MT; MT_ITERATE; MT_SELECT } from './modules/mito'
include { SPADES; BUSCO_LINEAGE; BUSCO; JELLYFISH; GENOMESCOPE } from './modules/assembly'
include { MITOM_DB; MITOM_NUC_DB; MITOM_ANNOTATE } from './modules/annotate'
include { GETORGANELLE_NUC as GETORGANELLE_NR; GETORGANELLE_NUC as GETORGANELLE_HISTONE; NR_ANNOTATE;
          HISTONE_ANNOTATE; NUC_MARKERS; NUC_COPIES; NUC_COMPARE } from './modules/nuclear'
include { TOPOLOGY; BLAST; SUMMARY } from './modules/report'

// Processes to run: the requested products plus everything they need (final_outputs / *_inputs in common.smk)
def needed(Map c, List st) {
    def nucDb = c.mitom_dir ? [] : ['MITOM_NUC_DB']
    def deps = [
        SUMMARY             : ['FASTP'] + (c.mito ? ['MT_SELECT'] : []) + (c.nr ? ['GETORGANELLE_NR'] : []) +
                              (c.assembly ? ['SPADES', 'BUSCO'] : []) + (c.assembly && c.nuc_markers ? ['NUC_COPIES'] : []) +
                              (c.genome_size ? ['GENOMESCOPE'] : []) + (c.annotate && c.mito ? ['MITOM_ANNOTATE'] : []) +
                              (c.topology_on && c.mito ? ['TOPOLOGY'] : []) + (c.blast && (c.mito || c.nr) ? ['BLAST'] : []),
        MEANGS              : ['FASTP'],
        GETORGANELLE_MT     : ['FASTP', 'GETORGANELLE_DB'] + (c.meangs ? ['MEANGS'] : []),
        MT_ITERATE          : ['FASTP', 'GETORGANELLE_MT', 'GETORGANELLE_DB'],
        MT_SELECT           : ['MT_ITERATE'] + (c.meangs ? ['MEANGS'] : []) + (c.assembly && c.mt_from_spades ? ['SPADES'] : []),
        SPADES              : ['FASTP'],
        BUSCO               : ['SPADES'] + (c.busco_lineage_dir ? [] : ['BUSCO_LINEAGE']),
        JELLYFISH           : ['FASTP'],
        GENOMESCOPE         : ['JELLYFISH'],
        MITOM_ANNOTATE      : ['MT_SELECT'] + (c.mitom_dir && c.mitos_refdir ? [] : ['MITOM_DB']),
        GETORGANELLE_NR     : ['FASTP'],
        GETORGANELLE_HISTONE: ['FASTP'],
        NR_ANNOTATE         : ['GETORGANELLE_NR'] + nucDb,
        HISTONE_ANNOTATE    : ['GETORGANELLE_HISTONE'] + nucDb,
        NUC_MARKERS         : ['SPADES'] + nucDb,
        NUC_COPIES          : ['NUC_MARKERS', 'BUSCO'],
        NUC_COMPARE         : ['NR_ANNOTATE', 'HISTONE_ANNOTATE', 'NUC_MARKERS'],
        TOPOLOGY            : ['MT_SELECT', 'FASTP'] + (c.nr ? ['GETORGANELLE_NR'] : []) + (c.annotate ? ['MITOM_ANNOTATE'] : []),
        BLAST               : (c.mito ? ['MT_SELECT'] : []) + (c.nr ? ['GETORGANELLE_NR'] : []) +
                              (c.mito && c.annotate && c.blast_mt_list.any { it != 'full' } ? ['MITOM_ANNOTATE'] : []) +
                              (c.nr && c.annotate && c.blast_nr_list.any { it != 'full' } ? ['NR_ANNOTATE'] : []),
    ]
    def targets
    if (st) {
        def byStep = [fastp      : ['FASTP'], meangs: ['MEANGS'], mito: ['MT_SELECT'], getorganelle: ['MT_SELECT'],
                      nr         : ['NR_ANNOTATE', 'HISTONE_ANNOTATE'], spades: ['SPADES'], busco: ['BUSCO'],
                      genome_size: ['GENOMESCOPE'], annotate: ['MITOM_ANNOTATE'], topology: ['TOPOLOGY'],
                      nuclear    : ['NUC_COPIES', 'NUC_COMPARE'], blast: ['BLAST'], summary: ['SUMMARY']]
        targets = st.collectMany { byStep[it] }
    } else {
        targets = ['SUMMARY'] + (c.mito ? ['MT_SELECT'] : []) +
                  (c.nr ? ['GETORGANELLE_NR', 'NR_ANNOTATE', 'HISTONE_ANNOTATE'] : []) +
                  (c.nr && c.assembly && c.nuc_markers ? ['NUC_COMPARE'] : [])
    }
    def need = [] as Set
    def add
    add = { p -> if (need.add(p)) (deps[p] ?: []).each { add(it) } }
    targets.each { add(it) }
    need
}

// `--steps summary` alone: the summary of an existing results directory, from the files that are there
def existingResults(Map c, List samples) {
    def out = file(params.outdir)
    def have = { pattern -> samples.every { s -> out.resolve(pattern.replace('{s}', s)).exists() } }
    c.mito = have('mito/{s}/{s}.mt.fasta')
    c.nr = have('rdna/{s}/{s}.nr.fasta')
    c.assembly = have('assembly/{s}/{s}.busco.txt')
    c.nuc_markers = have('nuclear/{s}/{s}.copy_number.tsv')
    c.genome_size = have('genome_size/{s}/{s}.genomescope.txt')
    c.annotate = have('annotation/{s}/{s}.mt.gb')
    c.topology_on = have('topology/{s}/{s}.topology.tsv')
    c.blast = out.resolve('blast/blast_topn.tsv').exists()
    def items = []
    def add = { key, pattern, on -> if (on) samples.each { s -> items << [key, s, out.resolve(pattern.replace('{s}', s))] } }
    add('fastp', 'qc/{s}/{s}.fastp.json', true)
    add('mt', 'mito/{s}/{s}.mt.json', c.mito)
    add('mt_fasta', 'mito/{s}/{s}.mt.fasta', c.mito)
    add('nr', 'rdna/{s}/{s}.nr.fasta', c.nr)
    add('scaffolds', 'assembly/{s}/{s}.scaffolds.fasta', c.assembly)
    add('busco', 'assembly/{s}/{s}.busco.txt', c.assembly)
    add('nuclear', 'nuclear/{s}/{s}.copy_number.tsv', c.assembly && c.nuc_markers)
    add('genomescope', 'genome_size/{s}/{s}.genomescope.txt', c.genome_size)
    add('annotation', 'annotation/{s}/{s}.mt.gb', c.annotate && c.mito)
    add('topology', 'topology/{s}/{s}.topology.tsv', c.topology_on && c.mito)
    if (c.blast && (c.mito || c.nr)) items << ['blast', '_', out.resolve('blast/blast_topn.tsv')]
    items
}

// [[key, sample, file], ...] -> one input tuple [meta, files], sorted so that -resume sees the same task
def collected(ch) {
    ch.toSortedList { a, b -> a[0] <=> b[0] ?: a[1] <=> b[1] }
      .map { l -> [l.collect { [it[0], it[1], it[2].name] }, l.collect { it[2] }] }
}

workflow {
    def c = cfg()
    def st = steps()
    def S = loadSamples()
    def names = S.keySet().sort()
    log.info "[wgs_skim] ${names.size()} samples: ${names.join(', ')}; results: ${file(params.outdir)}"

    if (st == ['summary']) {
        SUMMARY(collected(Channel.fromList(existingResults(c, names))), names)
        return
    }
    def need = needed(c, st)
    log.info "[wgs_skim] steps: ${need.sort().join(' ')}"
    def none = Channel.fromList(names).map { [it, []] }        // per-sample placeholder for an optional input

    def reads = Channel.fromList(names.collect { s -> [s, file(S[s][0]), file(S[s][1])] })
    FASTP(reads)
    def clean = FASTP.out.reads.map { s, r1, r2, j -> [s, r1, r2, j, cleanGb(j)] }   // + clean Gb (memory)
    def r12 = clean.map { s, r1, r2, j, gb -> [s, r1, r2] }

    // mitogenome
    def seed = none, mt = Channel.empty(), asm = Channel.empty(), ann = Channel.empty()
    def goDb = need.contains('GETORGANELLE_DB') ? GETORGANELLE_DB().collect() : Channel.empty()
    if (need.contains('MEANGS')) {
        MEANGS(r12)
    }
    if (need.contains('MT_SELECT')) {
        seed = c.meangs ? MEANGS.out.seed : Channel.fromList(names).map { [it, file("${projectDir}/../resources/no_seed.fasta")] }
        GETORGANELLE_MT(r12.join(seed), goDb)
        MT_ITERATE(clean.join(GETORGANELLE_MT.out.round0), goDb)
    }
    if (need.contains('SPADES')) {
        SPADES(clean.map { s, r1, r2, j, gb -> [s, r1, r2, j, gb, S[s][2] ?: ''] })
        asm = SPADES.out.asm
    }
    if (need.contains('MT_SELECT')) {
        def spadesForMt = c.assembly && c.mt_from_spades ? asm.map { s, f, d -> [s, [f, d]] } : none
        MT_SELECT(MT_ITERATE.out.result.join(seed).join(spadesForMt))
        mt = MT_SELECT.out.result
    }
    if (need.contains('MITOM_ANNOTATE')) {
        def db = need.contains('MITOM_DB') ? MITOM_DB().collect() : Channel.value([])
        MITOM_ANNOTATE(mt.map { s, f, j -> [s, f] }, db)
        ann = MITOM_ANNOTATE.out.gb
    }

    // whole-genome assembly, genome size
    if (need.contains('BUSCO')) {
        def lineage = c.busco_lineage_dir ? Channel.value(c.busco_lineage_dir) : BUSCO_LINEAGE(c.busco_lineage).collect()
        BUSCO(asm.map { s, f, d -> [s, f] }, lineage)
    }
    if (need.contains('GENOMESCOPE')) {
        JELLYFISH(r12)
        GENOMESCOPE(JELLYFISH.out.histo)
    }

    // nuclear markers
    def nucDb = need.contains('MITOM_NUC_DB') ? MITOM_NUC_DB().collect() : Channel.value([])
    def nrFasta = Channel.empty()
    if (need.contains('GETORGANELLE_NR')) {
        GETORGANELLE_NR(clean, Channel.value(['rdna', file(c.nr_seed), file(c.nr_label)]))
        nrFasta = GETORGANELLE_NR.out.fasta
    }
    if (need.contains('NR_ANNOTATE')) NR_ANNOTATE(nrFasta, nucDb)
    if (need.contains('HISTONE_ANNOTATE')) {
        GETORGANELLE_HISTONE(clean, Channel.value(['histone', file(c.histone_seed), file(c.histone_label)]))
        HISTONE_ANNOTATE(GETORGANELLE_HISTONE.out.fasta, nucDb)
    }
    if (need.contains('NUC_MARKERS')) NUC_MARKERS(asm.map { s, f, d -> [s, f] }, nucDb)
    if (need.contains('NUC_COPIES')) {
        NUC_COPIES(NUC_MARKERS.out.markers.map { s, f, cons, d -> [s, f, cons] }.join(BUSCO.out.summary))
    }
    if (need.contains('NUC_COMPARE')) {
        NUC_COMPARE(NR_ANNOTATE.out.gb.map { s, g, d -> [s, d] }
            .join(HISTONE_ANNOTATE.out.gb.map { s, g, d -> [s, d] })
            .join(NUC_MARKERS.out.markers.map { s, f, cons, d -> [s, d] }))
    }

    // topology, BLAST, summary
    if (need.contains('TOPOLOGY')) {
        def nrOpt = c.nr ? nrFasta : none
        def annOpt = c.annotate ? ann : none
        TOPOLOGY(mt.join(nrOpt).join(r12).join(annOpt)
            .map { s, f, j, nr, r1, r2, a -> [s, [f] + (nr ? [nr] : []), r1, r2, a, j] })
    }
    def blastTsv = Channel.empty()
    if (need.contains('BLAST')) {
        def q = Channel.empty()
        if (c.mito) q = q.mix(mt.map { s, f, j -> ['mt', s, f] })
        if (need.contains('MITOM_ANNOTATE') && c.blast_mt_list.any { it != 'full' }) q = q.mix(ann.map { s, g -> ['gb', s, g] })
        if (c.nr) q = q.mix(nrFasta.map { s, f -> ['nr', s, f] })
        if (need.contains('NR_ANNOTATE') && c.annotate && c.blast_nr_list.any { it != 'full' })
            q = q.mix(NR_ANNOTATE.out.gb.map { s, g, d -> ['nr_gb', s, g] })
        BLAST(collected(q))
        blastTsv = BLAST.out.tsv
    }
    if (need.contains('SUMMARY')) {
        def items = FASTP.out.reads.map { s, r1, r2, j -> ['fastp', s, j] }
        if (c.mito) items = items.mix(mt.flatMap { s, f, j -> [['mt', s, j], ['mt_fasta', s, f]] })
        if (c.nr) items = items.mix(nrFasta.map { s, f -> ['nr', s, f] })
        if (c.assembly) items = items.mix(asm.map { s, f, d -> ['scaffolds', s, f] },
                                          BUSCO.out.summary.map { s, f, d -> ['busco', s, f] })
        if (c.assembly && c.nuc_markers) items = items.mix(NUC_COPIES.out.tsv.map { s, f -> ['nuclear', s, f] })
        if (c.genome_size) items = items.mix(GENOMESCOPE.out.summary.map { s, f -> ['genomescope', s, f] })
        if (c.annotate && c.mito) items = items.mix(ann.map { s, g -> ['annotation', s, g] })
        if (c.topology_on && c.mito) items = items.mix(TOPOLOGY.out.tsv.map { s, f -> ['topology', s, f] })
        if (c.blast && (c.mito || c.nr)) items = items.mix(blastTsv.map { ['blast', '_', it] })
        SUMMARY(collected(items), names)
    }
}

workflow.onComplete {
    log.info "[wgs_skim] ${workflow.success ? 'done' : 'FAILED'}: results in ${file(params.outdir)}" +
             (workflow.success ? '' : " (fix the error, then run again with -resume: finished tasks are reused)")
}

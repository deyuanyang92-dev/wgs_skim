// Shared functions: the Nextflow counterpart of workflow/rules/common.smk (taxon presets, steps, memory, samples).
import groovy.json.JsonOutput
import groovy.json.JsonSlurper

def taxa() {
    [
        annelid   : [meangs_class: 'A-worms', genetic_code: 5, clade: 'Annelida-segmented-worms', busco_lineage: 'metazoa_odb12'],
        mollusc   : [meangs_class: 'Mollusca', genetic_code: 5, clade: 'Mollusca', busco_lineage: 'mollusca_odb12'],
        fish      : [meangs_class: 'Chordata', genetic_code: 2, clade: 'Chordata', busco_lineage: 'actinopterygii_odb12'],
        crustacean: [meangs_class: 'Arthropoda', genetic_code: 5, clade: 'Arthropoda', busco_lineage: 'arthropoda_odb12'],
    ]
}

def stepAliases() {
    [qc: 'fastp', 'meangs-getorganelle': 'mito', mitogenome: 'mito', rdna: 'nr', assembly: 'spades',
     jellyfish: 'genome_size', genomescope: 'genome_size', mitoz: 'annotate', mitom: 'annotate', nuc: 'nuclear',
     histone: 'nr', histones: 'nr', report: 'summary']
}

def stepsAll() {
    ['fastp', 'meangs', 'mito', 'getorganelle', 'nr', 'spades', 'busco', 'genome_size', 'annotate', 'topology',
     'blast', 'nuclear', 'summary']
}

// YAML / --key true|false|yes|no arrive as Boolean or String
def flag(v) {
    v instanceof Boolean ? v : (v?.toString()?.toLowerCase() in ['true', 'yes', '1', 'on'])
}

def csvList(v) {
    (v ?: '').toString().split(',').collect { it.trim() }.findAll { it }
}

def scriptsDir() {
    projectDir.parent.resolve('workflow/scripts')
}

def resourcesDir() {
    projectDir.parent.resolve('resources')
}

def steps() {
    def s = (params.steps ?: '').toString().trim().split(/[,\s]+/).findAll { it }.collect { stepAliases().get(it, it) }
    def bad = s.findAll { !(it in stepsAll()) }
    if (bad) error "unknown step '${bad[0]}'; choose from: ${stepsAll().join(', ')}"
    s
}

// params with the taxon preset and `steps` applied (as common.smk does to `config`)
def cfg() {
    def c = new LinkedHashMap(params)
    def preset = taxa()[c.taxon?.toString()]
    if (!preset) error "unknown taxon '${c.taxon}'; choose from: ${taxa().keySet().join(', ')}"
    preset.each { k, v -> if (c[k] in [null, '']) c[k] = v }
    ['meangs', 'mito', 'nr', 'assembly', 'nuc_markers', 'genome_size', 'annotate', 'blast', 'mt_from_spades'].each {
        c[it] = flag(c[it])
    }
    c.blast_mt_list = csvList(c.blast_mt)
    c.blast_nr_list = csvList(c.blast_nr)
    def topo = c.topology instanceof Boolean ? (c.topology ? 'yes' : 'no') : c.topology.toString().toLowerCase()
    c.topology_mode = [true: 'yes', false: 'no'].get(topo, topo)
    def want = steps() as Set
    if (want) {
        def genesWanted = (c.blast_mt_list + c.blast_nr_list).any { it != 'full' }
        c.annotate = !want.intersect(['annotate', 'topology']).isEmpty() || ('blast' in want && genesWanted)
        c.topology_mode = 'topology' in want ? 'yes' : 'no'
        c.blast = 'blast' in want
        c.mito = !want.intersect(['mito', 'getorganelle', 'annotate', 'topology', 'blast']).isEmpty()
        c.nr = !want.intersect(['nr', 'blast', 'nuclear']).isEmpty()
        c.assembly = !want.intersect(['spades', 'busco', 'nuclear']).isEmpty()
        c.nuc_markers = 'nuclear' in want
        c.genome_size = 'genome_size' in want
        if ('getorganelle' in want && want.intersect(['mito', 'meangs']).isEmpty())
            c.meangs = false                                   // GetOrganelle with its own animal_mt seeds
    }
    c.topology_on = c.topology_mode != 'no'
    c.db = file((c.db_dir ?: resourcesDir()).toString().replaceFirst(/^~/, System.getProperty('user.home')))
    c.nr_seed = c.nr_seed ?: resourcesDir().resolve('rdna/seed.fasta').toString()
    c.nr_label = c.nr_label ?: resourcesDir().resolve('rdna/label.fasta').toString()
    c.histone_seed = c.histone_seed ?: resourcesDir().resolve('histone/seed.fasta').toString()
    c.histone_label = c.histone_label ?: resourcesDir().resolve('histone/label.fasta').toString()
    def lin = c.busco_lineage.toString().replaceFirst(/^~/, System.getProperty('user.home'))
    c.busco_lineage_dir = file(lin).isDirectory() ? file(lin) : null      // null: downloaded once into db_dir
    c
}

// {sample: [r1, r2, scaffolds]} with workflow/scripts/samples.py (the same detection as the Snakemake workflow)
def loadSamples() {
    def code = '''
import json, sys
sys.path.insert(0, sys.argv[1])
from samples import find_pairs, read_sheet, sheet_columns
if sys.argv[2] == "sheet":
    pairs, given = read_sheet(sys.argv[3]), sheet_columns(sys.argv[3], ["scaffolds"])
else:
    pairs, given = find_pairs(sys.argv[3]), {}
print(json.dumps({s: [r1, r2, given.get(s, {}).get("scaffolds", "")] for s, (r1, r2) in pairs.items()}))
'''
    def (mode, path) = params.samples ? ['sheet', params.samples] : ['dir', params.reads_dir]
    if (!path) error "give the reads directory: --reads_dir /path/to/rawdata (or --samples sheet.tsv)"
    def p = [params.python, '-c', code, scriptsDir().toString(), mode, path.toString()].execute()
    def out = new StringBuilder(), err = new StringBuilder()
    p.waitForProcessOutput(out, err)
    if (p.exitValue() != 0) error "sample detection failed: ${err.toString().trim().readLines()[-1]}"
    def lines = out.toString().trim().readLines()
    if (lines.size() > 1) lines[0..-2].each { log.warn(it) }   // find_pairs prints its warnings (unpaired R1 files)
    new JsonSlurper().parseText(lines[-1]) as Map
}

// A Snakemake env of the same YAML already in the conda prefix (`wgs_skim install`) is used as is; otherwise
// Nextflow creates one from the YAML in conda.cacheDir.
def condaEnv(String name) {
    def yaml = projectDir.parent.resolve("workflow/envs/${name}.yaml")
    def prefix = System.getenv('WGS_SKIM_CONDA_PREFIX') ?: params.conda_prefix
    if (prefix && file(prefix).isDirectory()) {
        def text = yaml.text
        def hit = file(prefix).listFiles().find {
            it.name.endsWith('_.yaml') && it.text == text && file(it.toString()[0..-6]).resolve('env_setup_done').exists()
        }
        if (hit) return hit.toString()[0..-6]
    }
    yaml.toString()
}

def cleanGb(json) {
    try {
        new JsonSlurper().parse(json.toFile()).summary.after_filtering.total_bases / 1e9
    } catch (Exception e) {
        0.0d
    }
}

// GetOrganelle memory grows with the reads it uses: ~1 GB per Gb (measured, see common.smk)
def getorganelleMem(gb, cap, attempt) {
    def capGb = cap?.toString()?.isNumber() ? cap.toString().toDouble() : 0
    def used = capGb ? Math.min(gb as double, capGb) : gb as double
    "${(Math.max(8, 1.2 * used) * 1024 * attempt) as long} MB"
}

def nucMaxGb() {
    csvList(params.nuclear_max_gb).collect { it.toDouble() }.max()
}

// SPAdes -m is a hard limit: estimated from clean data, x1.5 on the retry
def spadesMem(gb, given, attempt) {
    if (given) return '1000 MB'
    def est = (gb as double) * (params.spades_gib_per_gb as double)
    "${(Math.max(params.spades_min_mem_gb as double, est) * Math.pow(1.5, attempt - 1) * 1024) as long} MB"
}

// shell lines that write the `snakemake` object of smk_shim.py
def smkJson(Map m) {
    def d = [threads: 1, wildcards: [:], params: [:], scriptdir: scriptsDir().toString()] + m
    "cat > snakemake.json <<'__SMK_JSON__'\n${JsonOutput.prettyPrint(JsonOutput.toJson(d))}\n__SMK_JSON__"
}

// collected inputs [[key, sample, file name], ...]: link them as <key>/<sample>/<name>, the per-sample layout that
// summary.py / blast_taxonomy.py read sample names from; scalar keys give one path instead of a list
def nestLinks(List meta) {
    meta.collect { k, s, n -> "mkdir -p '${k}/${s}' && ln -sf \"\$PWD/${n}\" '${k}/${s}/${n}'" }.join('\n')
}

def nestInputs(List meta, List scalar = []) {
    def d = [:]
    meta.each { k, s, n -> d.get(k, []) << "${k}/${s}/${n}".toString() }
    d.collectEntries { k, v -> [k, k in scalar ? v[0] : v] }
}

// staged paths (relative, with any stageAs directory) of a path input that may be one file or a list
def names(x) {
    (x instanceof Collection ? x : [x]).findAll { it }.collect { it.toString() }
}

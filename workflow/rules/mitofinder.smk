# MitoFinder: the mitochondrial contigs found in (and annotated from) an assembly, for every sample. Sources
# (mitofinder_assembly, several allowed; results per source: mitofinder/<source>/<sample>/):
#   spades      the SPAdes whole-genome assembly           (combination spades-mitofinder)
#   rna         the rnaSPAdes transcripts                  (combination rnaspades-mitofinder: mt genes from RNA-seq)
#   trinity     the Trinity transcripts                    (combination trinity-mitofinder)
#   given       an assembly made elsewhere (sheet column scaffolds, wgs_skim mitofinder -i)
#   megahit / metaspades   MitoFinder assembles the reads itself (the first mitofinder_max_pairs read pairs)
# Reference (-r): mitofinder_ref, else the sample's own MitoM annotation of the GetOrganelle mitogenome (as the
# user's mitofinder-batch runs did). MitoFinder 1.4.2 from a Singularity image (mitofinder_sif; bioconda 1.4.1
# runs its circularization check without arguments, so its "Circularization" is meaningless - docs/tools/MitoFinder.md).
MF_SOURCES = ["spades", "rna", "trinity", "given", "megahit", "metaspades"]
MF_TRANSCRIPTS = {"rna", "trinity"}


def mitofinder_sources(sample):
    """the sources MitoFinder runs on for this sample"""
    asked = [x.strip() for x in str(config["mitofinder_assembly"]).split(",") if x.strip()]
    for x in asked:
        if x not in MF_SOURCES:
            raise WorkflowError(f"unknown mitofinder_assembly '{x}'; choose from: {', '.join(MF_SOURCES)}")
    if not asked:                                       # automatic: what this run assembles
        want = set(STEPS)
        if sample in GIVEN_ASM:
            asked = ["given"]
        elif want & {"rnaspades", "transcriptome", "trinity"} and "spades" not in want:
            asked = [x for x, steps in (("rna", {"rnaspades", "transcriptome"}),
                                        ("trinity", {"trinity", "transcriptome"})) if want & steps]
        else:
            asked = ["spades"] if config["assembly"] else ["megahit"]
    if sample in GIVEN_ASM:                             # a given assembly stands in for the SPAdes one
        asked = ["given" if x == "spades" else x for x in asked]
    return list(dict.fromkeys(x for x in asked if x != "given" or sample in GIVEN_ASM))


def mitofinder_targets(pattern):
    return [pattern.format(src=src, sample=s) for s in SAMPLES for src in mitofinder_sources(s)]


def mitofinder_ref(wildcards):
    if config["mitofinder_ref"]:
        return config["mitofinder_ref"]
    if config["mito"] and config["annotate"]:
        return rules.mitom_annotate.output.gb.format(sample=wildcards.sample)
    raise WorkflowError("MitoFinder needs a reference mitogenome GenBank: --mitofinder-ref ref.gb "
                        "(or run it together with the mitogenome step: its MitoM annotation is the reference)")


def mitofinder_assembly(wildcards):
    s, src = wildcards.sample, wildcards.src
    if src in ("spades", "given"):
        return rules.spades.output.scaffolds.format(sample=s)
    if src == "rna":
        return rules.spades_mode.output.fasta.format(sample=s, mode="rna")
    return rules.trinity.output.fasta.format(sample=s)


def mitofinder_min_contig(src):
    """MitoFinder's 1000 bp suits genome assemblies; transcripts of single mitochondrial genes are shorter"""
    return config["mitofinder_min_contig_rna"] if src in MF_TRANSCRIPTS else config["mitofinder_min_contig"]


rule mitofinder_input:
    """what MitoFinder searches: a whole-genome assembly cut to its deep scaffolds (>= mitofinder_cov_factor x the
    weighted median coverage; the mitogenome is far deeper than the nuclear genome) plus the scaffolds the reference
    hits (low-coverage skims), transcripts by length only (scripts/mitofinder_prefilter.py)"""
    input:
        asm=mitofinder_assembly,
        ref=mitofinder_ref,
    output:
        temp(f"{OUT}mitofinder/{{src}}/{{sample}}/input.fasta"),
    params:
        min_len=lambda w: mitofinder_min_contig(w.src),
        factor=lambda w: 0 if w.src in MF_TRANSCRIPTS else config["mitofinder_cov_factor"],
    log:
        f"{OUT}logs/mitofinder/{{src}}/{{sample}}.input.log",
    threads: 4
    resources:
        mem_mb=4000,
    wildcard_constraints:
        src="spades|given|rna|trinity",
    conda:
        "../envs/blast.yaml"
    shell:
        "python {workflow.basedir}/scripts/mitofinder_prefilter.py {input.asm} {output} --min-len {params.min_len} "
        "--factor {params.factor} --ref {input.ref} --threads {threads} > {log} 2>&1"


def mitofinder_inputs(wildcards):
    s, src = wildcards.sample, wildcards.src
    d = {"ref": mitofinder_ref(wildcards)}
    if src in ("megahit", "metaspades"):
        d.update(r1=rules.fastp.output.r1.format(sample=s), r2=rules.fastp.output.r2.format(sample=s),
                 json=rules.fastp.output.json.format(sample=s))
    else:
        d["assembly"] = rules.mitofinder_input.output[0].format(src=src, sample=s)
    return d


def mitofinder_mem_mb(wildcards, input, attempt):
    if "assembly" in input.keys():
        return 8000 * attempt                      # annotation only: a BLAST of the assembly against the reference
    return int(config["mitofinder_mem_gb"] * 1024 * 1.5 ** (attempt - 1))


rule mitofinder:
    """MitoFinder on one sample and source -> <s>.mitofinder.fasta / .gb (all mitochondrial contigs), Final_Results/"""
    input:
        unpack(mitofinder_inputs),
    output:
        fasta=f"{OUT}mitofinder/{{src}}/{{sample}}/{{sample}}.mitofinder.fasta",
        gb=f"{OUT}mitofinder/{{src}}/{{sample}}/{{sample}}.mitofinder.gb",
    params:
        out=lambda w, output: str(Path(output.fasta).parent),
        data=lambda w, input: (f"-1 {input.r1} -2 {input.r2} --assembler {w.src}" if "r1" in input.keys()
                               else f"-a {input.assembly}"),
        pairs=config["mitofinder_max_pairs"],
        code=config["genetic_code"],
        trna=config["mitofinder_trna"],
        min_contig=lambda w: mitofinder_min_contig(w.src),
        sif=f"--sif {config['mitofinder_sif']}" if config["mitofinder_sif"] else "",
    resources:
        mem_mb=mitofinder_mem_mb,
    log:
        f"{OUT}logs/mitofinder/{{src}}/{{sample}}.log",
    benchmark:
        f"{OUT}benchmarks/mitofinder/{{src}}/{{sample}}.tsv"
    threads: 8
    wildcard_constraints:
        src="|".join(MF_SOURCES),
    conda:
        "../envs/mitofinder.yaml"
    shell:
        "{PY} {workflow.basedir}/scripts/mitofinder_run.py -s {wildcards.sample} -r {input.ref} -o {params.out} "
        "{params.data} --source {wildcards.src} --max-pairs {params.pairs} --code {params.code} --trna {params.trna} "
        "--min-contig {params.min_contig} -p {threads} -m $(( {resources.mem_mb} / 1024 )) {params.sif} --log {log}"


rule mitofinder_summary:
    """one source, every sample: contigs, genes, tRNAs, gene lengths; per-contig .infos; one FASTA per gene"""
    input:
        lambda w: [rules.mitofinder.output.fasta.format(src=w.src, sample=s) for s in SAMPLES
                   if w.src in mitofinder_sources(s)],
    output:
        tsv=f"{OUT}mitofinder/{{src}}/mitofinder_summary.tsv",
    params:
        dirs=lambda w, input: [str(Path(f).parent) for f in input],
        version="1.4.2 (Singularity)" if config["mitofinder_sif"] else "1.4.1 (bioconda: circularization not checked)",
    log:
        f"{OUT}logs/mitofinder/{{src}}/summary.log",
    resources:
        mem_mb=1000,
    wildcard_constraints:
        src="|".join(MF_SOURCES),
    conda:
        "../envs/python.yaml"
    shell:
        "python {workflow.basedir}/scripts/mitofinder_summary.py -o {output.tsv} --source {wildcards.src} "
        "--version '{params.version}' {params.dirs} > {log} 2>&1"


rule mitofinder_blast:
    """contamination check of one source: every gene of every sample against nt (the BLAST step's search and
    GenBank metadata, from blast2metadata); genes whose best hit is another family than the sample's other genes
    are flagged (transcriptomes: gut contents, parasites, cross-contaminated libraries)"""
    input:
        tsv=rules.mitofinder_summary.output.tsv,
    output:
        samples=f"{OUT}mitofinder/{{src}}/blast/samples_check.tsv",
        genes=f"{OUT}mitofinder/{{src}}/blast/genes_check.tsv",
        hits=f"{OUT}mitofinder/{{src}}/blast/blast_topn.tsv",
    params:
        genes_dir=lambda w, input: str(Path(input.tsv).parent / "genes"),
        out=lambda w, output: str(Path(output.samples).parent),
        genes=config["mitofinder_blast_genes"],
        db=f"--db {config['blast_db']}" if config["blast_db"] else "",
        email=f"--email {config['email']}" if config["email"] else "",
        topn=config["blast_topn"],
        min_pident=config["mitofinder_min_pident"],
    log:
        f"{OUT}logs/mitofinder/{{src}}/blast.log",
    benchmark:
        f"{OUT}benchmarks/mitofinder/{{src}}/blast.tsv"
    threads: 8
    resources:
        mem_mb=48000,                     # a local nt / marker database is memory-mapped (as rule blast)
    wildcard_constraints:
        src="|".join(MF_SOURCES),
    conda:
        "../envs/blast.yaml"
    shell:
        "python {workflow.basedir}/scripts/mitofinder_check.py --genes-dir {params.genes_dir} -o {params.out} "
        "--genes {params.genes} {params.db} {params.email} --topn {params.topn} --threads {threads} "
        "--min-pident {params.min_pident} > {log} 2>&1"

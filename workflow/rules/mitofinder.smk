# MitoFinder: the mitochondrial contigs found in (and annotated from) an assembly, for every sample. Sources
# (mitofinder_assembly, several allowed). Results per source, one directory each:
#   mitofinder_<source>/<sample>/   one sample's MitoFinder result
#   mitofinder_<source>/            mitofinder_summary.tsv, mitofinder_contigs.tsv, genes/, contigs/, contamination/
#   spades      the SPAdes whole-genome assembly           (combination spades-mitofinder)
#   rnaspades   the rnaSPAdes transcripts                  (combination rnaspades-mitofinder: mt genes from RNA-seq)
#   trinity     the Trinity transcripts                    (combination trinity-mitofinder)
#   assembly    an assembly made elsewhere (sheet column scaffolds, wgs_skim mitofinder -i)
#   megahit / metaspades   MitoFinder assembles the reads itself (the first mitofinder_max_pairs read pairs)
# Reference (-r): mitofinder_ref, else the sample's own MitoM annotation of the GetOrganelle mitogenome (as the
# user's mitofinder-batch runs did). MitoFinder 1.4.2 from a Singularity image (mitofinder_sif; bioconda 1.4.1
# runs its circularization check without arguments, so its "Circularization" is meaningless - docs/tools/MitoFinder.md).
MF_SOURCES = ["spades", "rnaspades", "trinity", "assembly", "megahit", "metaspades"]
MF_ALIAS = {"rna": "rnaspades", "given": "assembly"}            # names of wgs_skim <= 0.1
MF_TRANSCRIPTS = {"rnaspades", "trinity"}


def mitofinder_sources(sample):
    """the sources MitoFinder runs on for this sample"""
    asked = [MF_ALIAS.get(x.strip(), x.strip()) for x in str(config["mitofinder_assembly"]).split(",") if x.strip()]
    for x in asked:
        if x not in MF_SOURCES:
            raise WorkflowError(f"unknown mitofinder_assembly '{x}'; choose from: {', '.join(MF_SOURCES)}")
    if not asked:                                       # automatic: what this run assembles
        want = set(STEPS)
        if sample in GIVEN_ASM:
            asked = ["assembly"]
        elif want & {"rnaspades", "transcriptome", "trinity"} and "spades" not in want:
            asked = [x for x, steps in (("rnaspades", {"rnaspades", "transcriptome"}),
                                        ("trinity", {"trinity", "transcriptome"})) if want & steps]
        else:
            asked = ["spades"] if config["assembly"] else ["megahit"]
    if sample in GIVEN_ASM:                             # a given assembly stands in for the SPAdes one
        asked = ["assembly" if x == "spades" else x for x in asked]
    return list(dict.fromkeys(x for x in asked if x != "assembly" or sample in GIVEN_ASM))


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
    if src in ("spades", "assembly"):
        return scaffolds_of(s)
    if src == "rnaspades":
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
        fasta=temp(f"{OUT}mitofinder_{{src}}/{{sample}}/input.fasta"),
        places=temp(f"{OUT}mitofinder_{{src}}/{{sample}}/input.fasta.map.tsv"),
    params:
        min_len=lambda w: mitofinder_min_contig(w.src),
        factor=lambda w: 0 if w.src in MF_TRANSCRIPTS else config["mitofinder_cov_factor"],
        # transcripts: candidates joined into one sequence (scripts/mitofinder_prefilter.py; false = each transcript
        # on its own, as MitoFinder alone)
        concat=lambda w: (f"--concat --code {config['genetic_code']} --name {w.sample}_mt_transcripts"
                          if w.src in MF_TRANSCRIPTS and config["mitofinder_concat_transcripts"] else ""),
    log:
        f"{OUT}pipeline_info/logs/mitofinder_{{src}}/{{sample}}.input.log",
    threads: 4
    resources:
        mem_mb=4000,
    wildcard_constraints:
        src="spades|assembly|rnaspades|trinity",
    conda:
        "../envs/blast.yaml"
    shell:
        "python {workflow.basedir}/scripts/mitofinder_prefilter.py {input.asm} {output.fasta} --min-len {params.min_len} "
        "--factor {params.factor} --ref {input.ref} --threads {threads} {params.concat} > {log} 2>&1"


def mitofinder_inputs(wildcards):
    s, src = wildcards.sample, wildcards.src
    d = {"ref": mitofinder_ref(wildcards)}
    if src in ("megahit", "metaspades"):
        d.update(clean_reads(wildcards))
    else:
        d["assembly"] = rules.mitofinder_input.output.fasta.format(src=src, sample=s)
        d["places"] = rules.mitofinder_input.output.places.format(src=src, sample=s)
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
        fasta=f"{OUT}mitofinder_{{src}}/{{sample}}/{{sample}}.mitofinder.fasta",
        gb=f"{OUT}mitofinder_{{src}}/{{sample}}/{{sample}}.mitofinder.gb",
    params:
        out=lambda w, output: str(Path(output.fasta).parent),
        data=lambda w, input: (f"-1 {input.r1} -2 {input.r2} --assembler {w.src}" if "r1" in input.keys()
                               else f"-a {input.assembly} --places {input.places}"),
        pairs=config["mitofinder_max_pairs"],
        code=config["genetic_code"],
        trna=config["mitofinder_trna"],
        min_contig=lambda w: mitofinder_min_contig(w.src),
        sif=f"--sif {config['mitofinder_sif']}" if config["mitofinder_sif"] else "",
    resources:
        mem_mb=mitofinder_mem_mb,
    log:
        f"{OUT}pipeline_info/logs/mitofinder_{{src}}/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/mitofinder_{{src}}/{{sample}}.tsv"
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
        tsv=f"{OUT}mitofinder_{{src}}/mitofinder_summary.tsv",
    params:
        dirs=lambda w, input: [str(Path(f).parent) for f in input],
        version="1.4.2 (Singularity)" if config["mitofinder_sif"] else "1.4.1 (bioconda: circularization not checked)",
    log:
        f"{OUT}pipeline_info/logs/mitofinder_{{src}}/summary.log",
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
        samples=f"{OUT}mitofinder_{{src}}/contamination/samples_check.tsv",
        genes=f"{OUT}mitofinder_{{src}}/contamination/genes_check.tsv",
        hits=f"{OUT}mitofinder_{{src}}/contamination/blast_topn.tsv",
    params:
        genes_dir=lambda w, input: str(Path(input.tsv).parent / "genes"),
        out=lambda w, output: str(Path(output.samples).parent),
        genes=config["mitofinder_blast_genes"],
        db=f"--db {config['blast_db']}" if config["blast_db"] else "",
        email=f"--email {config['email']}" if config["email"] else "",
        topn=config["blast_topn"],
        min_pident=config["mitofinder_min_pident"],
        cache=str(DB / "genbank_metadata.tsv"),
    log:
        f"{OUT}pipeline_info/logs/mitofinder_{{src}}/blast.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/mitofinder_{{src}}/blast.tsv"
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
        "--min-pident {params.min_pident} --cache {params.cache} > {log} 2>&1"

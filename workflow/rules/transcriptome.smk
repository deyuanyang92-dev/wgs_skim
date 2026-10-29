# SPAdes in any of its modes and Trinity, with BUSCO on every assembly - the batch assemblies of the user's
# spadesv3.2.py / spadesv3.5.5.py (wgs, rna, meta, bio, sewage, plasmid, metaviral, metaplasmid; --careful,
# --only-assembler, -k) and trinity_batch3.1.py (--full_cleanup). The genome mode is rule `spades` (assembly.smk:
# its scaffolds feed BUSCO, the nuclear markers and the mitogenome); every other mode is rule spades_mode:
#   config spades_modes: "genome,rna,meta"  ->  assembly/<s>/<s>.scaffolds.fasta, <s>.rna.fasta, <s>.meta.fasta
# BUSCO (rule busco_assembly): -m transcriptome for rna / rnaviral / trinity, -m genome otherwise.
# Steps: spades (the modes), trinity, busco (all assemblies), rna = rnaSPAdes + Trinity + BUSCO.

# SPAdes mode -> (flag, result file in the SPAdes directory, BUSCO mode); results as spadesv3.5.5.py lists them
SPADES_MODES = {
    "isolate": ("--isolate", "scaffolds.fasta", "genome"),
    "sc": ("--sc", "scaffolds.fasta", "genome"),
    "meta": ("--meta", "scaffolds.fasta", "genome"),
    "rna": ("--rna", "transcripts.fasta", "transcriptome"),
    "rnaviral": ("--rnaviral", "transcripts.fasta", "transcriptome"),
    "plasmid": ("--plasmid", "scaffolds.fasta", "genome"),
    "metaplasmid": ("--metaplasmid", "scaffolds.fasta", "genome"),
    "metaviral": ("--metaviral", "scaffolds.fasta", "genome"),
    "bio": ("--bio", "gene_clusters.fasta", "genome"),
    "sewage": ("--sewage", "scaffolds.fasta", "genome"),
    "corona": ("--corona", "scaffolds.fasta", "genome"),
}
BUSCO_MODE = {m: v[2] for m, v in SPADES_MODES.items()} | {"trinity": "transcriptome"}
for m in OTHER_MODES:
    if m not in SPADES_MODES:
        raise WorkflowError(f"unknown SPAdes mode '{m}'; choose from: genome, {', '.join(SPADES_MODES)}")


def mode_mem_mb(gib_per_gb, floor_gb):
    """mem_mb from the clean data (fastp.json) x gib_per_gb, at least floor_gb, x1.5 on each retry"""
    def f(wildcards, input, attempt):
        return int(max(floor_gb, clean_gb(input.json) * gib_per_gb) * 1.5 ** (attempt - 1) * 1024)
    return f


def spades_kmers(mode):
    k = str(config["rnaspades_kmers"] if mode in ("rna", "rnaviral") else config["spades_kmers"]).strip()
    return f"-k {k}" if k and k != "auto" else ""


def spades_extra(mode):
    opts = [("--careful", config["spades_careful"] and mode in ("isolate", "sc", "plasmid")),
            ("--only-assembler", config["spades_only_assembler"])]
    return " ".join(x for x, on in opts if on)


rule spades_mode:
    """SPAdes in a non-genome mode (rna, meta, plasmid ...)"""
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
        json=rules.fastp.output.json,
    output:
        fasta=f"{OUT}assembly/{{sample}}/{{sample}}.{{mode}}.fasta",
    params:
        out=lambda w, output: str(Path(output.fasta).parent / f"spades_{w.mode}"),
        flag=lambda w: SPADES_MODES[w.mode][0],
        result=lambda w: SPADES_MODES[w.mode][1],
        k=lambda w: spades_kmers(w.mode),
        extra=lambda w: spades_extra(w.mode),
    resources:
        mem_mb=mode_mem_mb(config["spades_mode_gib_per_gb"], 16),
    log:
        f"{OUT}logs/spades_mode/{{sample}}.{{mode}}.log",
    benchmark:
        f"{OUT}benchmarks/spades_mode/{{sample}}.{{mode}}.tsv"
    threads: 16
    wildcard_constraints:
        mode="|".join(SPADES_MODES),
    conda:
        "../envs/spades.yaml"
    shell:
        """
        mem=$(( {resources.mem_mb} / 1024 ))
        if [ -f {params.out}/params.txt ]; then
            spades.py -o {params.out} --restart-from last -t {threads} -m $mem > {log} 2>&1
        else
            spades.py {params.flag} -1 {input.r1} -2 {input.r2} -o {params.out} -t {threads} -m $mem \
                {params.k} {params.extra} > {log} 2>&1
        fi
        if [ -s {params.out}/{params.result} ]; then cp {params.out}/{params.result} {output.fasta}
        else echo "[wgs_skim] SPAdes {wildcards.mode}: no {params.result} (nothing found)" >> {log}; : > {output.fasta}; fi
        rm -rf {params.out}/K* {params.out}/tmp {params.out}/misc {params.out}/corrected   # intermediates
        """


rule trinity:
    """Trinity transcripts (--full_cleanup: only Trinity.fasta and its gene-to-transcript map stay)"""
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
        json=rules.fastp.output.json,
    output:
        fasta=f"{OUT}assembly/{{sample}}/{{sample}}.trinity.fasta",
    params:
        out=lambda w, output: str(Path(output.fasta).parent / "trinity"),
        extra=config["trinity_extra"],
    resources:
        mem_mb=mode_mem_mb(config["trinity_gib_per_gb"], 16),
    log:
        f"{OUT}logs/trinity/{{sample}}.log",
    benchmark:
        f"{OUT}benchmarks/trinity/{{sample}}.tsv"
    threads: 16
    conda:
        "../envs/trinity.yaml"
    shell:
        """
        mem=$(( {resources.mem_mb} / 1024 ))G
        Trinity --seqType fq --left {input.r1} --right {input.r2} --max_memory $mem --CPU {threads} \
            --output {params.out}.Trinity_tmp --full_cleanup {params.extra} > {log} 2>&1
        mv {params.out}.Trinity_tmp.Trinity.fasta {output.fasta}
        mv {params.out}.Trinity_tmp.Trinity.fasta.gene_trans_map {output.fasta}.gene_trans_map 2>/dev/null || true
        """


rule busco_assembly:
    """BUSCO completeness of an assembly of another mode (-m transcriptome for RNA-seq assemblies)"""
    input:
        fasta=f"{OUT}assembly/{{sample}}/{{sample}}.{{asm}}.fasta",
        lineage=busco_lineage_dir,
    output:
        summary=f"{OUT}assembly/{{sample}}/{{sample}}.{{asm}}.busco.txt",
    params:
        out=lambda w, output: str(Path(output.summary).parent),
        name=lambda w: f"busco_{w.asm}",
        mode=lambda w: BUSCO_MODE[w.asm],
    resources:
        mem_mb=16000,
    log:
        f"{OUT}logs/busco_assembly/{{sample}}.{{asm}}.log",
    benchmark:
        f"{OUT}benchmarks/busco_assembly/{{sample}}.{{asm}}.tsv"
    threads: 16
    wildcard_constraints:
        asm="|".join(BUSCO_MODE),
    conda:
        "../envs/busco.yaml"
    shell:
        """
        if [ ! -s {input.fasta} ]; then echo "no assembly" > {output.summary}; exit 0; fi
        rm -rf {params.out}/{params.name}
        busco -i {input.fasta} -l {input.lineage} -m {params.mode} -o {params.name} --out_path {params.out} \
            -c {threads} --offline -f > {log} 2>&1
        cp {params.out}/{params.name}/short_summary.*.txt {output.summary}
        """

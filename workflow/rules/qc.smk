# Reads. qc: true -> fastp cleans them once (fastp/<sample>/); every later tool reads the clean files from there.
# qc: false (the reads are already clean, e.g. a fastp/ directory of an earlier run) -> no copy, no fastp/ directory:
# the tools read the given files by their absolute paths; read_stats only counts them (the data size sets SPAdes /
# GetOrganelle / Trinity memory and the GetOrganelle data cap) into pipeline_info/read_stats/.


def clean_reads(wildcards):
    """r1, r2 (the reads every tool uses) and json (their read / base counts, fastp format)"""
    s = wildcards.sample
    if config["qc"]:
        return {"r1": rules.fastp.output.r1.format(sample=s), "r2": rules.fastp.output.r2.format(sample=s),
                "json": rules.fastp.output.json.format(sample=s)}
    r1, r2 = SAMPLES[s]
    return {"r1": os.path.abspath(r1), "r2": os.path.abspath(r2), "json": rules.read_stats.output.json.format(sample=s)}


def reads_json(sample):
    return clean_reads(type("W", (), {"sample": sample})())["json"]


rule fastp:
    """fastp cleaning: adapters (paired-end detection), quality, length"""
    input:
        unpack(reads),
    output:
        r1=f"{OUT}fastp/{{sample}}/{{sample}}.R1.fq.gz",
        r2=f"{OUT}fastp/{{sample}}/{{sample}}.R2.fq.gz",
        html=f"{OUT}fastp/{{sample}}/{{sample}}.fastp.html",
        json=f"{OUT}fastp/{{sample}}/{{sample}}.fastp.json",
    resources:
        mem_mb=2000,
    log:
        f"{OUT}pipeline_info/logs/fastp/{{sample}}.log",
    benchmark:
        f"{OUT}pipeline_info/benchmarks/fastp/{{sample}}.tsv"
    threads: 8
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --detect_adapter_for_pe "
        "-w {threads} -h {output.html} -j {output.json} > {log} 2>&1"


rule read_stats:
    """qc: false - read and base counts of the reads as given (fastp with every filter off writes no reads)"""
    input:
        unpack(reads),
    output:
        json=f"{OUT}pipeline_info/read_stats/{{sample}}/{{sample}}.json",
    resources:
        mem_mb=2000,
    log:
        f"{OUT}pipeline_info/logs/read_stats/{{sample}}.log",
    threads: 4
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp -i {input.r1} -I {input.r2} -A -G -Q -L -w {threads} -h /dev/null -j {output.json} > {log} 2>&1"

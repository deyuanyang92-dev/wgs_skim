rule fastp:
    input:
        unpack(reads),
    output:
        r1=f"{OUT}qc/{{sample}}/{{sample}}.R1.fq.gz",
        r2=f"{OUT}qc/{{sample}}/{{sample}}.R2.fq.gz",
        html=f"{OUT}qc/{{sample}}/{{sample}}.fastp.html",
        json=f"{OUT}qc/{{sample}}/{{sample}}.fastp.json",
    resources:
        mem_mb=2000,
    log:
        f"{OUT}logs/fastp/{{sample}}.log",
    benchmark:
        f"{OUT}benchmarks/fastp/{{sample}}.tsv"
    threads: 8
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --detect_adapter_for_pe "
        "-w {threads} -h {output.html} -j {output.json} > {log} 2>&1"

rule fastqc:
    input:
        r1="results/fastq/{sample}_R1.fastq.gz",
        r2="results/fastq/{sample}_R2.fastq.gz",
    output:
        zip1="results/qc/fastqc/{sample}_R1_fastqc.zip",
        zip2="results/qc/fastqc/{sample}_R2_fastqc.zip",
    params:
        outdir="results/qc/fastqc",
    threads: 2
    log:
        "logs/fastqc/{sample}.log",
    conda:
        "../envs/qc.yaml"
    shell:
        "fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} &> {log}"


rule fastp:
    input:
        r1="results/fastq/{sample}_R1.fastq.gz",
        r2="results/fastq/{sample}_R2.fastq.gz",
    output:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        json="results/qc/fastp/{sample}.fastp.json",
        html="results/qc/fastp/{sample}.fastp.html",
    threads: 4
    log:
        "logs/fastp/{sample}.log",
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            --detect_adapter_for_pe -w {threads} \
            -j {output.json} -h {output.html} &> {log}
        """


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_R{mate}_fastqc.zip", sample=SAMPLES, mate=[1, 2]),
        expand("results/qc/fastp/{sample}.fastp.json", sample=SAMPLES),
        expand("results/salmon/{sample}/quant.sf", sample=SAMPLES),
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc -f -o results/qc -n multiqc_report.html results/qc results/salmon &> {log}"

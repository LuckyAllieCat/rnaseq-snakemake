rule salmon_quant:
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        index="resources/salmon_index",
    output:
        "results/salmon/{sample}/quant.sf",
    params:
        outdir="results/salmon/{sample}",
    threads: 4
    log:
        "logs/salmon/{sample}.log",
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} \
            --seqBias --gcBias -p {threads} -o {params.outdir} &> {log}
        """

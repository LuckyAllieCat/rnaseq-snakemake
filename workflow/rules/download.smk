rule download_subsample:
    output:
        r1="results/fastq/{sample}_R1.fastq.gz",
        r2="results/fastq/{sample}_R2.fastq.gz",
    params:
        url1=lambda wc: ena_url(SAMPLE_TABLE[wc.sample]["run"], 1),
        url2=lambda wc: ena_url(SAMPLE_TABLE[wc.sample]["run"], 2),
        lines=config["subsample_reads"] * 4,
    log:
        "logs/download/{sample}.log",
    conda:
        "../envs/download.yaml"
    shell:
        # `head` closes the pipe early on purpose, so curl's SIGPIPE is ignored
        # and the result is checked afterwards instead.
        """
        (curl -sSL --retry 3 {params.url1} | gzip -dc | head -n {params.lines} || true) \
            | gzip > {output.r1} 2> {log}
        (curl -sSL --retry 3 {params.url2} | gzip -dc | head -n {params.lines} || true) \
            | gzip > {output.r2} 2>> {log}
        n1=$(gzip -dc {output.r1} | wc -l)
        n2=$(gzip -dc {output.r2} | wc -l)
        echo "R1 lines: $n1, R2 lines: $n2" >> {log}
        [ "$n1" -gt 0 ] && [ "$n1" -eq "$n2" ]
        """

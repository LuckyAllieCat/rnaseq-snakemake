rule get_transcriptome:
    output:
        "resources/transcripts.fa.gz",
    params:
        url=config["reference"]["transcriptome_url"],
    log:
        "logs/reference/get_transcriptome.log",
    conda:
        "../envs/download.yaml"
    shell:
        "curl -sSL --retry 3 -o {output} {params.url} 2> {log}"


rule tx2gene:
    input:
        "resources/transcripts.fa.gz",
    output:
        "resources/tx2gene.tsv",
    # GENCODE headers: >ENST|ENSG|OTTHUMG|OTTHUMT|tx_name|gene_name|length|biotype|
    shell:
        """
        gzip -dc {input} | grep '^>' | sed 's/^>//' \
            | awk -F'|' 'BEGIN {{OFS="\\t"; print "transcript_id","gene_id","gene_name"}} {{print $1,$2,$6}}' \
            > {output}
        """


rule salmon_index:
    input:
        "resources/transcripts.fa.gz",
    output:
        directory("resources/salmon_index"),
    threads: 8
    log:
        "logs/reference/salmon_index.log",
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon index -t {input} -i {output} --gencode -p {threads} &> {log}"

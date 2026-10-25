rule deseq2:
    input:
        quant=expand("results/salmon/{sample}/quant.sf", sample=SAMPLES),
        tx2gene="resources/tx2gene.tsv",
        samples=config["samples"],
    output:
        results="results/deseq2/deseq2_results.tsv",
        counts="results/deseq2/normalized_counts.tsv",
        dds="results/deseq2/dds.rds",
        pca="results/deseq2/pca.png",
        ma="results/deseq2/ma_plot.png",
        volcano="results/deseq2/volcano.png",
        heatmap="results/deseq2/top_genes_heatmap.png",
    params:
        design=config["deseq2"]["design"],
        contrast=config["deseq2"]["contrast"],
        alpha=config["deseq2"]["alpha"],
    log:
        "logs/deseq2.log",
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2.R"


rule enrichment:
    input:
        "results/deseq2/deseq2_results.tsv",
    output:
        ora="results/enrichment/go_bp_ora.tsv",
        ora_plot="results/enrichment/go_bp_ora.png",
        gsea="results/enrichment/hallmark_gsea.tsv",
        gsea_plot="results/enrichment/hallmark_gsea.png",
    params:
        padj=config["enrichment"]["padj_cutoff"],
        lfc=config["enrichment"]["lfc_cutoff"],
    log:
        "logs/enrichment.log",
    conda:
        "../envs/enrichment.yaml"
    script:
        "../scripts/enrichment.R"

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# --- inputs -----------------------------------------------------------------
samples <- read.delim(snakemake@input[["samples"]], comment.char = "#",
                      stringsAsFactors = FALSE)
rownames(samples) <- samples$sample

files <- snakemake@input[["quant"]]
names(files) <- basename(dirname(files))
samples <- samples[names(files), , drop = FALSE]

design <- as.formula(snakemake@params[["design"]])
for (v in all.vars(design)) samples[[v]] <- factor(samples[[v]])

contrast <- unlist(snakemake@params[["contrast"]])
alpha <- snakemake@params[["alpha"]]

tx2gene <- read.delim(snakemake@input[["tx2gene"]])
gene_names <- unique(tx2gene[, c("gene_id", "gene_name")])

# --- transcript -> gene, DESeq2 --------------------------------------------
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, 1:2])
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = design)

# Drop genes with too few counts to be testable.
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds)
saveRDS(dds, snakemake@output[["dds"]])

res <- results(dds, contrast = contrast, alpha = alpha)
summary(res)
# Shrunken fold changes for ranking/plots; keep the Wald stat for GSEA.
shr <- lfcShrink(dds, contrast = contrast, res = res, type = "ashr")

out <- data.frame(gene_id = rownames(res),
                  baseMean = res$baseMean,
                  log2FoldChange = shr$log2FoldChange,
                  lfcSE = shr$lfcSE,
                  stat = res$stat,
                  pvalue = res$pvalue,
                  padj = res$padj)
out <- merge(gene_names, out, by = "gene_id", all.y = TRUE)
out <- out[order(out$padj), ]
write.table(out, snakemake@output[["results"]], sep = "\t", quote = FALSE,
            row.names = FALSE)

norm <- counts(dds, normalized = TRUE)
write.table(data.frame(gene_id = rownames(norm), norm, check.names = FALSE),
            snakemake@output[["counts"]], sep = "\t", quote = FALSE,
            row.names = FALSE)

# --- plots ------------------------------------------------------------------
vsd <- vst(dds, blind = FALSE)

pca <- plotPCA(vsd, intgroup = all.vars(design)) + theme_bw()
ggsave(snakemake@output[["pca"]], pca, width = 6, height = 4.5, dpi = 150)

png(snakemake@output[["ma"]], width = 900, height = 700, res = 150)
plotMA(shr, alpha = alpha, main = paste(contrast[2], "vs", contrast[3]))
dev.off()

vol <- out[!is.na(out$padj), ]
vol$significant <- vol$padj < alpha
volcano <- ggplot(vol, aes(log2FoldChange, -log10(padj), colour = significant)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) +
  geom_text(data = head(vol, 15), aes(label = gene_name),
            colour = "black", size = 2.5, vjust = -0.6, check_overlap = TRUE) +
  labs(title = paste(contrast[2], "vs", contrast[3])) +
  theme_bw()
ggsave(snakemake@output[["volcano"]], volcano, width = 6, height = 5, dpi = 150)

top <- head(out$gene_id[!is.na(out$padj)], 30)
mat <- assay(vsd)[top, ]
mat <- mat - rowMeans(mat)
rownames(mat) <- gene_names$gene_name[match(top, gene_names$gene_id)]
pheatmap(mat, annotation_col = as.data.frame(colData(vsd)[, all.vars(design), drop = FALSE]),
         filename = snakemake@output[["heatmap"]], width = 7, height = 8)

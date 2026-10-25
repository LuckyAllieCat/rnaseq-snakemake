log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(ggplot2)
})
set.seed(42)

res <- read.delim(snakemake@input[[1]], stringsAsFactors = FALSE)
res <- res[!is.na(res$padj), ]
res$ensembl <- sub("\\..*$", "", res$gene_id)  # strip GENCODE version suffix

padj_cut <- snakemake@params[["padj"]]
lfc_cut <- snakemake@params[["lfc"]]

placeholder_plot <- function(path, msg) {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = msg) + theme_void()
  ggsave(path, p, width = 6, height = 3, dpi = 150)
}

# --- over-representation: GO biological process -----------------------------
sig <- res$ensembl[res$padj < padj_cut & abs(res$log2FoldChange) > lfc_cut]
cat(length(sig), "significant genes used for ORA\n")

ego <- enrichGO(gene = sig, universe = res$ensembl, OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH",
                readable = TRUE)
ora <- if (is.null(ego)) data.frame() else as.data.frame(ego)
write.table(ora, snakemake@output[["ora"]], sep = "\t", quote = FALSE,
            row.names = FALSE)
if (nrow(ora) > 0) {
  ggsave(snakemake@output[["ora_plot"]], dotplot(ego, showCategory = 15),
         width = 8, height = 6, dpi = 150)
} else {
  placeholder_plot(snakemake@output[["ora_plot"]], "No enriched GO terms")
}

# --- GSEA: MSigDB Hallmark, genes ranked by DESeq2 Wald statistic -----------
# msigdbr >= 10 renamed `category` to `collection`.
hallmark <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr(species = "Homo sapiens", category = "H")
)
term2gene <- unique(hallmark[, c("gs_name", "gene_symbol")])

ranked <- res[!is.na(res$gene_name) & !is.na(res$stat), ]
ranked <- ranked[order(-abs(ranked$stat)), ]
ranked <- ranked[!duplicated(ranked$gene_name), ]
gene_list <- sort(setNames(ranked$stat, ranked$gene_name), decreasing = TRUE)

gsea <- GSEA(gene_list, TERM2GENE = term2gene, pvalueCutoff = 1,
             seed = TRUE, verbose = FALSE)
gsea_df <- as.data.frame(gsea)
write.table(gsea_df, snakemake@output[["gsea"]], sep = "\t", quote = FALSE,
            row.names = FALSE)

top <- head(gsea_df[gsea_df$p.adjust < 0.25, ], 20)
if (nrow(top) > 0) {
  top$Description <- sub("^HALLMARK_", "", top$Description)
  p <- ggplot(top, aes(NES, reorder(Description, NES), fill = NES > 0)) +
    geom_col() +
    scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "steelblue"),
                      labels = c(`TRUE` = "up", `FALSE` = "down"), name = NULL) +
    labs(x = "Normalized enrichment score", y = NULL,
         title = "Hallmark GSEA (FDR < 0.25)") +
    theme_bw()
  ggsave(snakemake@output[["gsea_plot"]], p, width = 7, height = 5, dpi = 150)
} else {
  placeholder_plot(snakemake@output[["gsea_plot"]], "No Hallmark sets at FDR < 0.25")
}

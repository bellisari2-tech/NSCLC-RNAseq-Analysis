# Differential expression analysis: NSCLC tumor vs normal

library(DESeq2)
library(ggplot2)

#Loading cleaned data from script 01
counts_filtered <- readRDS("data/processed/counts_filtered.rds")
condition <- readRDS("data/processed/condition.rds")
annot <- readRDS("data/processed/annot.rds")

#Build DESeq2 object and run analysis
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = data.frame(condition = condition),
                              design = ~ condition)
dds <- DESeq(dds)
red <- results(dds)
summary(res)

#QC: PCA (VST-transformed data)
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "condition")

#QC: MA plot (raw LFC, diagnostic only)
plotMA(res, ylim = c(-5, 5))

#Shrinking l2fc (moderates unstable estimates)
res_shrunk <- lfcShrink(dds, coef = "condition_tumor_vs_normal", type = "apeglm")
summary(res_shrunk)

#Annotate results with gene symbols
res_shrunk_df <- as.data.frame(res_shrunk)
res_shrunk_df$GeneID <- rownames(res_shrunk_df)
res_shrunk_annotated <- merge(res_shrunk_df, annot[, c("GeneID", "Symbol")], by = "GeneID")

#Volcano plot (shrunk LFC)
res_shrunk_annotated$significant <- ifelse(
  !is.na(res_shrunk_annotated$padj) & res_shrunk_annotated$padj < 0.05 & abs(res_shrunk_annotated$log2FoldChange) > 1, "yes", "no"
)

ggplot(res_shrunk_annotated, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.05, size = 1) +
  scale_color_manual(values = c("no" = "grey70", "yes" = "firebrick")) +
  theme_minimal() +
  labs(title = "Volcano plot (shrunk LFC): NSCLC tumor vs normal",
       x = "log2 Fold Change (shrunk)", y = "-log10(p-value)")

#Save annotated results table for downstream analysis
dir.create("data/processed", showWarnings = FALSE)
write.csv(res_shrunk_annotated, "data/processed/deseq2_results_annotated.csv", row.names = FALSE)
saveRDS(dds, "data/processed/dds.rds")

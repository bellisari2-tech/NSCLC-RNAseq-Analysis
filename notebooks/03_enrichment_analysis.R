#Checking possible normalization problem
library(DESeq2)
dds <- readRDS("data/processed/dds.rds")
summary(sizeFactors(dds))

condition <- readRDS("data/processed/condition.rds")
tapply(sizeFactors(dds), condition, summary)

#Everything ok

#loading DESeq2 results
res <- read.csv("data/processed/deseq2_results_annotated.csv")

#Applying filter
sig_strict <- subset(res, !is.na(padj) & padj < 0.01 & abs(log2FoldChange) > 2)
nrow(sig_strict)
up_genes <- subset(sig_strict, log2FoldChange > 0)
down_genes <- subset(sig_strict, log2FoldChange < 0)
nrow(up_genes)
nrow(down_genes)

#Functional enrichment analysis: up/down genes (NSCLC)
library(gprofiler2)

up_symbols <- up_genes$Symbol
down_symbols <- down_genes$Symbol

gost_up <- gost(query = up_symbols, organism = "hsapiens", sources = c("GO:BP", "KEGG"),
                custom_bg = res$Symbol)
gost_down <- gost(query = down_symbols, organism = "hsapiens", sources = c("GO:BP", "KEGG"),
                  custom_bg = res$Symbol)

head(gost_up$result[order(gost_up$result$p_value), c("term_name", "p_value", "source")], 15)
head(gost_down$result[order(gost_down$result$p_value), c("term_name", "p_value", "source")], 15)

#Enrichment plot for visualization

p_up <- gostplot(gost_up, capped = TRUE, interactive = FALSE)

p_down <- gostplot(gost_down, capped = TRUE, interactive = FALSE)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
p_up_final <- publish_gostplot(p_up, highlight_terms = head(gost_up$result[order(gost_up$result$p_value), "term_id"], 10)) +
  ggtitle("Enrichment: Upregulated genes (NSCLC tumor vs normal)") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.background = element_rect(fill = "white", color = NA))

p_down_final <- publish_gostplot(p_down, highlight_terms = head(gost_down$result[order(gost_down$result$p_value), "term_id"], 10)) +
  ggtitle("Enrichment: Downregulated genes (NSCLC tumor vs normal)") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.background = element_rect(fill = "white", color = NA))

ggsave("results/figures/enrichment_up_plot.png", plot = p_up_final, width = 12, height = 7, bg = "white")
ggsave("results/figures/enrichment_down_plot.png", plot = p_down_final, width = 12, height = 7, bg = "white")


#Saving results
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(gost_up$result[, colnames(gost_up$result) != "parents"], 
          "results/tables/enrichment_up.csv", row.names = FALSE)

write.csv(gost_down$result[, colnames(gost_down$result) != "parents"], 
          "results/tables/enrichment_down.csv", row.names = FALSE)

#04 Validation checks

#Verifying merge() integrity for gene ID's
res_raw <- read.csv("data/processed/deseq2_results_annotated.csv")
nrow(res_raw)

annot <- readRDS("data/processed/annot.rds")
nrow(annot)
length(unique(annot$GeneID))

#Checking some known marker genes
marker_genes <- c("EGFR", "KRT5", "KRT14", "NAPSA", "MKI67", "TP63")

res_raw[res_raw$Symbol %in% marker_genes, c("Symbol", "log2FoldChange", "padj")]

#Checking DESeq2 aumatic filter
library(DESeq2)
dds <- readRDS("data/processed/dds.rds")
res_original <- results(dds)

metadata(res_original)$filterThreshold

#Checking multiple test correction method
?gost

#Saving sessioninfo()
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

#
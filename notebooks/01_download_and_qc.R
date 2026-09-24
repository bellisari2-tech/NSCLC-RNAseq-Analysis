# 01 Download and QC: GSE81089 (NSCLC vs Normal tissue)

library(GEOquery)

#Loading Raw Counts
counts <- read.delim("data/raw/GSE81089_raw_counts_GRCh38.p13_NCBI.tsv",
                     row.names = 1, check.names = FALSE)
dim(counts) #Expected 39376 genes x 218 samples


#Loading annotation (Entrez ID -> Symbol)
annot <- read.delim("data/raw/Human.GRCh38.p13.annot.tsv",
                    check.names = FALSE)

#Metadata from GEO (Tumor vs Regular)
gse <- getGEO("GSE81089", GSEMatrix = TRUE)
pheno <- pData(gse[[1]])

#'source_name_ch1' is reliable for tumor/normal
#('tumor (t) or normal (n):ch1 is not: it contains patient ID)
table(pheno$source_name_ch1)
#expected: 199 "Human NSCLC tissue", 19 "Human non-malignant tissue"

#Verifying alignment sample-counts <-> pheno
stopifnot(all(colnames(counts) == rownames(pheno)))

#Creating condition variables
condition <- factor(pheno$source_name_ch1,
                    levels = c("Human non-malignant tissue", "Human NSCLC tissue"),
                    labels = c("normal", "tumor"))
table(condition)

#Filtering out low expression genes
keep <- rowSums(counts >= 10) >= 10
counts_filtered <- counts[keep, ]
dim(counts_filtered)

#Saving clean files -> next script
dir.create("data/processed", showWarnings = FALSE)
saveRDS(counts_filtered, "data/processed/counts_filtered.rds")
saveRDS(pheno, "data/processed/pheno.rds")
saveRDS(condition, "data/processed/condition.rds")
saveRDS(annot, "data/processed/annot.rds")

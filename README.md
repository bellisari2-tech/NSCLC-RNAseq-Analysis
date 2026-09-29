# NSCLC RNA-seq Differential Expression & Functional Enrichment Analysis

Differential gene expression and functional enrichment analysis of RNA-seq data from Non-Small Cell Lung Cancer (NSCLC) tumor tissue vs adjacent normal tissue.

## Background

This project analyzes public RNA-seq data to identify genes and biological processes that differ between NSCLC tumor and normal lung tissue, using a standard bulk RNA-seq differential expression workflow (DESeq2) followed by functional enrichment analysis (GO Biological Process, KEGG).

## Data

- **Source**: [GSE81089](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81089) 
  (Djureinovic et al.), NCBI GEO
- **Samples**: 218 total — 199 NSCLC tumor, 19 matched/adjacent normal lung tissue
- **Data type**: bulk RNA-seq raw gene-level counts (NCBI-generated, GRCh38.p13), 
  processed with TopHat 2.0.8b / Cufflinks 2.1.1 by the original submitters

## Pipeline

Three sequential R scripts ('notebooks/'):

1. '01_download_and_qc.R' - Load raw counts and gene annotation, retrieve sample metadata via GEOquery, verify sample alignment between counts matrix and metadata, filter low expression genes (>= 10 reads in >= 10 samples).

2. '02_deseq2_analysis.R - Differential expression analysis with DESeq2, PCA and MA-plot QC, log-fold-change shrinkage (apeglm), gene symbol annotation, volcano plot.

3. '03_enrichment_analysis.R' - Functional enrichment via g:Profiler, run separately on up- and downregulated genes (padj < 0.01, |log2FC| > 2), using the set of tested genes as the statistical background.

Environment details: see 'results/sessionInfo.txt'

#Key findings

- Differential expression: 3084 genes pass strict thresholds - 2581 upregulated, 503 downregulated in tumor vs normal.

- Histological heterogeneity in tumor samples: PCA on variance-stabilized expression shows that normal samples cluster tightly, while tumor samples form multiple distinct sub-clusters. this structure aligns with the histological subtype annotation available in the sample metadata (three unlabeled numeric codes in GEO), suggesting NSCLC histological heterogeneity - not batch effects or QC issues - drives part of the transciptomic variability.

- Upregulated in tumor -> keratinization/cornification program: The top enriched terms are dominated by keratinization, cornified envelope formation, and epidermal/epithelial differentiaton. This is independently supported by strong, highly significant upregulation of classic squamous-lineage markers such as KRT5, KRT14, and TP63 (all padj < 1e-12), and is plausibly linked to the histological heterogeneity observed in the PCA - consistent with a squamous-cell-carcinoma component in the cohort, though this specific link was not directly verified against per-sample histology labels.

- Downregulated in tumor -> vascular development terms: enrichment is dominated by angiogenesis and blood vessel development/morphogenesis terms. This is counterintuitive given the well established role of tumor-induced angiogenesis in cancer, and is more likely explained by two non-exclusive, factors rather than reduced vascularization per se: (1) a cell-composition effect typical of bulk RNA-seq, where stromal/endothelial signal is diluted by the dominant mass of proliferating tumor epithelial cells; (2) tumor vasculature, while angiogenic, is structually immature and disorganized compare to healthy vasculature, and may not activate the same "normal developmental" gene programs captured by these GO terms. Confirming this would require single-cell RNA-seq or computational deconvolution.

- Validation via known markers: MKI67 (proliferation marker) is significantly upregulated, consistent with increased tumor proliferation. EGFR shows modest but significant upregulation, consistent with its nown but heterogenous role in NSCLC (large efect EGFR alterations are typically genomic mutation/amplification rather than simple expression changes). NAPSA (adenocarcinoma marker) shows a non significant trend, consistent with histological heterogeneity diluting to a subtype specific signal in an aggregated tumor-vs-normal comparison.

## Limitations

- Unpaired design: although 19 patients have both tumor and matched normal samples, the full analysis uses an unpaired '~ condition' design (most patients lack a matched normal sample). A paired design ('~ patient + condition') on the matched subset was not implemented but would control for inter-patient variability. 
- Imbalanced group sizes: 199 tumor vs 19 normal, reduces precision of the normal-group baseline estimate.
- Histological heterogeneity is not modeled: the '~ condition' design treats all tumor histological subtypes as one group; histology-stratified analysis was not performed.
- Bulk RNA-seq cell composition confounding: as discussed above for the vascular signal, bulk tissue expression reflects a mixture of cell types whose relative proportions differ between tumor and normal tissue, which can affect interpretation of stromal/immune/vascular gene signals. Undocumented numeric histology codes in GEO metadata prevented direct, quantitative confirmation of the histology-keratinization link.


## Reproducing this analysis

1. Open 'lung_cancer_project.Rproj' in RStudio
2. Download the raw data manually from the GEO section (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81089) into 'data/raw/' (raw counts matrix and gene annotation table).
3. Run the three scripts in 'notebooks/' in numerical order.
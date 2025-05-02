# bulk-RNA-seq-pipeline
Scripts for bulk RNA sequencing analysis using FastQC, HISAT2, HTSeq, DESeq2, and more.
bulk-RNA-seq-pipeline/
│
├── 1_fastqc/                     # FastQC quality check
├── 2_alignment/                 # Alignment with HISAT2
├── 3_counting/                  # Gene-level quantification with HTSeq
├── 4_differential_expression/   # DESeq2-based differential expression analysis

Deseq2_workflow.R: This script performs normalization of raw count data using the median-of-ratios method implemented in DESeq2. It generates principal component analysis (PCA) plots for sample clustering and outputs differential expression results as CSV files for downstream analysis.

Deseq2_Workflow_2.sh: This script organizes and processes DESeq2 output files for downstream analysis. It moves quality control (QC) plots and volcano plots into dedicated folders, converts DESeq2 result CSVs into tab-delimited text files, extracts normalized read counts, and trims results to retain only statistical data. It appends gene symbols to DESeq2 results using HTSeq counts, filters for statistically significant genes (p ≤ 0.05), and further separates them into upregulated (positive log2FC) and downregulated (negative log2FC) gene lists. Finally, it creates clean lists of gene symbols only for easy integration into other analyses. This modular and organized approach ensures efficient post-DESeq2 analysis and publication-ready data organization.

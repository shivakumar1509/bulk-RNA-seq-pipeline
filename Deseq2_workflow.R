#Load the DESeq2 library into R
library('DESeq2')

# Set your working directory to a directory containing your counts files. It's a good idea to create a directory containing only counts files
setwd('/Users/shivakumar06/Desktop/7_Deseq2_ready_2_30SR')

# Make a variable called directory and specify the entire path to your counts files 
directory <- "/Users/shivakumar06/Desktop/7_Deseq2_ready_2_30SR"

# Make a variable called sampleFiles. Grep all files with the ".txt" extention from the "directory" variable ( the directory variable specifies a path to counts files). This will combine your individual count files produced from HTSeq-count for example, into a single data matrix. For DESeq 1, you have to combine your count files manually. DESeq does this for you.
sampleFiles <- grep(".txt",list.files(directory),value=TRUE)

# View contents of sampleFiles
sampleFiles

# Specify the sample conditions for the newly created samplesFiles and store as a new variable named "sampleCondition". In this example, I have 4 different sample conditions (qc, stmv, tmv, and uninfected).
sampleCondition<-c("Bl6", "Bl6", "Bl6",
                  "R2054Q", "R2054Q","R2054Q")


# Create a new sample table for DESeq2. Notice that we use the "sampleFiles" and "sampleCondition" variables created previously.
sampleTable <- data.frame(sampleName=sampleFiles, fileName=sampleFiles,condition=sampleCondition)

# View contents of sampleTable
sampleTable

# Make a data table from the HTSeq count files using the sampleTable information, specifying the counts directory and setting the table design to condition.
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory,design=~condition)

# Optional. Count Genes with non-zero in all samples.
GeneCounts <- counts(ddsHTSeq)
idx.nz <- apply(GeneCounts,1,function(x){all(x>0)})
sum(idx.nz)

# Estimate Size factor
ddsHTSeqsf<-estimateSizeFactors(ddsHTSeq)
#ddsHTSeqsf<-estimateSizeFactors(ddsHTSeq, controlGenes=11548)
sizeFactors(ddsHTSeqsf)

# Load gplots library
library("ggplot2")

# Load geneplotter library
library('geneplotter')

# Load the RColorBrewer library
library("RColorBrewer")

# Plot densities of counts for the different samples to assess their distribution
color = colorRampPalette(rev(brewer.pal(n = 9, name = "Set1")))(24)
pdf("count_density_allsamples.pdf", width=11, height=7.5)
par(oma=c(1,1,1,1), mar = c(5,5,2,17))
multiecdf(counts(ddsHTSeqsf, normalized=T)[idx.nz,], col=color, xlab="Mean Counts", xlim=c(0,200), main = "Count Density Plot",legend = list("topright", inset=c(-.47,0), legend=sampleFiles, fill=color, title="Samples", bg = "white", xpd=TRUE))
dev.off()

pdf("allsamples_multidensity_allsamples.pdf", width=11, height=7.5)
par(oma=c(1,1,1,1), mar = c(5,5,2,17))
multidensity(counts(ddsHTSeqsf, normalized=T)[idx.nz,], col=color, xlab="mean counts", xlim=c(0,100), main = "Multidensity Plot",legend = list("topright", inset=c(-.47,0), legend=sampleFiles, fill=color, title="Samples", bg = "white", xpd=TRUE))
dev.off()

# Designate the uninfected sample group as the control
colData(ddsHTSeqsf)$condition <- relevel(ddsHTSeqsf$condition, ref="Bl6")

# Run DESeq2 which will do the following
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
dds <-DESeq(ddsHTSeqsf)

# The rlogTransformation function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size. The rlog transformation produces a similar variance stabilizing effect as ΓÇÿvarianceStabilizingTransformationΓÇÖ, though ΓÇÿrlogΓÇÖ is more robust in the case when the size factors vary widely. The transformation is useful when checking for outliers or as input for machine learning techniques such as clustering or linear discriminant analysis.
rld <-rlogTransformation(dds,blind=TRUE)
rld

# This function calculates a variance stabilizing transformation (VST) from the fitted dispersion-mean relation(s) and then transforms the count data (normalized by division by the size factors or normalization factors), yielding a matrix of values which are now approximately homoskedastic (having constant variance along the range of mean values). The transformation also normalizes with respect to library size.
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vsd

# Load the gplots library
library("gplots")

# Heatmaps for 50 most highly expressed genes (not necessarily the biggest fold change). Change the values in [1:50] to increase or decrease the number of genes in the heatmap
# The data is of┬áraw counts (*.heatmap1)
hmcol <- colorRampPalette(brewer.pal(9,"Reds"))(100)
pdf("50_most_expressed_genes_DESeq2_heatmap1_allsamples.pdf", width=9, height=11)
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
heatmap.2(counts(dds,normalized=TRUE)[select,],col=hmcol, Rowv=FALSE, Colv=FALSE, scale="none", dendrogram="none",trace="none", margin=c(12,8), lhei=c(2,8), lwid=c(2,4))
dev.off()

pdf("100_most_expressed_genes_DESeq2_heatmap1_allsamples.pdf", width=9, height=22)
select2 <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:100]
heatmap.2(counts(dds,normalized=TRUE)[select2,],col=hmcol, Rowv=FALSE, Colv=FALSE, scale="none", dendrogram="none",trace="none", margin=c(12,8), lhei=c(2,16), lwid=c(2,4))
dev.off()


# The data is from the regularized log transformation (rld) (*.heatmap2)
pdf("50_most_rld_genes_DESeq2_heatmap2_allsamples.pdf", width=9, height=11)
heatmap.2(assay(rld)[select,],col=hmcol,Rowv=FALSE, Colv=FALSE, scale="none",dendrogram="none", trace="none",margin=c(12,8), lhei=c(2,8), lwid=c(2,4))
dev.off()

pdf("100_most_rld_genes_DESeq2_heatmap2_allsamples.pdf", width=9, height=22)
heatmap.2(assay(rld)[select2,],col=hmcol,Rowv=FALSE, Colv=FALSE, scale="none",dendrogram="none", trace="none",margin=c(12,8), lhei=c(2,16), lwid=c(2,4))
dev.off()

# The data is from the variance stabilizing transformation (*.heatmap3)
pdf("50_most_vsd_genes_DESeq2_heatmap3_allsamples.pdf", width=9, height=11)
heatmap.2(assay(vsd)[select,],col=hmcol,Rowv=FALSE, Colv=FALSE, scale="none",dendrogram="none", trace="none",margin=c(12,8), lhei=c(2,8), lwid=c(2,4))
dev.off()

pdf("100_most_vsd_genes_DESeq2_heatmap3_allsamples.pdf", width=9, height=22)
heatmap.2(assay(vsd)[select2,],col=hmcol,Rowv=FALSE, Colv=FALSE, scale="none",dendrogram="none", trace="none",margin=c(12,8), lhei=c(2,16), lwid=c(2,4))
dev.off()

# Calculate the sample to sample distances using the rld transformed data to make a dendrogram to look at the clustering of the samples
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), 
                                       paste(condition,sampleFiles,sep=":"))
hc <- hclust(distsRL)
pdf("deseq2_sample_clustering_allsamples.pdf", width=11, height=11)
heatmap.2(mat,Rowv=as.dendrogram(hc), symm=TRUE,trace="none", col=rev(hmcol),margin=c(16,16))
dev.off()

# Perform a Principal Component Analysis and print the PCA plot from the rld transformed data.

# First PCA plot (rld)
pcavar <- plotPCA(rld, intgroup = c("condition"))
pcavar$data$name <- gsub(".txt", "", pcavar$data$name)

plot.of.interst <- ggplot(pcavar$data, aes(x = PC1, y = PC2, colour = group, label = name)) +
  geom_point(aes(size = 3)) +
  ggrepel::geom_text_repel() +
  xlim(-70, 70) +
  ylim(-70, 70) +
  theme_minimal()  # Optional for better aesthetics

pdf("PCA_plots_rld.pdf", width = 11, height = 6)
print(plot.of.interst)
dev.off()

# Second PCA plot (vsd)
pcavar2 <- plotPCA(vsd, intgroup = c("condition"))
pcavar2$data$name <- gsub(".txt", "", pcavar2$data$name)

plot.of.interst2 <- ggplot(pcavar2$data, aes(x = PC1, y = PC2, colour = group, label = name)) +
  geom_point(aes(size = 3)) +
  ggrepel::geom_text_repel() +
  xlim(-70, 70) +
  ylim(-70, 70) +
  theme_minimal()  # Optional for better aesthetics

pdf("PCA_plots_vsd.pdf", width = 11, height = 6)
print(plot.of.interst2)
dev.off()

###########################################################################
#Reorganize Files
dir.create("../8_Volcano_and_QC_Plots")
dir.create("../8_Volcano_and_QC_Plots/QC_Plots")
dir.create("../8_Volcano_and_QC_Plots/Volcano_Plots")
dir.create("../9_Deseq2_Results_CSV")
setwd('/Users/shivakumar6/Desktop/9_Deseq2_Results_CSV/')
###########################################################################
diff_exp_anal_function <- function(GroupA=A, GroupB=B) {
  
  # Differential Expression Analysis for comparing the qc sample to the uninfected control sample. Inspect the counts files to confirm the direction differential expression.
  #Control or reference library shoulbe listed second!
  ref_GroupA_GroupB <- results(dds, contrast=c("condition", GroupA, GroupB))
  
  #generate results file and add normalized read counts to the results file
  resdata.GroupA_GroupB <- merge(as.data.frame(ref_GroupA_GroupB), as.data.frame(counts(dds, normalized=TRUE)),by="row.names", sort=FALSE)
  
  #Look at the file
  head (resdata.GroupA_GroupB)
  
  #Print file
  write.csv(resdata.GroupA_GroupB, file= paste(GroupA, "_vs_", GroupB, ".deseq2.csv", sep=""))
  
  # Make MA plots of the DE results and print or save to a pdf file.  DEG padj. <= 0.1. This is the default. Change the alpha value to change the level of significance e.g., plotMA(res_qc, alpha = 0.05, ylim=c(-2,2),main="DESeq2")
  pdf(file= paste(GroupA, "_vs_", GroupB,".pdf"), width = 7, height = 5)
  #plotMA(ref_GroupA_GroupB,ylim=c(-10,10),main="DESeq2", xlab = "Mean of Normalized Counts")
  plotMA(ref_GroupA_GroupB,ylim=c(-5,5), xlab = "Mean of Normalized Counts", ylab = expression(paste("Log"[2]," Fold-Change")))
  dev.off()
}


A <- "R2054Q"
B <- "Bl6"
diff_exp_anal_function()


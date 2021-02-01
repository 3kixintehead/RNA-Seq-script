# DE Analysis of RNA-seq data.
# File:    GSE125554_zikv_counts.csv
# Source:  Molecular Alterations of Extracellular Matrix in the Brain of Newborns with Congenital Zika Syndrome
# Author:  Bryce Watson
# Date:    01-18-2021

#Developed with guidance from Simon Cockell: https://github.com/sjcockell/lockdown-learning
#Developed with guidance from DESeq2 Vignette: https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Load Libraries

library(BiocManager)
library(DESeq2)
library(Biobase)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(apeglm)
library(pheatmap)
library(genefilter)
library(plotly)
library(tibble)
library(rmarkdown)

#Import count data
setwd("~/R/R_Projects/DESeq2") # Set path to project location.
zika_counts <- file.path("./Zika_DESeq/GSE125554_zika_counts.csv")
cts <- read.csv(zika_counts, row.names=1, stringsAsFactors=FALSE)

# DeSeq Parameters
countdata <- as.matrix(cts)
zika_meta <- file.path("./Zika_DESeq/zika_metadata.csv")
col_data <- read.csv(zika_meta)

# Make DeSeqDataset
dds = DESeqDataSetFromMatrix(countData = cts,
                             colData = col_data,
                             design = ~ Condition)

#dds prefiltering
ddsf <- dds[ rowSums(counts(dds)) > 1, ]

#DeSeq Dataset

dds <- DESeq(ddsf)
res <- results(dds)
res_df <- as.data.frame(res) #Converts to data.frame instead of DeSeq Dataframe

#Normalization
vst = varianceStabilizingTransformation(dds)
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind=FALSE) # Used optionally instead of vsd.

## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"

#Boxplot
boxplot(assay(vst))

#Examine plot of p-values
hist(res$pvalue, breaks=40, col='gray')

#Shrink Log-fold change
resultsNames(dds) # No need to execute unless getting new coef for lfcShrink().
res <- lfcShrink(dds, coef="Condition_infected_vs_control", type="apeglm")

#MA Plot
plotMA(res, alpha=0.05, main='MA Plot', ylim=c(-14,6), cex=.5,
       colNonSig = "gray",
       colSig = "orangered",
       colLine = "grey60")

#PCA
plotPCA(vst, intgroup='Condition') +
  geom_text_repel(aes(label=name)) +
  ggtitle("Principle Component Analysis") +
  theme(legend.position="none", plot.title = element_text(size = rel(1.5), hjust = 0.5))
  
#Genes of Interest Annotation
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 25)
mat  <- assay(vsd)[ topVarGenes, ]
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mat  <- mat - rowMeans(mat)
gns <- getBM(c("hgnc_symbol", "ensembl_gene_id"), "ensembl_gene_id", row.names(mat), mart=mart, useCache = FALSE)
row.names(mat)[match(gns[,2], row.names(mat))] <- gns[,1]
#Make Heatmap
anno <- as.data.frame(colData(vsd)[, c("Sample","Condition")])
pheatmap(mat, annotation_col = anno)

#Filtering for more precise stat measures (used in volcanoplot)
filter_df <- res_df[complete.cases(res_df),] #If every column has data, returns TRUE
filter_df=filter_df[order(filter_df$padj),]
filter_df$test = filter_df$padj < 0.05 & abs(filter_df$log2FoldChange) > 1


filter_df_padj <- filter_df[filter_df$padj < 0.05,] #Filters by p-value < 0.05
filter_df_padj=filter_df_padj[order(filter_df_padj$padj),]

#Find most significant genes
topSigGenes <- filter_df_padj[1:20,]
filter_df_padj <- data.frame(filter_df_padj[,-1], row.names=filter_df_padj[,1])
gns_padj <- getBM(c("hgnc_symbol", "ensembl_gene_id"), "ensembl_gene_id", row.names(filter_df_padj), mart=mart, useCache =FALSE)
filter_df_padj$rownames <- gns_padj$hgnc_symbol

#GGplot Volcanoplot
filter_df_pval = rownames_to_column(filter_df_pval, var='Gene')
filter_df = rownames_to_column(filter_df, var='Gene')

ggplot(filter_df, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color=test), size=1.5, alpha=0.4) +
  scale_color_manual(values=c('violetred', 'gray', 'orangered')) +
  xlim(-15, 20) +
  ggtitle('Volcano Plot') +
  labs(y=expression('-Log'[10]*' P'[adj]), x=expression('Log'[2]*' fold change')) +
  geom_text_repel(data=topSigGenes, force=5,aes(x = log2FoldChange, y = -log10(padj),label=Gene))+
  geom_point(data=topSigGenes,aes(x = log2FoldChange, y = -log10(padj), color='black',alpha=0.4))+
  theme_minimal() +
  theme(legend.position="none", plot.title = element_text(size = rel(1.5), hjust = 0.5))

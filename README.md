# RNA-Seq-script
Script developed as an exploratory script for RNA-seq experiments performed for the Van Doorslaer lab at the University of Arizona. Updates to script will be uploaded when available.

Key items that need to be updated when running a new analysis:

#Import Count Data

setwd() - Simple, but don't forget!
x_counts - Specify .csv, .tsv, or other file. 
cts - Use read.delim() function for non .csv.

#DESeq parameters

x_meta - file.path() for metadata file used by DESeq2.

#Shrink Log-fold change

resultsNames(dds) - Can be used to fetch coef for lfcshrink().

#Genes of Interest Annotation

topVarGenes - Change # based on scope of analysis.
mart - Change dataset and ID if different species.
gns - Change symbols, IDs used if necessary.

#Make Heatmap

anno - Add proper conditions for analysis.

#Volcanoplot

Filtering limits can be changed if necessary.
geom_text_repel - Use data subset highlighted for significance.

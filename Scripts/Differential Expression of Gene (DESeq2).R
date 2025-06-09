# Library load

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(GEOquery)
library(tidyverse)
library(org.Hs.eg.db)

# read counts data
counts_data <- read.csv("Data/GSE161969_raw_counts.tsv", sep = '\t', row.names = 1)
head(counts_data)     # first few values of counts data

# read the sample information. This is also called metadata. metadata is available along with raw data.
# here I created a tsv file named design which has two column "Sample" and "Group" and did not obtain metadata
sample_info <- read.csv("Data/design.csv", row.names = 1)
sample_info
dim(sample_info)  # dimension

# set factor level
factors <- factor(sample_info$Group)
groups <- unique(sample_info$Group)

groups
groups <- rev(groups)

groups
sample_info$Group <- factors
sample_info$Group

# create DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = sample_info, design = ~Group)

# set the reference for the group factor
dds$Group <- relevel(dds$Group, ref = "control")

# filter out low count genes
# keep genes with at least N counts >= 10, where N = size of smallest group
keep <- rowSums(counts(dds) >= 10) >= min(table(sample_info$Group))
dds <- dds[keep,]   


# perform statistical tests
dds <- DESeq(dds, test = "Wald", sfType='poscount')

# get the deseq result
deseq_result <- results(dds)

# convert result into data frame
deseq_result <- as.data.frame(deseq_result)
head(deseq_result)
dim(deseq_result)
names(deseq_result)   # column names of deseq result


# deseq_result table has row names that are gene ID(ENTREZID)
# So I created a new column with those geneID
deseq_result$GeneName <- row.names(deseq_result)
names(deseq_result)
head(deseq_result)

# reorganize table content
deseq_result <- subset(deseq_result,
                       select = c("GeneName","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean"))
names(deseq_result)


# This RNA-seq data used ENTREZ ID, but it is hard to understand
# genes by that ID, so I will convert those ID to Gene names
# Map Entrez IDs to gene symbols
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = as.character(deseq_result$GeneName),
  columns = c("SYMBOL"),
  keytype = "ENTREZID"
)

# Merge the gene symbols into the count matrix
deseq_annotate <- deseq_result %>%
  left_join(symbol_map, by = c("GeneName" = "ENTREZID"))

# Optional: move SYMBOL column to first
deseq_annotate <- deseq_annotate %>%
  relocate(SYMBOL)

# Optional: rename column
colnames(deseq_annotate)[1] <- "Gene"



# remove column with EntrezID
deseq_annotate$GeneName <- NULL

# export deseq result. I can use deseq_annotate variable name, but for the sake of simplicity I am renaming it
deseq_result <- deseq_annotate
write.table(deseq_result, file = "Data/deseq_result.all.tsv", row.names=F, sep = "\t")


# extract de genes with padj < 0.05 and log2FoldChange <= -1 or >= 1. log2FC <= -2 or >= 2 is also used, and it will give less deg or more significant deg

deg <- subset(deseq_result, padj<0.05 & abs(log2FoldChange) >= 1)
dim(deg) # means dimension of deg
deg <- deg[order(deg$padj),]
head(deg)


# finally export deg table
write.table(deg, file = "Data/deseq_deg.tsv", row.names = F, sep = "\t")






















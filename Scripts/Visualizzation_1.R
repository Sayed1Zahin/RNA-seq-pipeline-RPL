library(pheatmap)
library(RColorBrewer)
library(tidyverse)


# plot dispersion estimate

plotDispEsts(dds, main = "GSE161969 Dispersion Estimate")

# create a histogram plot of p-values
hist(deseq_result$padj, breaks = seq(0,1,length = 21), col = "gray", border = "white",
     xlab = "", ylab = "", ylim = c(0,3500), main = " GSE161969 Frequencies of padg-values")


# volcano plot
# set colors
old.pal <- palette(c("#00BFFF", "#FF3030"))

# set margin size
par(mar=c(4,4,2,1), cex.main = 1.5)

# set title
title = paste(groups[1], "vs", groups[2])
# plot values

plot(deseq_result$log2FoldChange, -log10(deseq_result$padj), main = title,
     xlab = "log2FC", ylab = "-log10(Padj)", pch = 20, cex = 0.5)


with(subset(deseq_result, padj < 0.05 & abs(log2FoldChange) >= 1),
     points(log2FoldChange, -log10(padj), pch = 20, col = (sign(log2FoldChange) +3)/2, cex = 1))
     

legend("bottomleft", title = paste("Padj<", 0.05, sep = ""),
       legend = c("down","up"), pch = 20, col = 1:2)


# variance stabilizing transformation'

vsd <- vst(dds, blind = FALSE)

# pca plot
# use transformed values to generate a pca plot
plotPCA(vsd, intgroup = c("Group"))


# heatmap

normalized_counts <- counts(dds, normalized = T)
head(normalized_counts)

transformed_counts <- log2(normalized_counts+1)
head(transformed_counts)

# extract de genes with padj < 0.05 and log2FoldChange <= -1 or >= 1. log2FC <= -2 or >= 2 is also used, and it will give less deg or more significant deg

deg <- subset(deseq_result, padj<0.05 & abs(log2FoldChange) >= 1)
dim(deg) # means dimension of deg
deg <- deg[order(deg$padj),]
head(deg)


top_hits <- deg$GeneName[1:10]
top_hits <- transformed_counts[top_hits,] 


# top_hits contain 10 genes where row name is ENTREZID
# create a new column with those geneID
# in order to do that, i will need to make top_hits as data frame

top_hits <- as.data.frame(top_hits)

# create a new column
top_hits$GeneName <- row.names(top_hits)
names(top_hits)

# Map Entrez IDs to gene symbols
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = as.character(top_hits$GeneName),
  columns = c("SYMBOL"),
  keytype = "ENTREZID"
)

# Merge the gene symbols into the count matrix
top_hits <- top_hits %>%
  left_join(symbol_map, by = c("GeneName" = "ENTREZID"))

# Define your gene symbols
gene_symbols <- c("HLA-J", "KCND2", "MX1", "LAMA2", "MMP8", 
                  "ISG15", "LFNG", "IFI44L", "IFIT1", "TNFRSF21")

# 'top_hits' is existing dataset with 10 rows
# Assign row names
rownames(top_hits) <- gene_symbols

top_hits
row.names <- top_hits$SYMBOL

# now exclude GeneName and SYMBOL column
top_hits <- top_hits[, 1:7]


# all these steps are needed just because our top_hits list had ENTREZID instead of Gene Name(SYMBOL)


# finally heat map
pheatmap(top_hits,cluster_rows = FALSE, cluster_cols = FALSE)










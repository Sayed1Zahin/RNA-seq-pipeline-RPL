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


with(subset(deseq_result, padj < 0.05 & abs(log2FoldChange) >= 2),
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

top <- deg$Gene
top[1:10]


top_hits <- deg$Gene[1:10]
top_hits <- transformed_counts[top_hits,] 










library(ggplot2)
library(ggrepel)
library(dplyr)


# load data
deseq_deg <- read.csv("Data/deseq_result.all.tsv", sep = '\t')
head(deseq_deg)
colnames(deseq_deg)

# create a column
deseq_deg$diffexpresssed <- "NO"

# set up condition for up and downregulated gene
deseq_deg$diffexpresssed[deseq_deg$log2FoldChange>0.58 & deseq_deg$padj<0.01] <- "UP"
deseq_deg$diffexpresssed[deseq_deg$log2FoldChange<0.58 & deseq_deg$padj<0.01] <- "DOWN"
head(deseq_deg)

# create a new column
deseq_deg$delabel <- NA

# creating volcano plot
ggplot(data = deseq_deg, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpresssed, label = delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values = c('blue','black','red'))+
  geom_vline(xintercept = c(-0.8,0.8), col = 'red')+
  geom_hline(yintercept = -log10(0.0001), col = 'red')+
  theme(text = element_text(size = 20))


# sort top genes
# Sort by p-value
deseq_deg <- deseq_deg[order(deseq_deg$pvalue), ]
# Take top 20 with non-NA gene names
top_genes <- head(deseq_deg[!is.na(deseq_deg$Gene), ], 20)

# Label only top genes
deseq_deg$delabel[deseq_deg$Gene %in% top_genes$Gene] <- top_genes$Gene



ggplot(data = deseq_deg, aes(x = log2FoldChange, y = -log10(pvalue), 
                             col = diffexpresssed, label = delabel)) +
  geom_point() +
  geom_text_repel(na.rm = TRUE, max.overlaps = 100) + # removes NA labels safely
  scale_color_manual(values = c('blue', 'black', 'red')) +
  geom_vline(xintercept = c(-0.8, 0.8), col = 'red') +
  geom_hline(yintercept = -log10(0.0001), col = 'red') +
  theme_minimal() +
  theme(text = element_text(size = 20))









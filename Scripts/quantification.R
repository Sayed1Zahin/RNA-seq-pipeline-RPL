
BiocManager::install("org.Hs.eg.db")
# load required packages
library(tidyverse)
library(GEOquery)
library(ggpubr)
library(openxlsx)
library(naniar)

# import data
counts_data <- read.csv("Data/GSE161969_fpkm.csv")
res <- getGEO(GEO = "GSE161969", GSEMatrix = T)
class(res)


# metadata
metadata <- pData(phenoData(res[[1]]))


# subset metadata
metadata %>% 
  head()
metadata_sub <- metadata %>% 
  select(c(2,8,10))

# data preprocessing
metadata_modified <- metadata_sub %>% 
  rename(tissue = source_name_ch1, character = characteristics_ch1) %>% 
  mutate(tissue = gsub("decidual ", "", tissue)) %>% 
  mutate(character = gsub("subject status: ", "", character))          


# Load necessary libraries
library(org.Hs.eg.db)
library(dplyr)  # for data manipulation

# Assume your count matrix is called `count_data`
# and Entrez IDs are the rownames
# Example:
# rownames(count_data) <- c("653635", "100996442", "729737")

# Convert rownames to a column so we can merge
counts_data$GeneID <- rownames(counts_data)

# Map Entrez IDs to gene symbols
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = as.character(counts_data$GeneID),
  columns = c("SYMBOL"),
  keytype = "ENTREZID"
)

# Merge the gene symbols into the count matrix
count_data_annotated <- counts_data %>%
  left_join(symbol_map, by = c("GeneID" = "ENTREZID"))

# Optional: move SYMBOL column to first
count_data_annotated <- count_data_annotated %>%
  relocate(SYMBOL)

# Optional: rename column
colnames(count_data_annotated)[1] <- "Gene"

# Set rownames again (optional)
rownames(count_data_annotated) <- count_data_annotated$EntrezID



count_data_annotated$GeneID <- NULL






#reshaping data
counts_data_long <- count_data_annotated %>% 
  pivot_longer(-Gene,
               names_to = "samples",
               values_to = "fpkm")


counts_final <- counts_data_long %>% 
  left_join(metadata_modified, by = c("samples" = "geo_accession"))

write.csv(counts_final, "Data/GSE161969_counts.csv", row.names = F)

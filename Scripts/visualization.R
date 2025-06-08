library(tidyverse)
data <- read.csv("Data/GSE161969_counts.csv")


# Convert count data to matrix if needed
counts_matrix <- as.matrix(data)

# Example coldata (sample metadata)
# Must match sample column names in counts_matrix
coldata <- data.frame(
  row.names = colnames(counts_matrix),
  condition = c("control", "recurrent pregnancy loss (RPL)")
    # example conditions
)

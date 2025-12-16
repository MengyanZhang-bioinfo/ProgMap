#
library(UpSetR)
library(dplyr)
library(tidyr)
pan_gene=read.table("./CPSig.txt",header=TRUE)
gene_matrix <- pan_gene %>%
  distinct(cancer, gene) %>%
  mutate(value = 1) %>%
  spread(key = cancer, value = value, fill = 0)

gene_matrix <- as.data.frame(gene_matrix)

upset(gene_matrix, 
      nintersects=39,
      sets = colnames(gene_matrix)[-1],
      sets.bar.color = "lightyellow",       
      order.by = "freq",                
      matrix.color = "black",         
      keep.order = TRUE, 
      empty.intersections=TRUE)                

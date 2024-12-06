### Filter out edges with weight of zero ### 

library(tidyverse)
network_AD <- read.table("network_AD_col_all.csv.gz", header = T, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
head(network_AD)
dim(network_AD)
network_AD_ <- network_AD %>% dplyr::filter(edge != 0)
dim(network_AD_)
write.table(network_AD_, file = "network_AD_edges.csv", sep ='\t', col.names = T, row.names = F)
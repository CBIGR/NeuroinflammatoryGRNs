### Filter out edges with weight of zero ### 

library(tidyverse, lib.loc = '/scratch/gent/vo/000/gvo00027/vsc44365/R_packages')
network_MDD <- read.table("network_MDD_col_all.csv.gz", header = T, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
head(network_MDD)
dim(network_MDD)
network_MDD_ <- network_MDD %>% dplyr::filter(edge != 0)
dim(network_MDD_)
write.table(network_MDD_, file = "network_MDD_edges.csv", sep ='\t', col.names = T, row.names = F)
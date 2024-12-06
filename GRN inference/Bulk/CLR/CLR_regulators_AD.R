### Script to select regulatory edges and top edges ###

library(tidyverse, lib.loc = '/scratch/gent/vo/000/gvo00027/vsc44365/R_packages')

# read in files
network_AD <- read.table("network_AD_edges.csv.gz",  header = T, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
head(network_AD)
regulators <- read.table('RegulatorsAD.txt')
regulators <- as.matrix(regulators)
regulators <- as.vector(regulators)

#select rows with regulators 
#in first step every row with TF is selected once, as well as those with 2TFs (both rows)
network_AD_ <- network_AD %>% dplyr::filter(gene1 %in% regulators)
dim(network_AD_)

#arrange according to edge weight and select top
network_AD_ <- network_AD_ %>% arrange(desc(edge))
head(network_AD_)
tail(network_AD_)
network_AD_ <- network_AD_ %>% slice(c(1:200000))
dim(network_AD_)

write.table(network_AD_, file = "network_AD_reg.csv", sep ='\t', col.names = T, row.names = F)
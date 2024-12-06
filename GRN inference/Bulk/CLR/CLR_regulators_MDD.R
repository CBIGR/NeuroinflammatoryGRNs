### Script to select regulatory edges and top edges ###

library(tidyverse, lib.loc = '/scratch/gent/vo/000/gvo00027/vsc44365/R_packages')

# read in files
network_MDD <- read.table("network_MDD_edges.csv.gz",  header = T, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
head(network_MDD)
regulators <- read.table('RegulatorsMDD.txt')
regulators <- as.matrix(regulators)
regulators <- as.vector(regulators)

#select rows with regulators
network_MDD_ <- network_MDD %>% dplyr::filter(gene1 %in% regulators)
dim(network_MDD_)

#arrange according to edge weight and select top
network_MDD_ <- network_MDD_ %>% arrange(desc(edge))
head(network_MDD_)
tail(network_MDD_)
network_MDD_ <- network_MDD_ %>% slice(c(1:200000))
dim(network_MDD_)

write.table(network_MDD_, file = "network_MDD_reg.csv", sep ='\t', col.names = T, row.names = F)
### Network inference with CLR using package MINET ###

# load MINET
library(minet, lib.loc = '/scratch/gent/vo/000/gvo00027/vsc44365/R_packages')

dataset_MDD <- read.csv("Processed_counts_MDD.csv", header = T, sep = '\t', row.names = NULL)
GeneSymbol <- dataset_MDD$GeneSymbol
dataset_MDD_m <- t(dataset_MDD[, -1])
dataset_MDD <- as.data.frame(dataset_MDD_m)
colnames(dataset_MDD) <- GeneSymbol

# Build mutual information matrix
mim_MDD <- build.mim(dataset_MDD, estimator = "mi.mm", disc = "equalfreq") # dataset must be dataframe
#Build network
network_MDD <- clr(mim_MDD) # weighted adjacency matrix of network
head(network_MDD)
dim(network_MDD)

# Make an according edge list (datframe) out of the network matrix (from Joke)
rownames_1 <- rep(row.names(network_MDD), length(colnames(network_MDD))) 
rownames_2 <- rep(colnames(network_MDD), each = length(row.names(network_MDD)))
network_MDD_col_all <- data.frame(rownames_1, rownames_2, c(network_MDD))
colnames(network_MDD_col_all) <- c('gene1', 'gene2', 'edge')
head(network_MDD_col_all)
dim(network_MDD_col_all)

# Write file
write.table(network_MDD, file = "network_MDD.csv", sep = '\t', row.names = T, col.names = T)
write.table(network_MDD_col_all, file = "network_MDD_col_all.csv", sep = '\t', row.names = F, col.names = T)
#after network inference, select regulatory links of known regulators
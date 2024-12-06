### Network inference with CLR using package MINET ###


# load MINET
library(minet, lib.loc = '/scratch/gent/vo/000/gvo00027/vsc44365/R_packages')

dataset_AD <- read.csv("Processed_counts_AD.csv", header = T, sep = '\t', row.names = NULL)
GeneSymbol <- dataset_AD$GeneSymbol
dataset_AD_m <- t(dataset_AD[, -1])
dataset_AD <- as.data.frame(dataset_AD_m)
colnames(dataset_AD) <- GeneSymbol

# Build mutual information matrix
mim_AD <- build.mim(dataset_AD, estimator = "mi.mm", disc = "equalfreq") #dataset must be dataframe
#Build network
network_AD <- clr(mim_AD) #weighted adjacency matrix of network
head(network_AD)
dim(network_AD)

# Make an according edge list (datframe) out of the network matrix (from Joke)
rownames_1 <- rep(row.names(network_AD), length(colnames(network_AD))) 
rownames_2 <- rep(colnames(network_AD), each = length(row.names(network_AD)))
network_AD_col_all <- data.frame(rownames_1, rownames_2, c(network_AD))
colnames(network_AD_col_all) <- c('gene1', 'gene2', 'edge')
head(network_AD_col_all)
dim(network_AD_col_all)

# Write file
write.table(network_AD, file = "network_AD.csv", sep = '\t', row.names = T, col.names = T)
write.table(network_AD_col_all, file = "network_AD_col_all.csv", sep = '\t', row.names = F, col.names = T)
# after network inference, select regulatory links of known regulators
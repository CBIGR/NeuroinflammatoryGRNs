### Network inference with GENIE3 ###

# load GENIE3
library(GENIE3)

# Make weight matrix
# exprMatr needs to be matrix
exprMatr_AD <- read.table('Processed_counts_AD.csv', header = T, sep = '\t')
GeneSymbol <- exprMatr_AD$GeneSymbol
exprMatr_AD <- as.matrix(exprMatr_AD[, -1])
rownames(exprMatr_AD) <- GeneSymbol
regulators <- read.table('RegulatorsAD.txt') # regulator list, must be vector
regulators <- as.matrix(regulators)
regulators <- as.vector(regulators)
weightMat_AD <- GENIE3(exprMatr_AD, nCores = 16, regulators = regulators)

# List of regulatory links 
linkList_AD <- getLinkList(weightMat_AD)
dim(linkList_AD)
linkList_AD_top <- getLinkList(weightMat_AD, reportMax = 200000)
head(linkList_AD) # weights do not have any statistical meaning

# Save weight matrix and list in file
write.table(weightMat_AD, file = 'weightMat_AD.csv', sep = '\t', row.names = T, col.names = T)
write.table(linkList_AD, file = 'LinkList_AD.csv', sep = '\t', row.names = F, col.names = T)
write.table(linkList_AD_top, file = 'LinkList_AD_top.csv', sep = '\t', row.names = F, col.names = T)

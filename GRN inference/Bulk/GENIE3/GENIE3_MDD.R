### Network inference with GENIE3 ###

#load GENIE3
library(GENIE3)

# Make weight matrix
# exprMatr needs to be matrix
exprMatr_MDD <- read.table('Processed_counts_MDD.csv', header = T, sep = '\t')
GeneSymbol <- exprMatr_MDD$GeneSymbol
exprMatr_MDD <- as.matrix(exprMatr_MDD[, -1])
rownames(exprMatr_MDD) <- GeneSymbol
regulators <- read.table('RegulatorsMDD.txt') # regulator list, must be vector
regulators <- as.matrix(regulators)
regulators <- as.vector(regulators)
weightMat_MDD <- GENIE3(exprMatr_MDD, nCores = 16, regulators = regulators)

# List of regulatory links 
linkList_MDD <- getLinkList(weightMat_MDD)
dim(linkList_MDD)
linkList_MDD_top <- getLinkList(weightMat_MDD, reportMax = 200000)
head(linkList_MDD) # weights do not have any statistical meaning

# Save weight matrix and list in file
write.table(weightMat_MDD, file = 'weightMat_MDD.csv', sep = '\t', row.names = T, col.names = T)
write.table(linkList_MDD, file = 'LinkList_MDD.csv', sep = '\t', row.names = F, col.names = T)
write.table(linkList_MDD_top, file = 'LinkList_MDD_top.csv', sep = '\t', row.names = F, col.names = T)

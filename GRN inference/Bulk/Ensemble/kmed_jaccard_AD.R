### Clustering with k-medoids on jaccard index with cluster package ### 


library(reshape2)
library(cluster)

Jac <- read.table('geneclust_out_AD', header = F, sep='\t', quote="", col.names= c("g1","g2","j"))
head(Jac)
# gene1 \t gene2 \t jaccard index

Jplus <- data.frame(Jac$g2, Jac$g1, Jac$j)
colnames(Jplus) <- colnames(Jac)
JJ <- rbind(Jac,Jplus)
JJ_T <- acast(JJ, g1~g2, value.var="j")
TT <- matrix(1, nrow = dim(JJ_T)[1], ncol = dim(JJ_T)[2])
diag(JJ_T) <- 1
JJ_T[is.na(JJ_T)] <- 0 #missing values to zero
Tf <- TT-JJ_T
rm(Jplus, JJ, JJ_T, TT)

# k-medoids
clustersAD <- pam(Tf, k = 155, cluster.only = T)
write.table(clustersAD, file = "kmed_modules_AD.txt", col.names = F, row.names = T, sep = "\t", quote = F)


### Clustering with k-medoids on jaccard index with cluster package ### 

library(reshape2, lib.loc = '/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/Rlib')
library(cluster, lib.loc = '/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/Rlib')

Jac <- read.table('geneclust_out_MDD', header = F, sep='\t', quote="", col.names= c("g1","g2","j"))
#head(Jac)
# gene1 \t gene2 \t jaccard index

Jplus <- data.frame(Jac$g2, Jac$g1, Jac$j)
colnames(Jplus) <- colnames(Jac)
JJ <- rbind(Jac,Jplus)
JJ_T <- acast(JJ, g1~g2, value.var="j")
TT <- matrix(1, nrow = dim(JJ_T)[1], ncol = dim(JJ_T)[2])
diag(JJ_T) <- 1
JJ_T[is.na(JJ_T)] <- 0 #missing values to zero
Tf <- TT-JJ_T
#rm(TT)

# k-medoids
clustersMDD <- pam(Tf, k = 156, cluster.only = T)
write.table(clustersMDD, file = "kmed_modules_MDD.txt", col.names = F, row.names = T, sep = "\t", quote = F)


### determine number of clusters to use in k-medoids ###

library(NbClust, lib.loc = '/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/Rlib')
library(tidyverse, lib.loc = '/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/Rlib')

countsAD <- read.table("/scratch/gent/vo/000/gvo00027/projects/CBIGR/21HPU_GRN_neuroinfl/Processed_counts_AD_GS.txt.gz", sep = "\t", header = T)
countsMDD <- read.table("/scratch/gent/vo/000/gvo00027/projects/CBIGR/21HPU_GRN_neuroinfl/Processed_counts_MDD_GS.txt.gz", sep = "\t", header = T)
countsAD <- column_to_rownames(countsAD, "GeneSymbol")
countsMDD <- column_to_rownames(countsMDD, "GeneSymbol")

print("start with AD")
nb_clust <- NbClust(data = countsAD, min.nc = 150, max.nc = 250,
                    method = "median", index = "ch")
nb_clust2 <- NbClust(data = countsAD, min.nc = 150, max.nc = 250,
                    method = "kmeans", index = "ch")
nb_clust3 <- NbClust(data = countsAD, min.nc = 150, max.nc = 250,
                     method = "median", index = "silhouette")
nb_clust4 <- NbClust(data = countsAD, min.nc = 150, max.nc = 250,
                    method = "kmeans", index = "silhouette")

ch <- nb_clust$Best.nc
ch2 <- nb_clust2$Best.nc
silhouette <- nb_clust3$Best.nc
silhouette2 <- nb_clust4$Best.nc
df <- data.frame(ch, ch2, silhouette, silhouette2)
print(df)
write.table(df, file = "Nb_clust_AD.txt", row.names = T, col.names = T, sep = "\t")

print("start with MDD")
nb_clust_MDD <- NbClust(data = countsMDD, min.nc = 150, max.nc = 250,
                    method = "median", index = "ch")
nb_clust2_MDD <- NbClust(data = countsMDD, min.nc = 150, max.nc = 250,
                        method = "kmeans", index = "ch")
nb_clust3_MDD <- NbClust(data = countsMDD, min.nc = 150, max.nc = 250,
                     method = "median", index = "silhouette")
nb_clust4_MDD <- NbClust(data = countsMDD, min.nc = 150, max.nc = 250,
                         method = "kmeans", index = "silhouette")

ch <- nb_clust_MDD$Best.nc
ch2 <- nb_clust2_MDD$Best.nc
silhouette <- nb_clust3_MDD$Best.nc
silhouette2 <- nb_clust4_MDD$Best.nc
df_MDD <- data.frame(ch, ch2, silhouette, silhouette2)
print(df_MDD)
write.table(df_MDD, file = "Nb_clust_MDD.txt", row.names = T, col.names = T, sep = "\t")
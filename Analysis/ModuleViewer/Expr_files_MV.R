### Script to get expression files in right format for ModuleViewer ###

Processed_counts_AD_GS <- read.delim("C:/UGent/Master 2/Master's dissertation/R/Processed_counts_AD_GS.txt")
head(Processed_counts_AD_GS)
colnames(Processed_counts_AD_GS) <- str_replace(colnames(Processed_counts_AD_GS),"\\.", "-")
write.table(Processed_counts_AD_GS, file = "Processed_counts_AD_GS_MV.txt", quote = F, sep = "\t")

Processed_counts_MDD_GS <- read.delim("C:/UGent/Master 2/Master's dissertation/R/Processed_counts_MDD_GS.txt")
head(Processed_counts_MDD_GS)
write.table(Processed_counts_MDD_GS, file = "Processed_counts_MDD_GS_MV.txt", quote = F, sep = "\t")

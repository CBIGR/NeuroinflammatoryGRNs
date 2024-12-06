### Script to put all ChIP-based enrichR target gene enrichment in the excel file Overview modules
# input file: mvf file without headers
# author: Hanne Puype
# date: 17/08/2023

library(tidyverse)

setwd("C:/Users/hannepu/OneDrive - UGent/Master 2/Master's dissertation/Analysis")

#### AD ####
ChIP_AD <- read.table("Chip_enrichr.txt", sep = "\t", header = F)
colnames(ChIP_AD) <- c("Module", "TGs", "TF")

ChIP_AD <- ChIP_AD %>% select(-TGs)
ChIP_AD <- unique(ChIP_AD)


ChIP_AD_wider <- ChIP_AD %>% pivot_wider(names_from = Module, values_from = TF)
ChIP_AD_wider_t <- t(ChIP_AD_wider)

ChIP_AD_wider_t <- rownames_to_column(as.data.frame(ChIP_AD_wider_t), "Module")
ChIP_AD_wider_t$Module <- as.numeric(ChIP_AD_wider_t$Module)

for (i in (1:155)) {
  if (!(i %in% ChIP_AD_wider_t$Module)) {
  ChIP_AD_wider_t[dim(ChIP_AD_wider_t)[1]+1,] <- c(i, "")
  }
} 

ChIP_AD_wider_t$Module <- as.numeric(ChIP_AD_wider_t$Module)
ChIP_AD_wider_t <- ChIP_AD_wider_t %>% arrange(Module)

ChIP_AD_wider_t$V1 <- as.character(ChIP_AD_wider_t$V1)
write.table(ChIP_AD_wider_t, file = "ChIP_all_AD.txt", row.names = F, col.names = F, quote = F, sep = "\t")


#### MDD ####
ChIP_MDD <- read.table("Chip_enrichr_MDD.txt", sep = "\t", header = F)
colnames(ChIP_MDD) <- c("Module", "TGs", "TF")

ChIP_MDD <- ChIP_MDD %>% select(-TGs)
ChIP_MDD <- unique(ChIP_MDD)

ChIP_MDD_wider <- ChIP_MDD %>% pivot_wider(names_from = Module, values_from = TF)
ChIP_MDD_wider_t <- t(ChIP_MDD_wider)

ChIP_MDD_wider_t <- rownames_to_column(as.data.frame(ChIP_MDD_wider_t), "Module")
ChIP_MDD_wider_t$Module <- as.numeric(ChIP_MDD_wider_t$Module)

for (i in (1:156)) {
  if (!(i %in% ChIP_MDD_wider_t$Module)) {
    ChIP_MDD_wider_t[dim(ChIP_MDD_wider_t)[1]+1,] <- c(i, "")
  }
} 

ChIP_MDD_wider_t$Module <- as.numeric(ChIP_MDD_wider_t$Module)
ChIP_MDD_wider_t <- ChIP_MDD_wider_t %>% arrange(Module)

ChIP_MDD_wider_t$V1 <- as.character(ChIP_MDD_wider_t$V1)
write.table(ChIP_MDD_wider_t, file = "ChIP_all_MDD.txt", row.names = F, col.names = F, quote = F, sep = "\t")

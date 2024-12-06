#### Custom gene set modules for GWAS enrichment ####

library(tidyverse)
kmed_modules_MV_AD <- read.delim("C:/Users/hannepu/OneDrive - UGent/Master 2/Master's dissertation/R/ModuleViewer/Ensemble/kmed_modules_MV_AD.txt", 
                                 header = FALSE)
head(kmed_modules_MV_AD)
kmed_modules_MV_AD[c('V')] <- str_split_fixed(kmed_modules_MV_AD$V2, pattern = '\\|', n = 160)
head(kmed_modules_MV_AD)
kmed_modules_MV_AD <- kmed_modules_MV_AD %>% dplyr::select(-("V2"))
kmed_modules_MV_AD['V2'] <- "description"
kmed_modules_MV_AD <- kmed_modules_MV_AD %>% dplyr::select(V1, V2, everything())

V1 <- kmed_modules_MV_AD$V1
V1 <- lapply(V1, trimws) #is vector
V1 <- unlist(V1)
kmed_modules_MV_AD$V1 <- V1

kmed_modules_MV_AD$V1 <- paste("Module", kmed_modules_MV_AD$V1, sep="_")

write.table(kmed_modules_MV_AD, file = "C:/PhD_22-23/Analyses/R/Modules_tab_AD.gmt", row.names = F, col.names = F, 
            quote = F, sep = "\t")

### MDD ###
kmed_modules_MV_MDD <- read.delim("C:/Users/hannepu/OneDrive - UGent/Master 2/Master's dissertation/R/ModuleViewer/Ensemble/kmed_modules_MV_MDD.txt", 
                                 header=FALSE)
head(kmed_modules_MV_MDD)

kmed_modules_MV_MDD[c('V')] <- str_split_fixed(kmed_modules_MV_MDD$V2, pattern = '\\|', n = 160)
head(kmed_modules_MV_MDD)
kmed_modules_MV_MDD <- kmed_modules_MV_MDD %>% dplyr::select(-("V2"))
kmed_modules_MV_MDD['V2'] <- "description"
kmed_modules_MV_MDD <- kmed_modules_MV_MDD %>% dplyr::select(V1, V2, everything())

V1 <- kmed_modules_MV_MDD$V1
V1 <- lapply(V1, trimws) #is vector
V1 <- unlist(V1)
kmed_modules_MV_MDD$V1 <- V1
write.table(kmed_modules_MV_MDD, file = "C:/PhD_22-23/Analyses/R/Modules_tab_MDD.gmt", row.names = F, 
            col.names = F, quote = F, sep = "\t")


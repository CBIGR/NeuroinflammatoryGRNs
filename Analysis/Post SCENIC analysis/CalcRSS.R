### RSS calculation ###
library(Seurat)
library(SCENIC)
library(tidyverse) 
library(SCopeLoomR)
library(BBmisc)

# AD
seuratObjAD <- readRDS("~/MyFiles/Master/single-cell/Seurat_AD.rds")
seuratObjAD <- SetIdent(seuratObjAD, value = "celltype")
cellAnn <- Idents(seuratObjAD)
loom <- open_loom(file.path = "~/MyFiles/Master/single-cell/SCENIC_output.loom")
regulonAUC <- get_regulons_AUC(loom)
close_loom(loom)
rss <- calcRSS(AUC = regulonAUC, cellAnnotation = cellAnn, cellTypes = NULL) 
RSS_AD <- as.data.frame(rss)
write.table(RSS_AD, file = "~/MyFiles/Master/single-cell/RSS_AD.csv", sep = ",", row.names = T, col.names = T, quote = F)

# MDD
seuratObjMDD <- readRDS("~/MyFiles/Master/single-cell/Seurat_MDD.rds")
dim(seuratObjMDD@assays$RNA)
loom <- open_loom(file.path = "~/MyFiles/Master/single-cell/SCENIC_output_MDD.loom")
regulonAUC <- get_regulons_AUC(loom)
dim(regulonAUC@assays@data@listData[["AUC"]])
close_loom(loom)
seuratObjMDD@assays[["RNA"]]@data@Dimnames[[2]] <- sub("/", ".", seuratObjMDD@assays[["RNA"]]@data@Dimnames[[2]]) # micro.macro is different
seuratObjMDD <- SetIdent(seuratObjMDD, value = "Celltype")
cellAnn <- Idents(seuratObjMDD)
rss <- calcRSS(AUC = regulonAUC, cellAnnotation = cellAnn, cellTypes = NULL) 
RSS_MDD <- as.data.frame(rss)
min(RSS_MDD)
write.table(RSS_MDD, file = "~/MyFiles/Master/single-cell/RSS_MDD.csv", quote = F, sep = ",", row.names = T, col.names = T)

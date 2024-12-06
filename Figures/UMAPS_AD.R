###############################################################################
# CREATE PLOTS FROM pySCENIC OUTPUT
###############################################################################

library("SCENIC") 
library("BBmisc")
#install.packages("remotes")
#install.packages("hdf5r")
#remotes::install_github("aertslab/SCopeLoomR")
library("SCopeLoomR")
library("Seurat")
library("dplyr")
library("ggplot2")
library(cowplot)
library(tidyverse)

path_to_loom <- "~/MyFiles/Master/" 

set.seed(123)


# Read data from loom file
loom <- open_loom(file.path = paste0(path_to_loom,"SCENIC_output.loom"))
regulonAUC <- get_regulons_AUC(loom)
close_loom(loom)

# Get AUC values in matrix + scale to range 0-1 for nicer heatmaps
AUCdata <- getAUC(regulonAUC) %>% BBmisc::normalize("range")

## EXPRESSION UMAP ##
seuratObj <- readRDS("~/MyFiles/Master/Seurat_AD.rds")
seuratObj_umapdata <- as.data.frame(seuratObj@reductions$umap@cell.embeddings) #retrieve UMAP embedding from the expression data 
seuratObj_umapdata$clusters<- seuratObj@meta.data$celltype #add extra column next to UMAP_1 and UMAP_2 which you want to use as colored clusters in your plot 
head(seuratObj_umapdata)
seuratObj_umapdata %>% group_by(clusters) %>% summarise(n())
p <- ggplot(data = seuratObj_umapdata, aes(x = UMAP_1, y = UMAP_2, color = clusters )) + 
  geom_point(shape = ".") +
  labs(color='Cell type', title = "UMAP of expression data") +
  theme_classic()
p  

p2 <- ggplot(data = seuratObj_umapdata, aes(x = UMAP_1, y = UMAP_2, color = clusters )) + 
  geom_point() +
  labs(color='Cell type', title = "UMAP of expression data") +
  theme_classic()
p2

plot_grid(p, p2, ncol = 2)

# Plot UMAP according to metadata
snRNA_metadta <- read_csv("~/MyFiles/Master/snRNA_metadta.csv")
snRNA_metadta <- snRNA_metadta %>% filter(...1 %in% colnames(seuratObj))
seuratObj[["Sex"]] <- snRNA_metadta$Sex
seuratObj[["Age"]] <- snRNA_metadta$Age
seuratObj_umapdata$batch<- seuratObj@meta.data$Batch
seuratObj_umapdata$batch <- as.factor(seuratObj_umapdata$batch)
seuratObj_umapdata$diagnosis<- seuratObj@meta.data$Diagnosis
seuratObj_umapdata$sample<- seuratObj@meta.data$Sample
seuratObj_umapdata$sex<- seuratObj@meta.data$Sex
seuratObj_umapdata$age<- seuratObj@meta.data$Age

p <- ggplot(data = seuratObj_umapdata, aes(x = UMAP_1, y = UMAP_2, color = batch )) + 
  geom_point(shape = ".") +
  labs(color='Batch', title = "UMAP of expression according to batch") +
  theme_classic()
p 
p2 <- ggplot(data = seuratObj_umapdata, aes(x = UMAP_1, y = UMAP_2, color = batch )) + 
  geom_point() +
  labs(color='Batch', title = "UMAP of expression according to batch") +
  theme_classic()
p2
plot_grid(p, p2, ncol = 2)

p <- ggplot(data = seuratObj_umapdata, aes(x = UMAP_1, y = UMAP_2, color = diagnosis )) + 
  geom_point(shape = ".") +
  labs(color='Diagnosis', title = "UMAP of expression according to diagnosis") +
  theme_classic()
p 
p2 <- ggplot(data = seuratObj_umapdata, aes(x = UMAP_1, y = UMAP_2, color = diagnosis )) + 
  geom_point() +
  labs(color='Diagnosis', title = "UMAP of expression according to diagnosis") +
  theme_classic()
p2
plot_grid(p, p2, ncol = 2)

p <- ggplot(data = seuratObj_umapdata, aes(x = UMAP_1, y = UMAP_2, color = sample )) + 
  geom_point(shape = ".") +
  labs(color='Sample', title = "UMAP of expression according to sample") +
  theme_classic()
p 
p2 <- ggplot(data = seuratObj_umapdata, aes(x = UMAP_1, y = UMAP_2, color = sample )) + 
  geom_point() +
  labs(color='Sample', title = "UMAP of expression according to sample") +
  theme_classic()
p2
plot_grid(p, p2, ncol = 2)

p <- ggplot(data = seuratObj_umapdata, aes(x = UMAP_1, y = UMAP_2, color = sex)) + 
  geom_point(shape = ".") +
  labs(color='Sex', title = "UMAP of expression according to sex") +
  theme_classic()
p 
p2 <- ggplot(data = seuratObj_umapdata, aes(x = UMAP_1, y = UMAP_2, color = sex)) + 
  geom_point() +
  labs(color='Sex', title = "UMAP of expression according to sex") +
  theme_classic()
p2
plot_grid(p, p2, ncol = 2)

p <- ggplot(data = seuratObj_umapdata, aes(x = UMAP_1, y = UMAP_2, color = age )) + 
  geom_point(shape = ".") +
  labs(color='Age', title = "UMAP of expression according to age") +
  theme_classic()
p 
p2 <- ggplot(data = seuratObj_umapdata, aes(x = UMAP_1, y = UMAP_2, color = age )) + 
  geom_point() +
  labs(color='Age', title = "UMAP of expression according to age") +
  theme_classic()
p2
plot_grid(p, p2, ncol = 2)

## GENE ACTIVITY UMAP ## 
umapAUC <- RunUMAP(t(AUCdata), dims = 13, n.neighbors = 30) #AUCdata is the regulon activity dataframe 
umapdata <- as.data.frame(umapAUC@cell.embeddings)
colnames(umapdata) <- c("umap1", "umap2") 
rownames(umapdata) <- colnames(AUCdata)
umapdata$cell_state <-  seuratObj_umapdata$clusters #add extra column next to UMAP_1 and UMAP_2 which you want to use as colored clusters in your plot 
p<-ggplot(data = umapdata, aes(x = umap1, y = umap2, color = cell_state )) + 
  geom_point(shape = ".") +  
  labs(color='Cell type', title = "UMAP of regulon data") +
  theme_classic()
p

p1 <- ggplot(data = umapdata, aes(x = umap1, y = umap2, color = cell_state )) + 
  geom_point() +  
  labs(color='Cell type', title = "UMAP of regulon data") +
  theme_classic()
p1

plot_grid(p, p1, ncol = 2)


## COLOR REGULON ACTIVITY ON UMAP ##
# https://github.com/satijalab/seurat/issues/3521 
network_AD <- read.table(paste0(path_to_loom,"edgelist_sc_AD.txt"), sep = "\t", header = T)
reg_IKZF1 <- network_AD %>% filter(regulator == "IKZF1")
reg_IRF8 <- network_AD %>% filter(regulator == "IRF8")
reg_NFATC2 <- network_AD %>% filter(regulator == "NFATC2")
reg_RUNX1 <- network_AD %>% filter(regulator == "RUNX1")
reg_TAL1 <- network_AD %>% filter(regulator == "TAL1")
reg_SPI1 <- network_AD %>% filter(regulator == "SPI1")

regulongenes <- list(reg_TAL1$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "TAL1_regulon")
p1 <- FeaturePlot(object = seuratObj, features = "TAL1_regulon1") + ggtitle("TAL1_regulon")
regulongenes <- list(reg_IKZF1$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "IKZF1_regulon")
p2 <- FeaturePlot(object = seuratObj, features = "IKZF1_regulon1") + ggtitle("IKZF1_regulon")
regulongenes <- list(reg_IRF8$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "IRF8_regulon")
p3 <- FeaturePlot(object = seuratObj, features = "IRF8_regulon1") + ggtitle("IRF8_regulon")
regulongenes <- list(reg_NFATC2$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "NFATC2_regulon")
p4 <- FeaturePlot(object = seuratObj, features = "NFATC2_regulon1") + ggtitle("NFATC2_regulon")
regulongenes <- list(reg_RUNX1$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "RUNX1_regulon")
p5 <- FeaturePlot(object = seuratObj, features = "RUNX1_regulon1") + ggtitle("RUNX1_regulon")
regulongenes <- list(reg_SPI1$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "SPI1_regulon")
p6 <- FeaturePlot(object = seuratObj, features = "SPI1_regulon1") + ggtitle("SPI1_regulon")

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)

# module 35: MYRF, SOX10, ZNF536, PROX1, CREB5, ST18, SOX8, RBPJ, ZBED3, BBX
reg_MYRF <- network_AD %>% filter(regulator == "MYRF") # no regulon
reg_SOX10 <- network_AD %>% filter(regulator == "SOX10")
reg_ZNF536 <- network_AD %>% filter(regulator == "ZNF536") 
reg_PROX1 <- network_AD %>% filter(regulator == "PROX1") # no regulon
reg_CREB5 <- network_AD %>% filter(regulator == "CREB5") # no regulon
reg_ST18 <- network_AD %>% filter(regulator == "ST18") # no regulon
reg_SOX8 <- network_AD %>% filter(regulator == "SOX8") 
reg_RBPJ <- network_AD %>% filter(regulator == "RBPJ") 
reg_ZBED3 <- network_AD %>% filter(regulator == "ZBED3") # no regulon
reg_BBX <- network_AD %>% filter(regulator == "BBX") # no regulon

regulongenes <- list(reg_SOX10$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "SOX10_regulon")
p1 <- FeaturePlot(object = seuratObj, features = "SOX10_regulon1") + ggtitle("SOX10_regulon")
regulongenes <- list(reg_ZNF536$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "ZNF536_regulon")
p2 <- FeaturePlot(object = seuratObj, features = "ZNF536_regulon1") + ggtitle("ZNF536_regulon")
regulongenes <- list(reg_SOX8$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "SOX8_regulon")
p4 <- FeaturePlot(object = seuratObj, features = "SOX8_regulon1") + ggtitle("SOX8_regulon")
regulongenes <- list(reg_RBPJ$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "RBPJ_regulon")
p5 <- FeaturePlot(object = seuratObj, features = "RBPJ_regulon1") + ggtitle("RBPJ_regulon")

plot_grid(p1, p2, p4, p5, ncol = 2, nrow = 2)

### CSI clusters on UMAPs ###
regulonclusters_AD <- read.table(paste0(path_to_loom,"regulonclusters15_AD.txt"), header = T)
regulonclusters_AD$regulon <- regulonclusters_AD$regulon %>% str_replace_all("_.*", "")
cluster1 <- regulonclusters_AD %>% filter(csi_cluster == 1)
cluster2 <- regulonclusters_AD %>% filter(csi_cluster == 2)
cluster3 <- regulonclusters_AD %>% filter(csi_cluster == 3)
cluster4 <- regulonclusters_AD %>% filter(csi_cluster == 4)
cluster5 <- regulonclusters_AD %>% filter(csi_cluster == 5)
cluster6 <- regulonclusters_AD %>% filter(csi_cluster == 6)
cluster7 <- regulonclusters_AD %>% filter(csi_cluster == 7)
cluster8 <- regulonclusters_AD %>% filter(csi_cluster == 8)
cluster9 <- regulonclusters_AD %>% filter(csi_cluster == 9)
cluster10 <- regulonclusters_AD %>% filter(csi_cluster == 10)
cluster11 <- regulonclusters_AD %>% filter(csi_cluster == 11)
cluster12 <- regulonclusters_AD %>% filter(csi_cluster == 12)
cluster13 <- regulonclusters_AD %>% filter(csi_cluster == 13)
cluster14 <- regulonclusters_AD %>% filter(csi_cluster == 14)
cluster15 <- regulonclusters_AD %>% filter(csi_cluster == 15)

reg_1 <- network_AD %>% filter(regulator %in% cluster1$regulon)
reg_2 <- network_AD %>% filter(regulator %in% cluster2$regulon)
reg_3 <- network_AD %>% filter(regulator %in% cluster3$regulon)
reg_4 <- network_AD %>% filter(regulator %in% cluster4$regulon)
reg_5 <- network_AD %>% filter(regulator %in% cluster5$regulon)
reg_6 <- network_AD %>% filter(regulator %in% cluster6$regulon)
reg_7 <- network_AD %>% filter(regulator %in% cluster7$regulon)
reg_8 <- network_AD %>% filter(regulator %in% cluster8$regulon)
reg_9 <- network_AD %>% filter(regulator %in% cluster9$regulon)
reg_10 <- network_AD %>% filter(regulator %in% cluster10$regulon)
reg_11 <- network_AD %>% filter(regulator %in% cluster11$regulon)
reg_12 <- network_AD %>% filter(regulator %in% cluster12$regulon)
reg_13 <- network_AD %>% filter(regulator %in% cluster13$regulon)
reg_14 <- network_AD %>% filter(regulator %in% cluster14$regulon)
reg_15 <- network_AD %>% filter(regulator %in% cluster15$regulon)

regulongenes <- list(reg_1$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_1")
p1 <- FeaturePlot(object = seuratObj, features = "CSI_module_11") + ggtitle("CSI_module_1")
regulongenes <- list(reg_2$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_2")
p2 <- FeaturePlot(object = seuratObj, features = "CSI_module_21") + ggtitle("CSI_module_2")
regulongenes <- list(reg_3$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_3")
p3 <- FeaturePlot(object = seuratObj, features = "CSI_module_31") + ggtitle("CSI_module_3")
regulongenes <- list(reg_4$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_4")
p4 <- FeaturePlot(object = seuratObj, features = "CSI_module_41") + ggtitle("CSI_module_4")
regulongenes <- list(reg_5$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_5")
p5 <- FeaturePlot(object = seuratObj, features = "CSI_module_51") + ggtitle("CSI_module_5")
regulongenes <- list(reg_6$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_6")
p6 <- FeaturePlot(object = seuratObj, features = "CSI_module_61") + ggtitle("CSI_module_6")
regulongenes <- list(reg_7$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_7")
p7 <- FeaturePlot(object = seuratObj, features = "CSI_module_71") + ggtitle("CSI_module_7")
regulongenes <- list(reg_8$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_8")
p8 <- FeaturePlot(object = seuratObj, features = "CSI_module_81") + ggtitle("CSI_module_8")
regulongenes <- list(reg_9$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_9")
p9 <- FeaturePlot(object = seuratObj, features = "CSI_module_91") + ggtitle("CSI_module_9")
regulongenes <- list(reg_10$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_10")
p10 <- FeaturePlot(object = seuratObj, features = "CSI_module_101") + ggtitle("CSI_module_10")
regulongenes <- list(reg_11$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_11")
p11 <- FeaturePlot(object = seuratObj, features = "CSI_module_111") + ggtitle("CSI_module_11")
regulongenes <- list(reg_12$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_12")
p12 <- FeaturePlot(object = seuratObj, features = "CSI_module_121") + ggtitle("CSI_module_12")
regulongenes <- list(reg_13$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_13")
p13 <- FeaturePlot(object = seuratObj, features = "CSI_module_131") + ggtitle("CSI_module_13")
regulongenes <- list(reg_14$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_14")
p14 <- FeaturePlot(object = seuratObj, features = "CSI_module_141") + ggtitle("CSI_module_14")
regulongenes <- list(reg_15$gene)
seuratObj <- AddModuleScore(object = seuratObj, features = regulongenes, name = "CSI_module_15")
p15 <- FeaturePlot(object = seuratObj, features = "CSI_module_151") + ggtitle("CSI_module_15")

plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, ncol = 4, nrow = 4)

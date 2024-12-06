#### Cell type proportions in samples of single-cell datasets AD and MDD 
# <https://bioinformatics.stackexchange.com/questions/11150/percentage-of-each-cluster-in-seurat>  

# Packages and data
library(Seurat)
library(ggplot2)
setwd("~/MyFiles/Master/single-cell")
seuratobj_AD <- readRDS("Seurat_AD.rds")
seuratobj_MDD <- readRDS("Seurat_MDD.rds")

table(seuratobj_AD$celltype, seuratobj_AD$Diagnosis)
seuratobj_AD <- SetIdent(seuratobj_AD, value = "celltype")

prop.table(table(Idents(seuratobj_AD), seuratobj_AD$Diagnosis), margin = 2)

pt <- table(Idents(seuratobj_AD), seuratobj_AD$Diagnosis)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
names(pt) <- c("Celltype", "Diagnosis", "Freq")
ggplot(pt, aes(x = Diagnosis, y = Freq, fill = Celltype)) +
  geom_col(position = "fill", width = 0.5) +
  labs(title = "Cell type proportions",
       x = "Diagnosis",
       y = "Proportion",
       fill = "") +
  scale_fill_discrete(type = c('#882255',"#D55E00","#E69F00","#F0E442",
                               '#DDCC77', '#999933','#CC6677')) +
  theme_light() +
  theme(axis.line.x = element_blank())
ggsave(filename = "Figures/Celltype_proportions_AD.pdf", device = "pdf")

######## MDD #########
table(seuratobj_MDD$Celltype, seuratobj_MDD$Diagnosis)
seuratobj_MDD <- SetIdent(seuratobj_MDD, value = "Celltype")

prop.table(table(Idents(seuratobj_MDD), seuratobj_MDD$Diagnosis), margin = 2)
seuratobj_MDD <- subset(seuratobj_MDD, subset = Celltype != "MIX")

pt <- table(Idents(seuratobj_MDD), seuratobj_MDD$Diagnosis)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
names(pt) <- c("Celltype", "Diagnosis", "Freq")
ggplot(pt, aes(x = Diagnosis, y = Freq, fill = Celltype)) +
  geom_col(position = "fill", width = 0.5) +
  labs(title = "Cell type proportions",
       x = "Diagnosis",
       y = "Proportion",
       fill = "") +
  scale_fill_discrete(type = c('#882255','#CC6677',"#D55E00","#E69F00","#F0E442",
                               '#DDCC77', '#999933')) +
  scale_x_discrete(limits = c("Suicide", "Control"), labels = c("MDD", "Control")) +
  theme_light() +
  theme(axis.line.x = element_blank())
ggsave(filename = "Figures/Celltype_proportions_MDD.pdf", device = "pdf")

# Statistical test for different proportion of astrocytes in MDD vs control
# <https://bioconductor.org/packages/release/bioc/vignettes/speckle/inst/doc/speckle.html> 
table(seuratobj_MDD$Sample_ID, seuratobj_MDD$Diagnosis)

#BiocManager::install("speckle")
library(speckle)
library(limma)
library(scater)
library(patchwork)
library(edgeR)
library(statmod)

# defauls is logit transformation 
propeller(seuratobj_MDD, clusters = seuratobj_MDD$Celltype, sample = seuratobj_MDD$Sample_ID, 
          group = seuratobj_MDD$Diagnosis)
# Astrocytes is significantly different! The other cell types are not
# BaselineProp.clusters BaselineProp.Freq PropMean.Control PropMean.Suicide PropRatio Tstatistic      P.Value         FDR
#            ASC       0.051680721      0.084174767       0.03185455      2.6424727   3.6999087 0.0006876235  0.005500988


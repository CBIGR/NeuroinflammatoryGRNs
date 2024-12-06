# single-cell MDD dataset GSE144136

#### single-cell preprocessing MDD dataset ###
# packages
library(tidyverse, lib.loc = '/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/Rlib')
library(Seurat)
library(edgeR)
library(biomaRt)
library(Matrix)

#load dataset 
rm(list=ls()); gc() 
counts_MDD.data <- ReadMtx(mtx = "GSE144136/GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz", 
                           cells = "GSE144136/GSE144136_CellNames.csv.gz", 
                           features = "GSE144136/GSE144136_GeneNames.csv.gz", feature.column = 1)

# create SeuratObject 
counts_MDD <- CreateSeuratObject(counts = counts_MDD.data)

# Examine QC metrics
print(c("median of nCount_RNA", median(counts_MDD@meta.data$nCount_RNA))) #UMI counts 2445
print(c("median of nFeature_RNA", median(counts_MDD@meta.data$nFeature_RNA))) #gene counts 1622

hist(counts_MDD$nCount_RNA)
hist(counts_MDD$nFeature_RNA)

VlnPlot(counts_MDD, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) 

FeatureScatter(counts_MDD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

# protein-coding genes
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", 
                host = 'www.ensembl.org')
genes_protcod <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", 
                                               "transcript_biotype"), 
                                filters = c("transcript_biotype", "chromosome_name"),
                                values = list("protein_coding", c(1:22)), mart = mart)
head(genes_protcod)
hgnc_symbol <- genes_protcod$hgnc_symbol

# QC subset
dim(counts_MDD)
counts_MDD <- subset(counts_MDD, subset = nFeature_RNA > 200 &
                         nFeature_RNA < 5000)
dim(counts_MDD)

# select protein-coding genes 
# put genes in vector
genes_MDD <- rownames(counts_MDD)
head(genes_MDD)
# select protein-coding genes 
genes_MDD <- genes_MDD[genes_MDD %in% hgnc_symbol]
length(genes_MDD)
# subset
counts_MDD <- subset(counts_MDD, features = genes_MDD)
dim(counts_MDD)

# Normalize, highly variable features and scaling
counts_MDD <- NormalizeData(counts_MDD)
counts_MDD <- FindVariableFeatures(counts_MDD, nfeatures = 8000)
variablefeatures_MDD <- VariableFeatures(counts_MDD)
length(variablefeatures_MDD)

regulators <- read.table("dbTF.csv", header = T)
dim(regulators)
regulators_MDD <- genes_MDD[genes_MDD %in% regulators$Gene_symbol]
length(regulators_MDD)
# combine vector variablefeatures and regulators
genestoselect_MDD <- c(variablefeatures_MDD, regulators_MDD)
genestoselect_MDD <- unique(genestoselect_MDD)
length(genestoselect_MDD)
counts_MDD <- subset(counts_MDD, features = genestoselect_MDD)
dim(counts_MDD)

# Scale counts
counts_MDD <- ScaleData(counts_MDD)
#saveRDS(counts_MDD, "Seurat_MDD.rds")
counts_MDD.data <- GetAssayData(object = counts_MDD, slot = "data")
dim(counts_MDD.data) #8638 77437
write.table(counts_MDD.data, file = "counts_MDD.txt", sep = "\t", row.names = T, col.names = T)

#PCA
counts_MDD <- RunPCA(counts_MDD)
ElbowPlot(counts_MDD, ndims = 50)
#DimPlot(counts_MDD)

# UMAP 
counts_MDD <- FindNeighbors(counts_MDD, dims = 1:15) #first 15 PCs
counts_MDD <- FindClusters(counts_MDD, resolution = 0.5) 
counts_MDD <- RunUMAP(counts_MDD, dims = 1:15)
DimPlot(counts_MDD, reduction = "umap")
# Run again for better UMAP with 50 PCs
counts_MDD <- FindNeighbors(counts_MDD, dims = 1:50) #first 50 PCs
counts_MDD <- FindClusters(counts_MDD, resolution = 1, algorithm = 4) 
counts_MDD <- RunUMAP(counts_MDD, dims = 1:50)
DimPlot(counts_MDD, reduction = "umap", group.by = "Celltype")

# metadata
CellNames_MDD <- read.csv("GSE144136/GSE144136_CellNames.csv.gz", header = F)
head(CellNames_MDD)
CellNames_MDD[, c("Celltype", "Rest")] <- CellNames_MDD$V1 %>% str_split_fixed("[.]", 2)
CellNames_MDD[, c("Sample_ID", "Diagnosis", "Batch", "UMI")] <- CellNames_MDD$Rest %>% str_split_fixed("_", 4)
head(CellNames_MDD)
CellNames_MDD <- column_to_rownames(CellNames_MDD, "V1")
CellNames_MDD <- CellNames_MDD %>% dplyr::select(-Rest)
CellNames_MDD$Celltype <- gsub("_.*","",CellNames_MDD$Celltype)
CellNames_MDD <- CellNames_MDD %>% mutate(Celltype = case_when(Celltype == "Micro/Macro" ~ "MICRO",
                                                               Celltype == "Astros" ~ "ASC",
                                                               Celltype == "Endo" ~ "END",
                                                               Celltype == "Ex" ~ "EX",
                                                               Celltype == "Inhib" ~ "INH",
                                                               Celltype == "Mix"~ "MIX",
                                                               Celltype == "OPCs" ~ "OPC",
                                                               Celltype == "Oligos" ~ "ODC",
                                                               T ~ Celltype))
counts_MDD <- AddMetaData(counts_MDD, CellNames_MDD)
head(counts_MDD@meta.data)
#saveRDS(counts_MDD, "Seurat_MDD.rds")

# single-cell AD dataset GSE174367

#### Single-cell preprocessing AD ####
# packages
library(tidyverse, lib.loc = '/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/Rlib')
library(Seurat)
library(edgeR)
library(biomaRt)
library(Matrix)

#load dataset 
rm(list=ls()); gc()
counts_AD.data <- ReadMtx(mtx = "snRNA_counts.mtx", 
                             features = "genes.csv",
                             cells = "barcodes_rna.csv", feature.column = 1)

# create SeuratObject 
counts_AD <- CreateSeuratObject(counts = counts_AD.data)

# Examine QC metrics
print(c("median of nCount_RNA", median(counts_AD@meta.data$nCount_RNA))) #UMI counts 6382
print(c("median of nFeature_RNA", median(counts_AD@meta.data$nFeature_RNA))) #gene counts 2470

hist(counts_AD$nCount_RNA)
hist(counts_AD$nFeature_RNA)

VlnPlot(counts_AD, features = c("nFeature_RNA", "nCount_RNA"),
               ncol = 2) 

FeatureScatter(counts_AD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

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
dim(counts_AD)
counts_AD <- subset(counts_AD, subset = nFeature_RNA > 200 &
                       nFeature_RNA < 7500)
dim(counts_AD)

# select protein-coding genes 
# put genes in vector
genes_AD <- rownames(counts_AD)
head(genes_AD)
# select protein-coding genes 
genes_AD <- genes_AD[genes_AD %in% hgnc_symbol]
length(genes_AD)
# subset
counts_AD <- subset(counts_AD, features = genes_AD)
dim(counts_AD)

# Normalize, highly variable features and scaling
counts_AD <- NormalizeData(counts_AD)
counts_AD <- FindVariableFeatures(counts_AD, nfeatures = 8000)
variablefeatures_AD <- VariableFeatures(counts_AD)
length(variablefeatures_AD)

regulators <- read.table("dbTF.csv", header = T)
dim(regulators)
regulators_AD <- genes_AD[genes_AD %in% regulators$Gene_symbol]
length(regulators_AD)

# combine vector variablefeatures and regulators
genestoselect_AD <- c(variablefeatures_AD, regulators_AD)
genestoselect_AD <- unique(genestoselect_AD)
length(genestoselect_AD)
counts_AD <- subset(counts_AD, features = genestoselect_AD)
dim(counts_AD)

# Scale counts 
counts_AD <- ScaleData(counts_AD)
#saveRDS(counts_AD, "Seurat_AD.rds")
counts_AD.data <- GetAssayData(object = counts_AD, slot = "data")
dim(counts_AD.data) #8655 59968
write.table(counts_AD.data, file = "counts_AD.txt", sep = "\t", row.names = T, col.names = T)

#PCA
counts_AD <- RunPCA(counts_AD)
ElbowPlot(counts_AD, ndims = 50)
#DimPlot(counts_AD)

# UMAP 
counts_AD <- FindNeighbors(counts_AD, dims = 1:13) #first 13 PCs
counts_AD <- FindClusters(counts_AD, resolution = 0.5) 
counts_AD <- RunUMAP(counts_AD, dims = 1:13)
DimPlot(counts_AD, reduction = "umap")

# metadata
snRNA_metadta <- read_csv("snRNA_metadta.csv")
snRNA_metadta <- snRNA_metadta %>% filter(...1 %in% colnames(counts_AD))
head(snRNA_metadta)
counts_AD[["celltype"]] <- snRNA_metadta$celltype
counts_AD[["Batch"]] <- snRNA_metadta$Batch
counts_AD[["cluster"]] <- snRNA_metadta$cluster
counts_AD[["Diagnosis"]] <- snRNA_metadta$Diagnosis
counts_AD[["Sample"]] <- snRNA_metadta$Sample.ID
head(counts_AD@meta.data)
#saveRDS(counts_AD, "Seurat_AD.rds")

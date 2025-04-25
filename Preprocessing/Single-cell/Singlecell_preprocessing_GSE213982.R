# single cell preprocessing MDD dataset
# packages
library(tidyverse)
library(Seurat)
library(edgeR)
library(biomaRt)
library(Matrix)
options(bitmapType = "cairo")
setwd("/kyukon/scratch/gent/vo/000/gvo00027/projects/CBIGR/21HPP_GRN_neuroinfl/singlecell/Female_samples_MDD")

#load dataset 
rm(list=ls()); gc() 
counts_MDD.data <- ReadMtx(mtx = "GSE213982/matrix.mtx.gz", 
                           cells = "GSE213982/barcodes.tsv.gz", 
                           features = "GSE213982/features.tsv.gz", feature.column = 1, skip.cell = 1, skip.feature = 1)

# create SeuratObject 
counts_MDD <- CreateSeuratObject(counts = counts_MDD.data)

# Filter for female samples 
head(counts_MDD@meta.data)
metadata_MDD <- counts_MDD@meta.data
metadata_MDD[,c("Sample", "Barcode", "Celltype", "Cellstate")] <- str_split_fixed(rownames(metadata_MDD), "\\.", 4)
#unique(metadata_MDD$Sample)[1:38]
counts_MDD$Sample <- metadata_MDD$Sample
counts_MDD$Celltype <- metadata_MDD$Celltype
counts_MDD$Cellstate <- metadata_MDD$Cellstate
Idents(counts_MDD) <- "Sample"
counts_MDD <- subset(counts_MDD, idents = unique(metadata_MDD$Sample)[1:38])
table(counts_MDD$Sample)

# Examine QC metrics
median(counts_MDD@meta.data$nCount_RNA) #UMI counts 2779
median(counts_MDD@meta.data$nFeature_RNA) #gene counts 1636

png("hist_counts_MDD_female.png")
hist(counts_MDD$nCount_RNA)
dev.off()
png("hist_features_MDD_female.png")
hist(counts_MDD$nFeature_RNA)
dev.off()

Idents(counts_MDD) <- "orig.ident"
plt <- VlnPlot(counts_MDD, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = F) 
ggsave(filename = "Vlnplot_MDD_female2.png", plot = plt, width = 7, height = 3.5)

plt2 <- FeatureScatter(counts_MDD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
ggsave(filename = "scatter_MDD_female.png", plot = plt2, width = 7, height = 3.5)

# protein-coding genes
#mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", 
#                host = 'www.ensembl.org')
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", mirror = "asia")
genes_protcod <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", 
                                               "transcript_biotype"), 
                                filters = c("transcript_biotype", "chromosome_name"),
                                values = list("protein_coding", c(1:22)), mart = mart)
#head(genes_protcod)
hgnc_symbol <- genes_protcod$hgnc_symbol

# QC subset
dim(counts_MDD)
max(counts_MDD$nFeature_RNA)
counts_MDD <- subset(counts_MDD, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA < 100000)
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

regulators <- read.table("../dbTF.csv", header = T)
dim(regulators)
regulators_MDD <- genes_MDD[genes_MDD %in% regulators$Gene_symbol]
length(regulators_MDD)
# combine vector variablefeatures and regulators
genestoselect_MDD <- c(variablefeatures_MDD, regulators_MDD)
genestoselect_MDD <- unique(genestoselect_MDD)
length(genestoselect_MDD)
counts_MDD_f <- subset(counts_MDD, features = genestoselect_MDD)
dim(counts_MDD_f)

counts_MDD <- ScaleData(counts_MDD, features = genestoselect_MDD)
saveRDS(counts_MDD, file = "Seurat_MDD_female.rds")
counts_MDD_f <- ScaleData(counts_MDD_f, features = genestoselect_MDD)
saveRDS(counts_MDD_f, file = "Seurat_MDD_M_F.rds")

# add metadata diagnosis
counts_MDD@meta.data <- counts_MDD@meta.data %>% mutate(Diagnosis = case_when(Sample == "F1" ~ "Suicide",
                                                                            Sample == "F2" ~ "Suicide",
                                                                            Sample == "F3" ~ "Suicide",
                                                                            Sample == "F4" ~ "Suicide",
                                                                            Sample == "F5" ~ "Suicide",
                                                                            Sample == "F6" ~ "Suicide",
                                                                            Sample == "F7" ~ "Control",
                                                                            Sample == "F8" ~ "Suicide",
                                                                            Sample == "F9" ~ "Suicide",
                                                                            Sample == "F10" ~ "Control",
                                                                            Sample == "F11" ~ "Suicide",
                                                                            Sample == "F12" ~ "Suicide",
                                                                            Sample == "F13" ~ "Control",
                                                                            Sample == "F14" ~ "Suicide",
                                                                            Sample == "F15" ~ "Suicide",
                                                                            Sample == "F16" ~ "Suicide",
                                                                            Sample == "F17" ~ "Suicide",
                                                                            Sample == "F18" ~ "Suicide",
                                                                            Sample == "F19" ~ "Suicide",
                                                                            Sample == "F20" ~ "Suicide",
                                                                            Sample == "F21" ~ "Control",
                                                                            Sample == "F22" ~ "Control",
                                                                            Sample == "F23" ~ "Control",
                                                                            Sample == "F24" ~ "Control",
                                                                            Sample == "F25" ~ "Suicide",
                                                                            Sample == "F26" ~ "Control",
                                                                            Sample == "F27" ~ "Suicide",
                                                                            Sample == "F28" ~ "Suicide",
                                                                            Sample == "F29" ~ "Control",
                                                                            Sample == "F30" ~ "Control",
                                                                            Sample == "F31" ~ "Control",
                                                                            Sample == "F32" ~ "Control",
                                                                            Sample == "F33" ~ "Control",
                                                                            Sample == "F34" ~ "Control",
                                                                            Sample == "F35" ~ "Control",
                                                                            Sample == "F36" ~ "Control",
                                                                            Sample == "F37" ~ "Control",
                                                                            Sample == "F38" ~ "Control",
                                                                            Sample == "M1" ~ "Suicide",
                                                                            Sample == "M2" ~ "Control",
                                                                            Sample == "M3" ~ "Control",
                                                                            Sample == "M4" ~ "Suicide",
                                                                            Sample == "M5" ~ "Suicide",
                                                                            Sample == "M6" ~ "Suicide",
                                                                            Sample == "M7" ~ "Control",
                                                                            Sample == "M8" ~ "Suicide",
                                                                            Sample == "M9" ~ "Control",
                                                                            Sample == "M10" ~ "Suicide",
                                                                            Sample == "M11" ~ "Suicide",
                                                                            Sample == "M12" ~ "Control",
                                                                            Sample == "M13" ~ "Control",
                                                                            Sample == "M14" ~ "Suicide",
                                                                            Sample == "M15" ~ "Control",
                                                                            Sample == "M16" ~ "Control",
                                                                            Sample == "M17" ~ "Suicide",
                                                                            Sample == "M18" ~ "Suicide",
                                                                            Sample == "M19" ~ "Control",
                                                                            Sample == "M20" ~ "Control",
                                                                            Sample == "M21" ~ "Control",
                                                                            Sample == "M22" ~ "Control",
                                                                            Sample == "M23" ~ "Suicide",
                                                                            Sample == "M24" ~ "Control",
                                                                            Sample == "M25" ~ "Control",
                                                                            Sample == "M26" ~ "Suicide",
                                                                            Sample == "M27" ~ "Control",
                                                                            Sample == "M28" ~ "Suicide",
                                                                            Sample == "M29" ~ "Control",
                                                                            Sample == "M30" ~ "Suicide",
                                                                            Sample == "M31" ~ "Control",
                                                                            Sample == "M32" ~ "Suicide",
                                                                            Sample == "M33" ~ "Suicide",
                                                                            Sample == "M34" ~ "Suicide"))


table(counts_MDD$Diagnosis)




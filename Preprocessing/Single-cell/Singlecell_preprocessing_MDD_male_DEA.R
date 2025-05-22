library(tidyverse)
library(Seurat)
library(edgeR)
library(biomaRt)
library(Matrix)
library(presto, lib.loc = "/kyukon/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/tmp/VZZ_presto/")

so <- readRDS("singlecell/Seurat_MDD.rds")
so

so <- FindVariableFeatures(so, nfeatures = 8000)
so

# protein-coding genes
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", 
                host = 'www.ensembl.org')
genes_protcod <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", 
                                               "transcript_biotype"), 
                                filters = c("transcript_biotype", "chromosome_name"),
                                values = list("protein_coding", c(1:22)), mart = mart)
#head(genes_protcod)
hgnc_symbol <- genes_protcod$hgnc_symbol

# select protein-coding genes 
# put genes in vector
genes <- rownames(so)
head(genes)

# select protein-coding genes 
genes <- genes[genes %in% hgnc_symbol]
length(genes)

# subset
so <- subset(so, features = genes)
dim(so)


so <- FindVariableFeatures(so, nfeatures = 8000)
so

regulators <- read.table("singlecell/dbTF.csv", header = T)
dim(regulators)

regulators <- genes[genes %in% regulators$Gene_symbol]
length(regulators)

# combine vector variablefeatures and regulators
genestoselect <- c(VariableFeatures(so), regulators)
genestoselect <- unique(genestoselect)
length(genestoselect)

so <- subset(so, features = genestoselect)
dim(so)

head(so@meta.data)

# Import metadata

metadata <- read.table("singlecell/Metadata_MDD.txt", sep = "\t", header=T, row.names=1)
metadata

so@meta.data <- cbind(so@meta.data[,1:5], metadata[Cells(so),])
so@meta.data

# Check cell labels
table(so$Celltype)

# Include diagnosis
table(so$Diagnosis)

# Perform DEA
so$Celltype_byDiagnosis <- paste0(so$Celltype, so$Diagnosis)
table(so$Celltype_byDiagnosis)

so<-SetIdent(so, value = "Celltype_byDiagnosis")

deg_tests <- list()
for (ct in levels(as.factor(so$Celltype))) {
  
  degs <- FindMarkers(so, ident.1 = paste0(ct, "Suicide"), 
                      ident.2 = paste0(ct, "Control"))
  
  degs <- degs[( abs(degs$avg_log2FC) > 0.25) & (degs$p_val_adj) < 0.05, ]
  
  deg_tests[[ct]] <- degs
}

# Check results

names(deg_tests)

for (ct in names(deg_tests)){
  deg_tests[[ct]]$contrast <- paste0(ct,"_Suicide_vs_Control")
  deg_tests[[ct]]$gene_name <- row.names(deg_tests[[ct]])
  
}

write.csv(bind_rows(deg_tests, .id = "column_label"), file = "singlecell/MDD_male_DEGs_MDDvC_byCelltype.csv", quote = F, row.names = F)


# With pseudobulking to reduce false positives
pseudo_mdd <- AggregateExpression(so, assays = "RNA", return.seurat = T, group.by = c("Diagnosis", "Sample_ID", "Celltype"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_mdd))

pseudo_mdd$celltype.condition <- paste0(pseudo_mdd$Celltype, pseudo_mdd$Diagnosis)

pseudo_mdd <- SetIdent(pseudo_mdd, value = "celltype.condition")

bulk.ast.de <- FindMarkers(object = pseudo_mdd, 
                           ident.1 = "ASCSuicide", 
                           ident.2 = "ASCControl",
                           test.use = "DESeq2")
head(bulk.ast.de, n = 15)



deg_tests <- list()
for (ct in levels(as.factor(so$Celltype))) {
  if (ct == "END"){
    next #because too little cells and results in bad pseudobulks (with zeros, that are not accepted by DESeq2)
  }
  
  message("Calculating for ", ct)
  degs <- FindMarkers(object = pseudo_mdd, 
                      ident.1 = paste0(ct,"Suicide"), 
                      ident.2 = paste0(ct,"Control"),
                      test.use = "DESeq2")
  
  # degs <- degs[( abs(degs$avg_log2FC) > 0.25) & (degs$p_val_adj < 0.05), ]
  
  deg_tests[[ct]] <- degs
}

for (ct in names(deg_tests)){
  deg_tests[[ct]]$contrast <- paste0(ct,"_Suicide_vs_Control")
  deg_tests[[ct]]$gene_name <- row.names(deg_tests[[ct]])
  
}


degs_total <- bind_rows(deg_tests, .id = "column_label")
dim(degs_total)

degs_total <- degs_total[degs_total$p_val_adj < 0.05,]
degs_total <- degs_total[! is.na(degs_total$column_label),]
dim(degs_total)

write.csv(degs_total, file = "singlecell/MDD_male_DEGs_MDDvC_byCelltype_pb.csv", quote = F, row.names = F)



########################## Repear with AD #################################################


so <- readRDS("singlecell/Seurat_AD1.rds")
so

so <- FindVariableFeatures(so, nfeatures = 8000)
so

# select protein-coding genes 
# put genes in vector
genes <- rownames(so)
head(genes)

# select protein-coding genes 
genes <- genes[genes %in% hgnc_symbol]
length(genes)

# subset
so <- subset(so, features = genes)
dim(so)


so <- FindVariableFeatures(so, nfeatures = 8000)
so

regulators <- read.table("singlecell/dbTF.csv", header = T)
dim(regulators)

regulators <- genes[genes %in% regulators$Gene_symbol]
length(regulators)

# combine vector variablefeatures and regulators
genestoselect <- c(VariableFeatures(so), regulators)
genestoselect <- unique(genestoselect)
length(genestoselect)

so <- subset(so, features = genestoselect)
dim(so)

head(so@meta.data)

# Import metadata

metadata <- read.table("singlecell/Metadata_MDD.txt", sep = "\t", header=T, row.names=1)
metadata

so@meta.data <- cbind(so@meta.data[,1:5], metadata[Cells(so),])
so@meta.data

# Check cell labels
table(so$Celltype)

# Include diagnosis
table(so$Diagnosis)

# Perform DEA
so$Celltype_byDiagnosis <- paste0(so$Celltype, so$Diagnosis)
table(so$Celltype_byDiagnosis)

so<-SetIdent(so, value = "Celltype_byDiagnosis")

deg_tests <- list()
for (ct in levels(as.factor(so$Celltype))) {
  
  degs <- FindMarkers(so, ident.1 = paste0(ct, "Suicide"), 
                      ident.2 = paste0(ct, "Control"))
  
  degs <- degs[( abs(degs$avg_log2FC) > 0.25) & (degs$p_val_adj) < 0.05, ]
  
  deg_tests[[ct]] <- degs
}

# Check results

names(deg_tests)

for (ct in names(deg_tests)){
  deg_tests[[ct]]$contrast <- paste0(ct,"_Suicide_vs_Control")
  deg_tests[[ct]]$gene_name <- row.names(deg_tests[[ct]])
  
}

write.csv(bind_rows(deg_tests, .id = "column_label"), file = "singlecell/MDD_male_DEGs_MDDvC_byCelltype.csv", quote = F, row.names = F)


# With pseudobulking to reduce false positives
pseudo_mdd <- AggregateExpression(so, assays = "RNA", return.seurat = T, group.by = c("Diagnosis", "Sample_ID", "Celltype"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_mdd))

pseudo_mdd$celltype.condition <- paste0(pseudo_mdd$Celltype, pseudo_mdd$Diagnosis)

pseudo_mdd <- SetIdent(pseudo_mdd, value = "celltype.condition")

bulk.ast.de <- FindMarkers(object = pseudo_mdd, 
                           ident.1 = "ASCSuicide", 
                           ident.2 = "ASCControl",
                           test.use = "DESeq2")
head(bulk.ast.de, n = 15)



deg_tests <- list()
for (ct in levels(as.factor(so$Celltype))) {
  if (ct == "END"){
    next #because too little cells and results in bad pseudobulks (with zeros, that are not accepted by DESeq2)
  }
  
  message("Calculating for ", ct)
  degs <- FindMarkers(object = pseudo_mdd, 
                      ident.1 = paste0(ct,"Suicide"), 
                      ident.2 = paste0(ct,"Control"),
                      test.use = "DESeq2")
  
  # degs <- degs[( abs(degs$avg_log2FC) > 0.25) & (degs$p_val_adj < 0.05), ]
  
  deg_tests[[ct]] <- degs
}

for (ct in names(deg_tests)){
  deg_tests[[ct]]$contrast <- paste0(ct,"_Suicide_vs_Control")
  deg_tests[[ct]]$gene_name <- row.names(deg_tests[[ct]])
  
}


degs_total <- bind_rows(deg_tests, .id = "column_label")
dim(degs_total)

degs_total <- degs_total[degs_total$p_val_adj < 0.05,]
degs_total <- degs_total[! is.na(degs_total$column_label),]
dim(degs_total)

write.csv(degs_total, file = "singlecell/MDD_male_DEGs_MDDvC_byCelltype_pb.csv", quote = F, row.names = F)

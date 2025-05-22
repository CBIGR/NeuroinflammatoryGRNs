library(tidyverse)
library(Seurat)
library(edgeR)
library(biomaRt)
library(Matrix)
library(presto, lib.loc = "/kyukon/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/tmp/VZZ_presto/")

so <- readRDS("/kyukon/scratch/gent/vo/000/gvo00027/projects/CBIGR/21HPP_GRN_neuroinfl/singlecell/Female_samples_MDD/Seurat_MDD_M_F.rds")
metadata_complete <- so@meta.data

# Save separately
write.csv(metadata_complete, "singlecell/MDD_M_F_metadata.csv", quote=F)

so <- readRDS("/kyukon/scratch/gent/vo/000/gvo00027/projects/CBIGR/21HPP_GRN_neuroinfl/singlecell/Female_samples_MDD/Seurat_MDD_female.rds")
so

# Check HVGs
head(VariableFeatures(so))

# Check metadata
colnames(so@meta.data)

# Check cell labels
table(so$Cellstate)
table(so$Celltype)
# Include diagnosis
so$Diagnosis <- metadata_complete[Cells(so), "Diagnosis"]
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

write.csv(bind_rows(deg_tests, .id = "column_label"), file = "singlecell/MDD_female_DEGs_MDDvC_byCelltype.csv", quote = F, row.names = F)


# With pseudobulking to reduce false positives

pseudo_mdd <- AggregateExpression(so, assays = "RNA", return.seurat = T, group.by = c("Diagnosis", "Sample", "Celltype"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_mdd))

pseudo_mdd$celltype.condition <- paste0(pseudo_mdd$Celltype, pseudo_mdd$Diagnosis)

pseudo_mdd <- SetIdent(pseudo_mdd, value = "celltype.condition")

bulk.ast.de <- FindMarkers(object = pseudo_mdd, 
                            ident.1 = "AstControl", 
                            ident.2 = "AstSuicide",
                            test.use = "DESeq2")
head(bulk.ast.de, n = 15)



deg_tests <- list()
for (ct in levels(as.factor(so$Celltype))) {
  
  degs <- FindMarkers(object = pseudo_mdd, 
                      ident.1 = paste0(ct,"Suicide"), 
                      ident.2 = paste0(ct,"Control"),
                      test.use = "DESeq2")
  
  degs <- degs[( abs(degs$avg_log2FC) > 0.25) & (degs$p_val_adj < 0.05), ]
  
  deg_tests[[ct]] <- degs
}

for (ct in names(deg_tests)){
  deg_tests[[ct]]$contrast <- paste0(ct,"_Suicide_vs_Control")
  deg_tests[[ct]]$gene_name <- row.names(deg_tests[[ct]])
  
}

write.csv(bind_rows(deg_tests, .id = "column_label"), file = "singlecell/MDD_female_DEGs_MDDvC_byCelltype_pb.csv", quote = F, row.names = F)

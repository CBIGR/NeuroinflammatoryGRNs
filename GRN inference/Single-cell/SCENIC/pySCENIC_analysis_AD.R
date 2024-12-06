######################################################################################################################

###############################################################################
# CREATE PLOTS FROM pySCENIC OUTPUT
###############################################################################

library("SingleCellExperiment")
library("SCENIC", lib.loc="/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/Rlib") ### installation see https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Setup.html#installation 
library("BBmisc")
library("SCopeLoomR", lib.loc="/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/Rlib") #####
library("Seurat")
library("dplyr")
library("ggplot2")
library("pheatmap")
library('openxlsx')
library("stringr")
library(data.table)

########################################
##### Getwd SAM2and3_WT_subset
########################################

#setwd("/scratch/gent/vo/000/gvo00027/projects/CBIGR/21JDS_singlecell_GBM/SCENIC/SCENIC_vsn/results/HPU_SCENIC_thesis/counts_MDD")

sampleName <- "HPU_counts_AD1" #Change for this analysis!!!
sampleFolder<-"/scratch/gent/vo/000/gvo00027/projects/CBIGR/21JDS_singlecell_GBM/SCENIC/SCENIC_vsn/results/HPU_SCENIC_thesis/counts_AD1/"

##add some subfolders
dir.create(paste0(sampleFolder,"results/"))

counts_AD1 <- fread("/kyukon/scratch/gent/vo/000/gvo00027/projects/CBIGR/21HPU_GRN_neuroinfl/singlecell/SCENIC/input/counts_AD1.txt", data.table=FALSE)
rownames(counts_AD1) <- counts_AD1[,1]
counts_AD1[,1] <- NULL
print(counts_AD1[1:5,1:5])
seuratObj <- CreateSeuratObject(counts_AD1)
cell_ann <- fread("/scratch/gent/vo/000/gvo00027/projects/CBIGR/21HPU_GRN_neuroinfl/singlecell/snRNA_metadta.csv", data.table=FALSE)

sum(rownames(seuratObj@meta.data) %in% str_replace(cell_ann[,1], "-", "\\."))
##### Read object --> JDS: same object as given as input to SCENIC
celltypes <- cell_ann[match(rownames(seuratObj@meta.data),str_replace(cell_ann[,1], "-", "\\.")),10]
seuratObj$cell_type <- celltypes

seuratObj <- SetIdent(seuratObj, value = "cell_type")

## Prepro 
######### 1. CREATE NEW CLUSTERS ##########
######################################################################################################################################
set.seed(123)

cellInfo <- as.data.frame(seuratObj@active.ident)
dim(cellInfo)
#1386    1

# Read data from loom file
loom <- open_loom(file.path = paste0(sampleFolder,"counts_AD1/scenic/counts_AD1/SCENIC_output.loom"))
regulonAUC <- get_regulons_AUC(loom)
thresholds <- get_regulon_thresholds(loom)
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
close_loom(loom)

dim(regulonAUC)
# 353 1386
length(thresholds)
# 380
length(regulons)
# 353


head(regulonAUC)
head(thresholds)
#head(regulons$`Ar(+)`)
#length(regulons$`Ar(+)`)


##Save regulon matrix
toExport<-as.data.frame(t(regulons_incidMat), stringsAsFactors=F)
toExport$gene<-rownames(toExport)
toExport<-toExport[,c('gene',setdiff(colnames(toExport),"gene"))]
tmpList<-list(toExport)
names(tmpList)<-"regulonMatrix"
#write.xlsx(tmpList, file = paste0(sampleFolder,"results/tableRegulons.xlsx"))


# FIX THRESHOLDS: names and values switched in output. pySCENIC/loom bug?
thresholds_values <- names(thresholds)
thresholds_names <- thresholds
thresholds <- thresholds_values
names(thresholds) <- thresholds_names
length(thresholds)
# 380
length(regulons)
# 353
head(thresholds)

#JDS: normally the thresholds and regulons should have the same dimensions: I will filter the regulons not in common out
sum(names(regulons) %in% names(thresholds))
which(names(thresholds) %in% names(regulons)) # all 353 first regulons in the threshold are present 
which(names(thresholds) %in% names(regulons)) # control 
thresholds <- thresholds[which(names(thresholds) %in% names(regulons))]
length(thresholds)
# 353
length(regulons)
# 353
head(thresholds)

# Get AUC values in matrix + scale to range 0-1 for nicer heatmaps
AUCdata <- getAUC(regulonAUC) %>% BBmisc::normalize("range")
dim(AUCdata)
# 353 1386
AUCdata[1:5,1:5]

# Define regulons per cell
regulonsCells <- setNames(lapply(names(thresholds), 
                                 function(x) {
                                   trh <- thresholds[x]
                                   names(which(getAUC(regulonAUC)[x,]>trh))
                                 }),names(thresholds))

length(regulonsCells)
#353

regulonActivity <- reshape2::melt(regulonsCells)
head(regulonActivity)

# Create binary activity (regulon On or Off based on thresholds)
binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
class(binaryRegulonActivity) <- "matrix"
dim(binaryRegulonActivity)
# 353 1386
binaryRegulonActivity[1:5,1:5]

# Only non-duplicated regulons # JDS: does not seem necessary 
sum(duplicated(rownames(binaryRegulonActivity)))
binaryRegulonActivity_nonDupl <- binaryRegulonActivity[which(rownames(binaryRegulonActivity) %in% onlyNonDuplicatedExtended(rownames(binaryRegulonActivity))),]
dim(binaryRegulonActivity_nonDupl)
# 353 1386

# Get minimum cells
minCells <- ncol(binaryRegulonActivity) * .01 # 1% cells

### All regulons ###
regulonSelection <- list()
regulonSelection[["All_regulons_with_duplicated"]] <- rownames(binaryRegulonActivity)
length(rownames(binaryRegulonActivity))
#353

### Active in > 1% cells ###
regMinCells <- names(which(rowSums(binaryRegulonActivity_nonDupl) > minCells))
regulonSelection[["Regs_active_morethan_1%_cells"]] <- regMinCells
length(regMinCells)
#337

### Correlation between regulons ###
reguCor <- cor(t(binaryRegulonActivity_nonDupl[regMinCells,]))
diag(reguCor) <- 0
nrow(reguCor)
#337

### Regulons correlated with other regulons and active in > 1% cells ###
corrRegs <- names(which(rowSums(abs(reguCor) > 0.30) > 0))
regulonSelection[["Regs_active_and_corr_0.3"]]  <- corrRegs
length(corrRegs)
#132

### Regulons NOT correlated with other regulons and/or active in < 1% cells ###
missingRegs <- rownames(binaryRegulonActivity_nonDupl)[which(!rownames(binaryRegulonActivity_nonDupl) %in% corrRegs)]
regulonSelection[["Regs_inactive_or_noncorr"]]  <- missingRegs
length(missingRegs)
#221
#353-221=132

### Write to Excel
names(regulonSelection)<-paste0(1:4,"_",names(regulonSelection))
# write.xlsx(regulonSelection, file = paste0(sampleFolder,"results/listsRegulons.xlsx"))


## Code Jordy Heatmaps ## 
idents <- sort(seuratObj@active.ident)
ordered <- names(idents)
cellInfo <- as.data.frame(idents)

nrClust <- length(levels(seuratObj@active.ident))
nrClust
#1
colors <- c("#ED8141")
colors <- colorRampPalette(colors)(nrClust)
names(colors) <- unique(cellInfo$idents)
colors <- list(idents = colors)

tables <- table(seuratObj@active.ident)
names(tables) <- NULL

### Setting gaps ###
gaps <- tables[1]
for (k in 1:nrClust) { #add values together
  value <- tables[k]
  gaps <- append(gaps, value+gaps[k-1])
}

for (i in seq_len(length(regulonSelection)))
{
  selRegs <- names(regulonSelection)[i]
  if(length(regulonSelection[[selRegs]])>1)
  {
    AUCMat <- AUCdata[regulonSelection[[selRegs]],,drop=FALSE]
    AUCMat <- AUCMat[,ordered]
    
    avg_exp_df <- data.frame(row.names = regulonSelection[[selRegs]]) #Create empty dataframe
    
    for (j in levels(seuratObj@active.ident)) {
      AUCMatTmp <- AUCMat[,names(seuratObj@active.ident[which(seuratObj@active.ident==j)])] #Subset for celpopulation
      AUCMatTmp <- as.matrix(AUCMatTmp)
      
      if (dim(AUCMatTmp)[2] == 1) {
        average_expressions <- AUCMatTmp
      } else {
        average_expressions <- apply(AUCMatTmp, 1, function(x){ #apply across row: average expression in regulon for that celpop
          sum(x)/length(colnames(AUCMatTmp))
        })
      }
      
      average_expressions <- as.numeric(average_expressions)
      avg_exp_df[,j] <- average_expressions
    }
    
    names(avg_exp_df) <- levels(seuratObj@active.ident)
    AUCMatT <- avg_exp_df
    
    ##    p <-pheatmap(AUCMatT, cluster_cols = F, scale = "row", 
    ##                 main = selRegs, treeheight_row = 0, 
    ##                 cellwidth = 10, cellheight = 10, 
    ##                 show_colnames = T, border_color = "black",
    ##                 fontsize_row = 4, fontsize_col = 4, angle_col = 315, silent = T)
    #JDS: WHEN ONLY ONE SAMPLE IS PRESENT THIS LINE HAS TO BE SLIGHTLY CHANGED
    p <-pheatmap(AUCMatT, cluster_cols = F, cluster_rows = F, scale = "none", 
                 main = selRegs, treeheight_row = 0, 
                 cellwidth = 10, cellheight = 10, 
                 show_colnames = T, border_color = "black",
                 fontsize_row = 4, fontsize_col = 4, angle_col = 315, silent = T)
    
    p.height <- 1 + 0.15 * length(rownames(AUCMatT))
    
    ggsave(filename = paste0(sampleFolder,"results/Heatmap_AUCRegulons_Averaged_",i,".pdf"), p, dpi = 600, height = p.height, width = 10, limitsize = F)
    
  }
}


### BINARY HEATMAPS for different selections
for(i in 1:4)
{
  print(paste0("Print heatmap ",i))
  selRegs <- names(regulonSelection)[i]
  if(length(regulonSelection[[selRegs]])>1)
  {
    binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]],,drop=FALSE]
    NMF::aheatmap(binaryMat, scale="none", revC=TRUE, main=selRegs,
                  annCol=cellInfo[colnames(binaryMat),, drop=FALSE],
                  annColor=NULL,
                  color = c("white", "black"),
                  filename=paste0(sampleFolder,"results/binaryRegulonActivity_Heatmap_",i,"_new.pdf"))
  }
}

### COLORED AUC HEATMAPS for different selections

library(RColorBrewer)
Colorset<-c(brewer.pal(n = 6, name = "Set3")[1],"Yellow",brewer.pal(n = 6, name = "Set3")[c(3:6)])

for(i in seq_len(length(regulonSelection)))
{
  print(paste0("Print heatmap ",i))
  selRegs <- names(regulonSelection)[i]
  
  names(Colorset)<-levels(cellInfo[colnames(AUCMat),, drop=FALSE]$idents) ##Works?????
  Colorset_list <- list(idents = Colorset)
  
  if(length(regulonSelection[[selRegs]])>1)
  {
    AUCMat <- AUCdata[regulonSelection[[selRegs]],,drop=FALSE]
    NMF::aheatmap(AUCMat, scale="none", revC=TRUE, main=selRegs,
                  annCol=cellInfo[colnames(AUCMat),, drop=FALSE],
                  annColor=Colorset_list, 
                  filename=paste0(sampleFolder,"results/RegulonActivity_Heatmap_",i,"_new.pdf"))
  }
}


# Boxplots #cells/regulon & #regulons/cell
pdf(file=paste0(sampleFolder,"results/boxplot_regulonsPerCell.pdf"))
par(mfrow=c(1,2))
boxplot(rowSums(binaryRegulonActivity_nonDupl), main="nCells per regulon",
        sub='number of cells \nthat have the regulon active',
        col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)
boxplot(colSums(binaryRegulonActivity_nonDupl), main="nRegulons per Cell",
        sub='number of regulons \nactive per cell',
        col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)
dev.off()


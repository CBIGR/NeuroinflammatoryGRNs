library(Seurat)
library(SeuratObject)
library(reticulate)
library(Matrix)
library(stringr)

setwd("/kyukon/scratch/gent/vo/000/gvo00027/projects/CBIGR/21HPP_GRN_neuroinfl/singlecell")

scipy <- import("scipy")

dir_name <- paste0("Female_samples_MDD_pyimport/")
dir.create(dir_name)

for (f in list.files(path = "./Female_samples_MDD")) {
  
  if (endsWith(f, ".rds") & startsWith(f, "Seurat_")) {
    message("converting ", f)
    
    f_name   <- str_split_1(f,".rds")[1]
    dir_name <- paste0("Female_samples_MDD_pyimport/",f_name, "_pyimport/")
    dir.create(dir_name)

    so <- readRDS(paste0("Female_samples_MDD/",f))
    mtx <- GetAssayData(so, assay = "RNA", layer = "data")

    scipy$sparse$save_npz(file=paste0(dir_name, f_name, ".npz" ), matrix=mtx)
    write(colnames(mtx), file= paste0(dir_name, f_name, "_cellnames.txt" ))
    write(row.names(mtx), file= paste0(dir_name, f_name, "_featurenames.txt" ))
    # saveRDS(object = mtx, file = paste0( str_split_1(f,".rds")[1], "_matrixonly.rds" ))
  }

}



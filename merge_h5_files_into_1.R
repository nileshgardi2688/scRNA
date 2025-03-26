install.packages("hdf5r")
library(hdf5r)
library(Seurat)
h5_files <- list.files(pattern = "*.h5")
h5_read <- lapply(h5_files, Read10X_h5)
h5_seurat <- lapply(h5_read, CreateSeuratObject)
h5_seurat  <- merge(h5_seurat[[1]], y = h5_seurat[2:length(h5_seurat)], 
                  add.cell.ids = c("1", "2", "3"), project = "project")
				  #Add more cell ids for each of your eight samples.
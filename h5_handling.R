library(rhdf5)
library(Seurat)

geo <- H5Fopen("~/Downloads/GSE190889_annoted_fullsize.h5")

h5ls(geo)

mat <- as.sparse(geo$matrix$counts)

barcodes <- geo$matrix$barcodes
features <- geo$matrix$genes

rownames(mat) <- features
colnames(mat) <- barcodes


seu <- CreateSeuratObject(counts = mat)

seu[["cell_type"]] <- geo$metadata$cell_type
seu[["health_status"]] <- geo$metadata$health_status
seu[["index"]] <- geo$metadata$index
seu[["internal_id"]] <- geo$metadata$internal_id
seu[["n_counts"]] <- geo$metadata$n_counts
seu[["n_genes"]] <- geo$metadata$n_genes
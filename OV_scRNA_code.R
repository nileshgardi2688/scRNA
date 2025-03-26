# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Set seed
set.seed(42)

# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

#####Unbiased
setwd("/home/pavan/Nilesh/single_cell_data/acral/")

# get data location
dirs <- list.dirs(path = 'Data/', recursive = F, full.names = F)

for(x in dirs){
    #name <- gsub('_filtered_feature_bc_matrix','', x)
    
    cts <- ReadMtx(mtx = paste0('Data/',x,'/matrix.mtx'),
                   features = paste0('Data/',x,'/genes.tsv'),
                   cells = paste0('Data/',x,'/barcodes.tsv'))
    
    # create seurat objects
    assign(x, CreateSeuratObject(counts = cts))
}


# merge datasets
ls()
merged_seurat <- merge(HRR581078, y = c(HRR581079, HRR581080,HRR581081,HRR581082,HRR581083,HRR581084,HRR581085),
      add.cell.ids = ls()[3:10],
      project = 'ACRAL')
	  
merged_seurat

#Merging Seurat Objects: Multiple Counts in Assays Layer
merged_seurat <- JoinLayers(merged_seurat)

#################################QC
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')
merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 500 & mitoPercent < 25)

merged_seurat_filtered

# create a sample column
a<-read.table("Patients.txt",header=T)
merged_seurat_filtered$Patient <- a[,2]


View(merged_seurat_filtered@meta.data)

############### split the RNA measurements into patient sepecific manner
merged_seurat_filtered[["RNA"]] <- split(merged_seurat_filtered[["RNA"]], f = merged_seurat_filtered$Patient)


# pre-process standard workflow
merged_seurat_filtered <- NormalizeData(merged_seurat_filtered)
merged_seurat_filtered<- FindVariableFeatures(merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(merged_seurat_filtered,resolution = c(0.1,0.3, 0.5))
#View(merged_seurat_filtered@meta.data)
Idents(merged_seurat_filtered) <- "RNA_snn_res.0.1"
merged_seurat_filtered<- RunUMAP(object = merged_seurat_filtered, dims = 1:20,reduction = "pca", reduction.name = "umap.unintegrated")
table(merged_seurat_filtered@meta.data$RNA_snn_res.0.1)


DimPlot(merged_seurat_filtered, reduction = "umap.unintegrated",group.by ="Patient")
write.table(merged_seurat_filtered@meta.data,"merged_seurat_filtered_metadata.txt",sep="\t")


# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap.unintegrated', group.by = 'Patient', label=TRUE)

#############################################Integration
merged_seurat_filtered <- IntegrateLayers(object = merged_seurat_filtered, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",verbose = FALSE)

# re-join layers after integration
merged_seurat_filtered[["RNA"]] <- JoinLayers(merged_seurat_filtered[["RNA"]])

merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, reduction = "integrated.cca", dims = 1:20)
merged_seurat_filtered <- FindClusters(merged_seurat_filtered, resolution = c(0.1,0.3, 0.5, 0.8))

Idents(merged_seurat_filtered) <- "RNA_snn_res.0.1"
merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:20, reduction = "integrated.cca")

# Visualization
p2<-DimPlot(merged_seurat_filtered, reduction = "umap", group.by = "Patient",label=TRUE)

##################Before and after integration
png("dimplot_before_after_integration_patient.png")
grid.arrange(p1, p2)
dev.off()

png("dimplot_res_0.1.png")
DimPlot(merged_seurat_filtered, reduction = 'umap',label=TRUE)
dev.off()


#######################################Manual annotation of clusters
Idents(merged_seurat_filtered) <- "RNA_snn_res.0.1"
markers.to.plot <- c("MITF","COL1A1","CD3D","C1QB","PECAM1","S100A8","CD79A","EPCAM","KRT14","DLL3","TYRP1")
png("Dotplot_clusterd_known_markers.png")
DotPlot(merged_seurat_filtered, features = markers.to.plot,dot.scale = 8)+ RotatedAxis()
dev.off()


####################################################Violin plots
#Malignant
png("Malignant_MITF.png")
VlnPlot(merged_seurat_filtered, features = c("MITF"),layer = "data", log = TRUE)
dev.off()
#CAF
png("CAF_COL1A1.png")
VlnPlot(merged_seurat_filtered, features = c("COL1A1"),layer = "data", log = TRUE)
dev.off()
#T
png("T_cell_CD3D.png")
VlnPlot(merged_seurat_filtered, features = c("CD3D"),layer = "data", log = TRUE)
dev.off()
#Macrophage
png("Macrophage_C1QB.png")
VlnPlot(merged_seurat_filtered, features = c("C1QB"),layer = "data", log = TRUE)
dev.off()
#Endothelial
png("Endothelial_PECAM1.png")
VlnPlot(merged_seurat_filtered, features = c("PECAM1"),layer = "data", log = TRUE)
dev.off()
#Neutriphils
png("Neutriphils_S100A8.png")
VlnPlot(merged_seurat_filtered, features = c("S100A8"),layer = "data", log = TRUE)
dev.off()
#B
png("B_cell_CD79A.png")
VlnPlot(merged_seurat_filtered, features = c("CD79A"),layer = "data", log = TRUE)
dev.off()
#Epithelial
png("Epithelial_EPCAM.png")
VlnPlot(merged_seurat_filtered, features = c("EPCAM"),layer = "data", log = TRUE)
dev.off()
#Keratinocytes
png("Keratinocytes_KRT14.png")
VlnPlot(merged_seurat_filtered, features = c("KRT14"),layer = "data",log = TRUE)
dev.off()
#DLL3
png("DLL3.png")
VlnPlot(merged_seurat_filtered, features = c("DLL3"),layer = "data",log = TRUE)
dev.off()
#TYRP1
png("TYRP1.png")
VlnPlot(merged_seurat_filtered, features = c("TYRP1"),layer = "data",log = TRUE)
dev.off()




#######################################Feature plots
#Malignant
png("Malignant_MITF_featureplot.png")
FeaturePlot(merged_seurat_filtered, features = c("MITF"), reduction = "umap")
dev.off()
#CAF
png("CAF_COL1A1_featureplot.png")
FeaturePlot(merged_seurat_filtered, features = c("COL1A1"),reduction = "umap")
dev.off()
#T
png("T_cell_CD3D_featureplot.png")
FeaturePlot(merged_seurat_filtered, features = c("CD3D"),reduction = "umap")
dev.off()
#Macrophage
png("Macrophage_C1QB_featureplot.png")
FeaturePlot(merged_seurat_filtered, features = c("C1QB"),reduction = "umap")
dev.off()
#Endothelial
png("Endothelial_PECAM1_featureplot.png")
FeaturePlot(merged_seurat_filtered, features = c("PECAM1"),reduction = "umap")
dev.off()
#Neutriphils
png("Neutriphils_S100A8_featureplot.png")
FeaturePlot(merged_seurat_filtered, features = c("S100A8"),reduction = "umap")
dev.off()
#B
png("B_cell_CD79A_featureplot.png")
FeaturePlot(merged_seurat_filtered, features = c("CD79A"),reduction = "umap")
dev.off()
#Epithelial
png("Epithelial_EPCAM_featureplot.png")
FeaturePlot(merged_seurat_filtered, features = c("EPCAM"),reduction = "umap")
dev.off()
#Keratinocytes
png("Keratinocytes_KRT14_featureplot.png")
FeaturePlot(merged_seurat_filtered, features = c("KRT14"),reduction = "umap")
dev.off()
png("DLL3_featureplot.png")
FeaturePlot(merged_seurat_filtered, features = c("DLL3"),reduction = "umap")
dev.off()
png("TYRP1_featureplot.png")
FeaturePlot(merged_seurat_filtered, features = c("TYRP1"),reduction = "umap")
dev.off()

write.table(FetchData(object = merged_seurat_filtered, vars = c("MITF","COL1A1","CD3D","C1QB","PECAM1","S100A8","CD79A","EPCAM","KRT14","DLL3","TYRP1"), layer = "data"),"After_integration_All_known_markers_DLL3_TYRP1_data.txt",sep="\t")


#####################################Annotation with singleR (on my own windows PC)

library(celldex)
library(SingleR)
# Access the raw counts
merged_seurat_filtered_counts <- GetAssayData(merged_seurat_filtered, layer = 'counts')

# Access the normalised data
merged_seurat_filtered_data <- GetAssayData(merged_seurat_filtered, layer= 'data')

# 2. Get reference dataset
ref <- celldex::HumanPrimaryCellAtlasData()
unique(ref$label.main)
ref <- ref[,grepl('DC|B_cell|Neutrophils|Epithelial_cells|T_cells|Endothelial_cells|Keratinocytes|Monocyte|Erythroblast|
+                  Macrophage|NK_cell|Fibroblast', ref$label.main)]


# 3. Run SingleR

ct_ann <- SingleR(test = merged_seurat_filtered_data, # we could also use sce or raw_counts
                  ref = ref, 
                  labels = ref$label.main,
                  de.method = 'wilcox')

merged_seurat_filtered$singleR.labels <- ct_ann$pruned.labels[match(rownames(merged_seurat_filtered@meta.data), rownames(ct_ann))]
DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'singleR.labels')


write.table(merged_seurat_filtered@meta.data,"merged_seurat_filtered_metadata_singleR.txt",sep="\t")


a<-(janitor::tabyl(merged_seurat_filtered@meta.data,RNA_snn_res.0.1,singleR.labels))
write.table(a,"overlap_clusters_singleR.txt",sep="\t")




##################################Not used
####################remove unwanted objects from enviornment
rm(HRR581078, HRR581079, HRR581080,HRR581081,HRR581082,HRR581083,HRR581084,HRR581085)
# perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Patient')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                       anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)

# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- FindNeighbors(object = seurat.integrated, dims = 1:20)
seurat.integrated <- FindClusters(object = seurat.integrated,resolution = c(0.1,0.3, 0.5))
View(seurat.integrated@meta.data)
write.table(seurat.integrated@meta.data,"seurat.integrated_metadata.txt",sep="\t")
Idents(seurat.integrated) <- "integrated_snn_res.0.1"
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:20)
table(seurat.integrated@meta.data$RNA_snn_res.0.1)
#Idents(seurat.integrated) <- "integrated_snn_res.0.3"
#seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:20)


p2 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Patient')

png("dimplot_before_after_integration_patient.png")
grid.arrange(p1, p2)
dev.off()

png("dimplot_res_0.1.png")
DimPlot(seurat.integrated, reduction = 'umap',label=TRUE)
dev.off()

#Malignant
VlnPlot(merged_seurat_filtered, features = c("MITF"),layer = "data", log = TRUE)
#CAF
VlnPlot(merged_seurat_filtered, features = c("COL1A1"),layer = "data", log = TRUE)
#T
VlnPlot(merged_seurat_filtered, features = c("CD3D"),layer = "data", log = TRUE)
#Macrophage
VlnPlot(merged_seurat_filtered, features = c("C1QB"),layer = "data", log = TRUE)
#Endothelial
VlnPlot(merged_seurat_filtered, features = c("PECAM1"),layer = "data", log = TRUE)
#Neutriphils
VlnPlot(merged_seurat_filtered, features = c("S100A8"),layer = "data", log = TRUE)
#B
VlnPlot(merged_seurat_filtered, features = c("CD79A"),layer = "data", log = TRUE)
#Epithelial
VlnPlot(merged_seurat_filtered, features = c("EPCAM"),layer = "data", log = TRUE)
#Keratinocytes
VlnPlot(merged_seurat_filtered, features = c("KRT14"),layer = "data", log = TRUE)


############################################################copykat usage
setwd("/home/pavan/Nilesh/single_cell_data/acral/Data")

HRR581078 <- ReadMtx(mtx = paste0('HRR581078/','/matrix.mtx'),
                   features = paste0('HRR581078/','/genes.tsv'),
                   cells = paste0('HRR581078/','/barcodes.tsv'))
				   
#HRR581078[1:10,1:10]

# create Seurat object
HRR581078.seurat <- CreateSeuratObject(counts = HRR581078)

#QC
HRR581078.seurat$mitoPercent <- PercentageFeatureSet(HRR581078.seurat, pattern='^MT-')
HRR581078.seurat_filtered <- subset(HRR581078.seurat, subset = nFeature_RNA > 500 & mitoPercent < 25)

#Matrix
HRR581078_rawmat <- as.matrix(HRR581078.seurat_filtered[["RNA"]]$counts)

write.table(HRR581078_rawmat,"HRR581078.rawmat.txt",sep="\t", quote = FALSE, row.names = TRUE)




plots <- VlnPlot(ifnb, features = c("CD3D"), group.by = "integrated_snn_res.0.1",pt.size = 0, combine = FALSE)











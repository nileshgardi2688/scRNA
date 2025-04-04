
```{r setup}

### load libraries
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)

theme_set(theme_cowplot())

#color scheme
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  Patient1 = "#E2D200",
  Patient2 = "#FD6467",
  Patient3 = "#B40F20",
  Patient4 = "royalblue4",
  Patient5 = "dodgerblue3",
  Patient6 = "deepskyblue",
  Patient7 = "cadetblue1",
  Patient8 = "darkolivegreen2",
  Patient9 = "chartreuse3",
  Patient10 = "darkgreen")

#filter parameters
nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 40
pHB_lower <- 0
pHB_upper <- 5

```

```{r data loading}

### sample list
samples <- read_excel("patients_metadata.xlsx") %>% .$sample_id_old

### import cellranger files from different data sets

for (i in seq_along(samples)){
  assign(paste0("scs_data", i), Read10X_h5(filename = paste0(samples[i], ".h5")))
}

### create seurat objects from cellranger files
for (i in seq_along(samples)){
  assign(paste0("seu_obj", i), CreateSeuratObject(counts = eval(parse(text = paste0("scs_data", i))), project = samples[i], min.cells = 3))
}

### merge data sets
seu_obj <- merge(seu_obj1, y = c(seu_obj2, seu_obj3, seu_obj4, seu_obj5, seu_obj6, seu_obj7, seu_obj8, seu_obj9, seu_obj10, seu_obj11, seu_obj12, seu_obj13), add.cell.ids = samples, project = "carcinoid")

### calculate mitochondrial, hemoglobin and ribosomal gene counts
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^MT-", col.name = "pMT")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^HBA|^HBB", col.name = "pHB")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^RPS|^RPL", col.name = "pRP")

```

```{r clear environment}

remove(seu_obj1)
remove(seu_obj2)
remove(seu_obj3)
remove(seu_obj4)
remove(seu_obj5)
remove(seu_obj6)
remove(seu_obj7)
remove(seu_obj8)
remove(seu_obj9)
remove(seu_obj10)
remove(seu_obj11)
remove(seu_obj12)
remove(seu_obj13)

remove(scs_data1)
remove(scs_data2)
remove(scs_data3)
remove(scs_data4)
remove(scs_data5)
remove(scs_data6)
remove(scs_data7)
remove(scs_data8)
remove(scs_data9)
remove(scs_data10)
remove(scs_data11)
remove(scs_data12)
remove(scs_data13)

```

```{r add clinical data}

metatable <- read_excel("patients_metadata.xlsx")

metadata <- FetchData(seu_obj, "orig.ident")
metadata$cell_id <- rownames(metadata)
metadata$sample_id_old <- metadata$orig.ident
metadata <- left_join(x = metadata, y = metatable, by = "sample_id_old")
rownames(metadata) <- metadata$cell_id

seu_obj <- AddMetaData(seu_obj, metadata = metadata)

```


```{r plot wrapper function}

qc_std_plot_helper <- function(x) x + 
  scale_color_viridis() +
  geom_point(size = 0.01, alpha = 0.3)

qc_std_plot <- function(seu_obj) {
  qc_data <- as_tibble(FetchData(seu_obj, c("nCount_RNA", "nFeature_RNA", "pMT", "pHB", "pRP")))
  plot_grid(
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pMT))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pHB))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pRP))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pMT, color = nFeature_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pHB, color = nFeature_RNA))) + 
      geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pRP, color = nFeature_RNA))) + 
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pMT, color = nCount_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pHB, color = nCount_RNA))) + 
      geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pRP, color = nCount_RNA))) + 
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    
    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nCount_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nFeature_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2),
    
    
    ggplot(gather(qc_data, key, value), aes(key, value)) +
      geom_violin() +
      facet_wrap(~key, scales = "free", ncol = 5),
    
    ncol = 3, align = "hv"
  )
}

```

```{r qc plots before filtering}

seu_obj_unfiltered <- seu_obj

saveRDS(seu_obj_unfiltered, file = "seurat_objects/all_unfiltered.RDS")



qc_std_plot(seu_obj_unfiltered)
ggsave2("qc_scatterplot_all.png", path = "output/qc", width = 30, height = 30, units = "cm")



qcparams <- c("nFeature_RNA", "nCount_RNA", "pMT", "pHB", "pRP")
for (i in seq_along(qcparams)){
  print(VlnPlot(object = seu_obj, features = qcparams[i], group.by = "sample_id_old", pt.size = 0,layer = "counts", log = TRUE))
}
for (i in seq_along(qcparams)){
  print(RidgePlot(object = seu_obj, features = qcparams[i], group.by = "sample_id_old",layer = "counts", log = TRUE))
}

VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "pMT"), pt.size = 0, group.by = "sample_id_old", ncol = 1,layer = "counts", log = TRUE)>ggsave2("qc_vlnplots_nfeature_ncount_pmt_per_patient.pdf", path = "output/qc", width = 30, height = 20, units = "cm")

```

```{r after filtering}

#seu_obj_unfiltered <- readRDS("seurat_objects/all_unfiltered.RDS")

seu_obj <- subset(seu_obj_unfiltered, subset = 
                    nFeature_RNA > nFeature_lower & 
                    nFeature_RNA < nFeature_upper & 
                    nCount_RNA > nCount_lower & 
                    nCount_RNA < nCount_upper & 
                    pMT < pMT_upper & 
                    pHB < pHB_upper)

qc_std_plot(seu_obj)

seu_obj_unfiltered
seu_obj

```

```{r normalization}
options(future.globals.maxSize = 1000 * 1024^2)
seu_obj <- SCTransform(seu_obj, verbose = T, vars.to.regress = c("nCount_RNA", "pMT"), conserve.memory = T)

saveRDS(seu_obj, file = "seurat_objects/all_SCTransform.RDS")

```

```{r dimensionality reduction and clustering}

seu_obj <- RunPCA(seu_obj)
ElbowPlot(seu_obj, ndims = 50)

#for (i in c(5, 10, 15, 20)) {
#  umaptest <- RunUMAP(seu_obj, dims = 1:i, verbose = F)
#  print(DimPlot(umaptest, reduction = "umap", group.by = "sample_id") + labs(title = paste0(i, " dimensions")))
#  print(FeaturePlot(umaptest, features = c("EPCAM", "PTPRC"), sort.cell = T))
#  print(FeaturePlot(umaptest, features = c("MARCO", "KIT"), sort.cell = T))
#  print(FeaturePlot(umaptest, features = c("FOXJ1", "AGER"), sort.cell = T))
#  print(FeaturePlot(umaptest, features = c("JCHAIN", "VWF"), sort.cell = T))
#  remove(umaptest)
#}

seu_obj <- RunUMAP(seu_obj, dims = 1:10, verbose = T)
seu_obj <- FindNeighbors(seu_obj, dims = 1:10)

for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  seu_obj <- FindClusters(seu_obj, resolution = i)
  print(DimPlot(seu_obj, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}

for (i in c("nFeature_RNA", "nCount_RNA", "pMT", "pHB", "pRP")) {
  print(FeaturePlot(seu_obj, features = i, coord.fixed = T, order = T))
}

DimPlot(seu_obj, group.by = "sample_id_old")
ggsave2("DimPlot_per_patient.pdf", path = "output/qc", width = 30, height = 20, units = "cm")
DimPlot(seu_obj, group.by = "SCT_snn_res.0.2", label = T)


```

```{r main cell type annotation}

mainmarkers <- c("PECAM1", "VWF", "ACTA2", "JCHAIN", "MS4A1", "PTPRC", "CD68", "KIT", "EPCAM", "CDH1", "KRT7", "KRT19")

for (i in seq_along(mainmarkers)) {
  FeaturePlot(seu_obj, features = mainmarkers[i], coord.fixed = T, order = T, cols = viridis(10))
  ggsave2(paste0("FeaturePlot_mainmarkers_", mainmarkers[i], ".png"), path = "output/annotation", width = 10, height = 10, units = "cm")
}

DotPlot(seu_obj, features = mainmarkers, group.by = "SCT_snn_res.0.2") +   coord_flip() +   scale_color_viridis()
ggsave2("DotPlot_mainmarkers.png", path = "output/annotation", width = 30, height = 8, units = "cm")

DimPlot(seu_obj, group.by = "SCT_snn_res.0.2", label = T, label.size = 5)
ggsave2("DimPlot_all_clusters.png", path = "output/annotation", width = 20, height = 20, units = "cm")

Idents(seu_obj) <- seu_obj$SCT_snn_res.0.2
annotation_curated_main <- read_excel("curated_annotation/curated_annotation_main.xlsx")
new_ids_main <- annotation_curated_main$main_cell_type
names(new_ids_main) <- levels(seu_obj)
seu_obj <- RenameIdents(seu_obj, new_ids_main)
seu_obj@meta.data$main_cell_type <- Idents(seu_obj)

```

```{r cell cycle scoring}

### add cell cycle, cc.genes loaded with Seurat

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

score_cc <- function(seu_obj) {
  seu_obj <- CellCycleScoring(seu_obj, s.genes, g2m.genes)
  seu_obj@meta.data$CC.Diff <- seu_obj@meta.data$S.Score - seu_obj@meta.data$G2M.Score
  return(seu_obj)
}

seu_obj <- score_cc(seu_obj)

FeatureScatter(seu_obj, "G2M.Score", "S.Score", group.by = "Phase", pt.size = .1) +
  coord_fixed(ratio = 1)

```

```{r subset, rescale and save RDS files}

saveRDS(seu_obj, file = "seurat_objects/all.RDS")

Idents(seu_obj) <- seu_obj@meta.data$main_cell_type

epi <- subset(seu_obj, idents = "Epithelial")
imm <- subset(seu_obj, idents = "Immune")
str <- subset(seu_obj, idents = "Stromal")

epi <- ScaleData(epi)
imm <- ScaleData(imm)
str <- ScaleData(str)

saveRDS(epi, file = "seurat_objects/epi.RDS")
saveRDS(imm, file = "seurat_objects/imm.RDS")
saveRDS(str, file = "seurat_objects/str.RDS")

```

```{r plots for figure 1}
DimPlot(seu_obj, group.by = "tissue_type", cols = use_colors, pt.size = 0.1)
ggsave2("DimPlot_Tumor_Normal.png", path = "output/fig1", width = 15, height = 15, units = "cm")

DimPlot(seu_obj, group.by = "patient_id", cols = use_colors, pt.size = 0.1)
ggsave2("DimPlot_Patients.png", path = "output/fig1", width = 15, height = 15, units = "cm")

DimPlot(seu_obj, group.by = "main_cell_type", cols = use_colors, pt.size = 0.1)
ggsave2("DimPlot_Main_cell_type.png", path = "output/fig1", width = 15, height = 15, units = "cm")

cell_types <- FetchData(seu_obj, vars = c("sample_id_old", "main_cell_type", "tissue_type")) %>% 
  mutate(main_cell_type = factor(main_cell_type, levels = c("Stromal", "Immune", "Epithelial")))

ggplot(data = cell_types) + 
  geom_bar(mapping = aes(x = sample_id_old, fill = main_cell_type, ), position = "fill", width = 0.75) +
  scale_fill_manual(values = use_colors) +
  coord_flip()
ggsave2("main_cell_type_count_relative.pdf", path = "output/fig1", width = 15, height = 30, units = "cm")

ggplot(data = cell_types) + 
  geom_bar(mapping = aes(x = sample_id_old, fill = main_cell_type, ), width = 0.75) +
  scale_fill_manual(values = use_colors) +
  coord_flip()
ggsave2("main_cell_type_count_absolute.pdf", path = "output/fig1", width = 15, height = 30, units = "cm")

```

write.table(FetchData(object = seu_obj, vars = c("DLL3","SEZ6","EPCAM","CDH17","SSTR2","CD276","TACSTD2","NECTIN4"), layer = "data"),"DLL3_others.txt",sep="\t")
write.table(seu_obj@meta.data,"metadata.txt",sep="\t")

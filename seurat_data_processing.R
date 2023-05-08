library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(dplyr)

setwd("GSE184878_RAW_SZ/")  # path to the directory where data is saved

dirs <- list.files(path = "count_matrix/", recursive = F, full.names = F)

# go in each file to read the count matrix and create Seurat Object

for(x in dirs){
  name <- gsub('_10xfiltered_feature_count_matrix.csv','', x)
  
  # read matrix, features, barcodes
  cts <- read.csv(file = paste0("count_matrix/",x))
  
  # create seurat objects from the count matrix, assign it to value of name
  assign(name, CreateSeuratObject(counts = cts))
}


ls()

merged.seurat <- merge(GSM5599210_C1, y=c(GSM5599211_C2, GSM5599212_C3, GSM5599213_C4, GSM5599214_S1, GSM5599215_S2,
                                          GSM5599216_S3), add.cell.ids = ls()[3:9], project = "SZ_markers")

merged.seurat   # 33538 features across 27793 samples within 1 assay
str(merged.seurat)


# Quality Control and Filtering -------------------------
View(merged.seurat@meta.data)

merged.seurat@meta.data <- merged.seurat@meta.data[-1,]  # first row is NULL


# to get patient_id and tissue name in a separate column - extract from rownames
#create a sample column, split the sample column
merged.seurat$sample <- rownames(merged.seurat@meta.data)
merged.seurat@meta.data <- separate(merged.seurat@meta.data, col = 'sample', into = c('patient','condition','barcode'), 
                                    sep = '_')



# to specifiy control and schizophrenia condition

for (i in merged.seurat@meta.data$condition) {
  cond_substr <- substr(i,1,1)
  if(cond_substr=="C"){
    merged.seurat@meta.data$condition[which(merged.seurat@meta.data$condition==i)] <- 'C'
  }
  else if(cond_substr=="S"){
    merged.seurat@meta.data$condition[which(merged.seurat@meta.data$condition==i)] <- 'S'
  }
}


unique(merged.seurat@meta.data$patient) # 7 patients
unique(merged.seurat@meta.data$condition) # 2 conditions

# quality control - check percentage of mitochondrial cells 
merged.seurat$mitoPercent<- PercentageFeatureSet(merged.seurat, pattern = '^MT-') #/or merged_seurat[['mitoPercent']]


# visualize features as violin plot
colnames(merged.seurat@meta.data) # check column names of metadata

# to visualize good quality cells - good number of genes and molecules together
FeatureScatter(merged.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ geom_smooth(method = 'lm')
#good quality cells should follow a straight line


# filtering merged data -------------------------
merged.seurat_filtered <- subset(merged.seurat, subset = nCount_RNA > 800 & nFeature_RNA > 500 & mitoPercent < 10)

merged.seurat_filtered  # 33538 features across 25275 samples within 1 assay


# perform standard workflow steps to find if there is any batch effects  -----------------------------
merged.seurat_filtered <- FindVariableFeatures(object = merged.seurat_filtered, selection.method = "vst", nfeatures = 2000)  # HVG - top 2000

merged.seurat_filtered <- NormalizeData(object = merged.seurat_filtered)  # log normalization

merged.seurat_filtered <- ScaleData(object = merged.seurat_filtered)  # scaling

merged.seurat_filtered <- RunPCA(object = merged.seurat_filtered)  # dimensionality reduction PCA

ElbowPlot(merged.seurat_filtered)  

merged.seurat_filtered <- FindNeighbors(object = merged.seurat_filtered) # find nearest neighbours

# check with a series of different resolution to find best one
merged.seurat_filtered <- FindClusters(object = merged.seurat_filtered, resolution = c(0.1,0.2,0.3,0.5,0.7))

View(merged.seurat_filtered@meta.data)

DimPlot(merged.seurat_filtered, group.by = "RNA_snn_res.0.1", label = TRUE, label.size = 5) # 7 clusters
DimPlot(merged.seurat_filtered, group.by = "RNA_snn_res.0.2", label = TRUE, label.size = 5) # 8 clusters


# Setting identity of clusters -----------
Idents(merged.seurat_filtered) # default 21 clusters - not required
Idents(merged.seurat_filtered) <- "RNA_snn_res.0.1"  # change it to 7 clusters of resolution 0.1
Idents(merged.seurat_filtered)


merged.seurat_filtered <- RunUMAP(object = merged.seurat_filtered, dims = 1:50) # non-linear dimensionality reduction

DimPlot(merged.seurat_filtered, reduction = "umap", label = TRUE) # 7 clusters

p1 <- DimPlot(merged.seurat_filtered, reduction = "umap", group.by = 'patient', label = TRUE, label.size = 5)  # group by patients - 7
p2 <- DimPlot(merged.seurat_filtered, reduction = "umap", group.by = 'condition', label = TRUE, label.size = 5)  # group by condition type - 2

grid.arrange(p1, p2, ncol = 2)

# batch effect found by condition - perform integration on condition ----------------

obj.list<- SplitObject(merged.seurat_filtered, split.by = 'condition')  # if batch effects coming from condition

obj.list  # list of objects split according to 2 conditions - S, C

# for each object do normalization and identify HVGs
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

#find integration anchors (canonical correlation analysis - CCA) - to integrate the data across different conditions
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

# Integrate the data using the anchors
seurat.integrated <- IntegrateData(anchorset = anchors)

# perform remaining workflow - scale, PCA, UMAP on new object seurat.integrated --------------------------
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
ElbowPlot(seurat.integrated)

seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)  # choosing 50 dimensions

#DimPlot(seurat.integrated, reduction = 'pca')
DimPlot(seurat.integrated, reduction = "umap", label = TRUE, label.size = 5)

p11 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'patient', label = TRUE, label.size = 5)
p22 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'condition', label = TRUE, label.size = 5)

grid.arrange(p11, p22, ncol = 2)

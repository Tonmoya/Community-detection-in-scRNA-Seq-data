# continued from seurat_data_processing.R


library(WGCNA)
library(purrr)
library(tibble)
library(Seurat)

levels(seurat.integrated) # 7 clusters

#Find marker genes using - Wilcoxon Rank SUm default, MAST, t-test - available in Seurat

markers.wrs <- FindAllMarkers(seurat.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)  # Wilcoxon Rank SUm default

markers.mast <- FindAllMarkers(seurat.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")  # MAST 

markers.t <- FindAllMarkers(seurat.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "t")  # t-test 

markers.intersect <- intersect(markers.wrs[,7], c(markers.mast[,7],markers.t[,7]))  # intersect the three marker gene lists to get common genes


## DATA IMPUTATION using SAVER -------------------------------------------------

library(SAVER)
library(flashClust)
library(igraph)

genes_mark <- markers.intersect   # taking the common marker genes
beta_mat <- as.data.frame(as.matrix(seurat.integrated@assays[["RNA"]]@data))
sub_mar <- beta_mat[genes_mark,]

# imputation using SAVER to remove zero inflated values
Beta_saver <- saver(sub_mar, ncores = 12, estimates.only = TRUE)    # rows = genes, columns = cells

datExpr <- as.data.frame(t(Beta_saver))  # rows = cells, columns = genes
dim(datExpr)
# 
# # only keep good genes
datExpr <- datExpr[,goodGenes(datExpr)]  # WGCNA function
dim(datExpr) # 25275 cells X 934 genes 


# Divide data into control and SZ
dat_control <- data.frame(matrix(nrow = 1,ncol = dim(datExpr)[2]))
dat_sz <- data.frame(matrix(nrow = 1,ncol = dim(datExpr)[2]))

colnames(dat_control) <- as.numeric(colnames(datExpr))   # because gene names in column are numeric, if string then no need
colnames(dat_sz) <- as.numeric(colnames(datExpr))

for (s in 1:dim(datExpr)[1]) {
  if(substr(rownames(datExpr)[s],12,12)=="C"){
    dat_control <- rbind(dat_control,datExpr[s,])
  }
  else if(substr(rownames(datExpr)[s],12,12)=="S"){
    dat_sz <- rbind(dat_sz,datExpr[s,])
  }
  
}

dat_control <- dat_control[-1,]  # if first row is NA
dat_sz <- dat_sz[-1,]  # if first row is NA

write.csv(dat_control, file = "dat_control.csv")
write.csv(dat_sz, file = "dat_sz.csv")

# Find correlation matrix for control and SZ datasets - gene-gene pearson correlation
corr_control <- cor(dat_control, method = "pearson")  

corr_sz <- cor(dat_sz, method = "pearson")

dim(corr_control)
dim(corr_sz)      
#colSums(corr_genes==0)   # no zero values found in columns


# Using the SZ data for further analysis
corr_genes <- corr_sz

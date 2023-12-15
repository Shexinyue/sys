##加载包
pacman::p_load(hdf5r, dplyr, Seurat, patchwork, SingleR, celldex, pheatmap, ggplot2, cowplot, harmony, stringr)

##导入seurat对象
load("D:/data/GSE167036-8/umap_seurat.RData")

##infercnv的簇信息
library(readr)
PA1_kmeans_df_s <- read_delim("D:/data/inferCNV/PA1/kmeans_df_s.txt", delim = "\t", 
                              escape_double = FALSE, trim_ws = TRUE)
PA2_kmeans_df_s <- read_delim("D:/data/inferCNV/PA2/kmeans_df_s.txt", delim = "\t", 
                              escape_double = FALSE, trim_ws = TRUE)

#PA2-8

PA1_malignant <- PA1_kmeans_df_s[which(PA1_kmeans_df_s$kmeans_class %in% c(2,3,5)),]$CA
PA2_malignant <- PA2_kmeans_df_s[which(PA2_kmeans_df_s$kmeans_class %in% c(3,4)),]$CA
PA3_malignant <- PA3_kmeans_df_s[which(PA3_kmeans_df_s$kmeans_class %in% c(1,7)),]$CA
PA4_malignant <- PA4_kmeans_df_s[which(PA4_kmeans_df_s$kmeans_class %in% c(1,4,5,6)),]$CA
PA5_malignant <- PA5_kmeans_df_s[which(PA5_kmeans_df_s$kmeans_class %in% c(1,2,4,6,7)),]$CA
PA6_malignant <- PA6_kmeans_df_s[which(PA6_kmeans_df_s$kmeans_class %in% c(1,2,5)),]$CA
PA7_malignant <- PA7_kmeans_df_s[which(PA7_kmeans_df_s$kmeans_class %in% c(1,2,6,7)),]$CA
PA8_malignant <- PA8_kmeans_df_s[which(PA8_kmeans_df_s$kmeans_class %in% c(1,2,4,5,6,7)),]$CA


###
malignant <- c(PA1_malignant,PA2_malignant,PA3_malignant,PA4_malignant,PA5_malignant,
                PA6_malignant, PA7_malignant, PA8_malignant)
data_seurat_malignant <- data_seurat[,which(colnames(data_seurat) %in% malignant)]
data_seurat_malignant <- subset(data_seurat_malignant, labels == "Epithelial_cells")

#保存malli_seurat
save(data_seurat_malignant,file = "D:/data/GSE167036-8/data_seurat_malignant.RData")


##
m <- data.frame(
  CB <- colnames(data_seurat),
  mali <- c(rep("Benign",length(colnames(data_seurat))))
)
m$mali[which(m$CA %in% colnames(data_seurat_malignant))] <- "malignant"

data_seurat <- AddMetaData(data_seurat, metadata = m$malli, col.name = "mallig")


#######
load("D:/data/GSE167036-8/data_seurat_mallignant.RData")

data_seurat_malignant <- NormalizeData(data_seurat_malignant, normalization.method = "LogNormalize", scale.factor = 10000)

#
data_seurat_malignant_PT <- data_seurat_malignant[,which(data_seurat_malignant@meta.data$Patient %in% 
                                           c("CA1","CA2","CA3","CA4","CA5","CA6","CA7","CA8"))]
data_seurat_malignant_LNM <- data_seurat_malignant[,which(data_seurat_malignant@meta.data$Patient %in% 
                                                           c("LN1","LN2","LN3","LN4","LN5","LN6","LN7","LN8"))]

#鉴定高变基因
#这一步的目的是鉴定出细胞与细胞之间表达量相差很大的基因，用于后续鉴定细胞类型，
#我们使用默认参数，即“vst”方法选取2000个高变基因。
data_seurat_malignant <- FindVariableFeatures(data_seurat_malignant, selection.method = "vst", nfeatures = 2000)

#数据缩放mean 0 and variance 1
data_seurat_malignant <- ScaleData(data_seurat_malignant)
#线性降维PCA
data_seurat_malignant <- RunPCA(data_seurat_malignant, features = VariableFeatures(object = data_seurat_malignant))
#诊断图

cur_seu <- data_seurat_malignant %>% SCTransform() %>% RunPCA(dims = 1:30)
pct <- cur_seu[["pca"]]@stdev / sum(cur_seu[["pca"]]@stdev)*100
cumu <- cumsum(pct)
pc.use <- min(which(cumu >90 &pct< 5)[1],sort(which((pct[1:length(pct) - 1] -pct[2:length(pct)]) >0.1),
                                              decreasing = T)[1]+1)

ElbowPlot(cur_seu)$data %>% ggplot() +
  geom_point(aes(x = dims, y = stdev)) +
  geom_vline(xintercept =  pc.use, color = "darkred") +
  theme_bw() + labs(title = "ELboe plot :quantitative approach")

#22,22,18
##Harmony去样本间批次
data_seurat_malignant <-RunHarmony(data_seurat_malignant, "Patient")

#Jack#KNN聚类
data_seurat_malignant <- FindNeighbors(data_seurat_malignant, dims = 1:22, reduction = "harmony")
data_seurat_malignant <- FindClusters(data_seurat_malignant, resolution = 0.5)

##非线性降维(UMAP/tSNE)
data_seurat_malignant <- RunUMAP(data_seurat_malignant, reduction = "harmony", dims = 1:22)
print(DimPlot(data_seurat_malignant, group.by = c("seurat_clusters", "Patient"),reduction = "umap", raster = F))

##在metadata填加tissue标签
tissue <- data_seurat_malignant@meta.data$Patient
tissue[which(tissue %in% c("CA1", "CA2", "CA3", "CA4", "CA5", "CA6", "CA7", "CA8"))] <- "PT"
tissue[which(tissue %in% c("LN1", "LN2", "LN3", "LN4", "LN5", "LN6", "LN7", "LN8"))] <- "LNM"
data_seurat_malignant <- AddMetaData(object = data_seurat_malignant, metadata = tissue, col.name = "tissue")
print(DimPlot(data_seurat_malignant, group.by = c("seurat_clusters", "Patient", "tissue"),reduction = "umap", raster = F))

save(data_seurat_malignant,file = "D:/data/GSE167036-8/data_seurat_malignant.RData")


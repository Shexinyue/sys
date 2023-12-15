
##读取h5数据
BiocManager::install("hdf5r")#安装包
library(hdf5r)
library(Seurat) 
data_sample <- Read10X_h5("C:/Users/pp/Desktop/新建文件夹/GSM6177599_NYU_BRCA0_Vis_processed_filtered_feature_bc_matrix.h5") 
data_seurat <- CreateSeuratObject(data_sample,project = "data_sample", min.cells = 3, min.features = 200) 
#matrix.mtx
data_10x <- Matrix::readMM(file = "C:/Users/pp/Downloads/GSM6177600_NYU_BRCA1_lib1_lib1_gene_expression.mtx")

library(dplyr)
library(Seurat)
library(patchwork)

#导入示例数据
pbmc.data <- Read10X("D:/Download/filtered_gene_bc_matrices/hg19")

#创建Seurat对象
#数据集中测到的少于200个基因的细胞（min.features = 200）和少于3个细胞覆盖的基因（min.cells = 3）被过滤掉
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

#质控
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#使用小提琴图可视化QC指标#nFeature_RNA代表每个细胞测到的基因数目，nCount代表每个细胞测到所有基因的表达量之和，percent.mt代表测到的线粒体基因的比例。
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#去除线粒体基因表达比例过高的细胞，和一些极值细胞。
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#上预处理

##标准化 全局缩放归一化方法“lognormalize”,结果log转换
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#鉴定高变基因
#这一步的目的是鉴定出细胞与细胞之间表达量相差很大的基因，用于后续鉴定细胞类型，
#我们使用默认参数，即“vst”方法选取2000个高变基因。
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# 查看最高变的10个基因
top10 <- head(VariableFeatures(pbmc), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#数据缩放mean 0 and variance 1
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#线性降维
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)#print gene
#查看PCA结果
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
#维度dim选择，选择偏高的，后面一致
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
#量化拐点识别筛选贡献小于5%方差和累计贡献90%方差的PC截断点作为拐点
#连续PC之间的差异方差贡献百分比变化小于0.1%的点作为曲线拐点
#ggplot2
cur_seu <- pbmc %>% SCTransform() %>% RunPCA()
pct <- cur_seu[["pca"]]@stdev / sum(cur_seu[["pca"]]@stdev)*100
cumu <- cumsum(pct)
pc.use <- min(which(cumu >90 &pct< 5)[1],sort(which((pct[1:length(pct) - 1] -pct[2:length(pct)]) >0.1),
                                              decreasing = T)[1]+1)
ElbowPlot(cur_seu)$data %>% ggplot() +
  geom_point(aes(x = dims, y = stdev)) +
  geom_vline(xintercept =  pc.use, color = "darkred") +
  theme_bw() + labs(title = "ELboe plot :quantitative approach")

##Harmony去样本间批次
library(harmony)
pbmc<-RunHarmony(pbmc, "Patient")

#Jack#KNN聚类
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:10)
#选择resolution
pbmc <- FindClusters(pbmc, resolution = 0.5,random.seed = 1)
#查看前5个细胞的分类ID
head(Idents(pbmc), 5)

##非线性降维(UMAP/tSNE)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:10)
# 显示在聚类标签
DimPlot(pbmc, reduction = "umap", label = TRUE, raster=FALSE)
# 使用TSNE聚类
pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:10)
# 显示在聚类标签
DimPlot(pbmc, reduction = "tsne", label = TRUE, raster=FALSE)

##FindMarkers 命令找差异表达基因(聚类标志cluster biomarkers)
# cluster 1的标记基因
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
#找出区分cluster 5与cluster 0和cluster 3的所有标记
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# 找出每个cluster的标记与所有剩余的细胞相比较，只报告阳性细胞
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))
#每个聚类前10个差异基因表达热图(如果小于10，则绘制所有标记)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#通过对比我们鉴定的marker gene与已发表的细胞类型特意的基因表达marker，
#可以定义我们划分出来的细胞类群。最后，给我们定义好的细胞类群加上名称
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

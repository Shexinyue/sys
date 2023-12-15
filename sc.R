
install.packages("pacman") # pacman包如果已经安装则不需要这句
pacman::p_load(hdf5r, dplyr, Seurat, patchwork, SingleR, celldex, pheatmap, ggplot2, cowplot, harmony, stringr)
##读取
data_sample <- Read10X_h5("C:/Users/pp/Desktop/新建文件夹/GSM6177599_NYU_BRCA0_Vis_processed_filtered_feature_bc_matrix.h5") #h5,正常Read10X
data_seurat <- CreateSeuratObject(data_sample,project = "data_sample", min.cells = 3, min.features = 200) #names.field,names.delim,meta.data

load("raw.seurat.obj.RData")
data_seurat <- raw.seurat.obj

setwd("D:/data/OMIX001111")
data_seurat <- readRDS("OMIX001111-20-01.int.total.rawRNA.RDS")

##预处理
#读取蛋白编码基因文件
library(readr)
protein_coding <- read_delim("D:/sc-RNAseq_code/protein_coding.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
#保留蛋白编码基因
data_seurat <- data_seurat[which(rownames(data_seurat) %in% protein_coding$gene_name),]
#质量控制
data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")

#线粒体基因所占比例, MT-开头对应线粒体, 占比越多, 表明死细胞越多
data_seurat <- subset(data_seurat, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 20)

##标准化 
#全局缩放归一化方法“lognormalize”,结果log转换
data_seurat <- NormalizeData(data_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#鉴定高变基因
data_seurat <- FindVariableFeatures(data_seurat, selection.method = "vst", nfeatures = 2000)

#数据缩放mean 0 and variance 1
data_seurat <- ScaleData(data_seurat)
#线性降维PCA
data_seurat <- RunPCA(data_seurat, features = VariableFeatures(object = data_seurat))
#诊断图
#量化拐点识别筛选贡献小于5%方差和累计贡献90%方差的PC截断点作为拐点/连续PC之间的差异方差贡献百分比变化小于0.1%的点作为曲线拐点
cur_seu <- data_seurat
pct <- cur_seu[["pca"]]@stdev / sum(cur_seu[["pca"]]@stdev)*100
cumu <- cumsum(pct)
pc.use <- min(which(cumu >90 &pct< 5)[1],sort(which((pct[1:length(pct) - 1] -pct[2:length(pct)]) >0.1),
                                              decreasing = T)[1]+1)
ElbowPlot(cur_seu)$data %>% ggplot() +
  geom_point(aes(x = dims, y = Stdev)) +
  geom_vline(xintercept =  pc.use, color = "darkred") +
  theme_bw() + labs(title = "ELboe plot :quantitative approach")

#dims=1:16
##Harmony去样本间批次
data_seurat <-RunHarmony(data_seurat, "Patient")

#Jack#KNN聚类
data_seurat <- FindNeighbors(data_seurat, dims = 1:16, reduction = "pca")
# 根据轮廓系数选取分辨率resolution, approxSilhouette函数
resolutions=seq(0.2,1.4,by = 0.2)
library(bluster)
data_seurat <- FindClusters(data_seurat, random.seed = 0, resolution = resolutions, verbose = FALSE)
silhouette <- c()
for (i in 1:length(resolutions)) {
  clusters <- data_seurat@meta.data[[paste0("RNA_snn_res.",resolutions[i])]]
  
  sil.approx <- approxSilhouette(as.matrix(data_seurat@reductions[["pca"]]@cell.embeddings), ##每个细胞在PC轴上的坐标
                                 clusters=clusters)
  sil.data <- as.data.frame(sil.approx)
  sil.avg <- mean(sil.data$width)
  silhouette <- c(silhouette, sil.avg)
  i=i+1
}
resolution_max <- resolutions[which.max(silhouette)]
#0.2
data_seurat <- FindClusters(data_seurat, resolution = 0.2, random.seed = 0)
#查看前5个细胞的分类ID
head(Idents(data_seurat), 5)
##非线性降维(UMAP/tSNE)
data_seurat <- RunUMAP(data_seurat, reduction = "pca", dims = 1:16)
DimPlot(data_seurat, reduction = "umap", raster = F)
# 显示在聚类标签
plot1 <- DimPlot(data_seurat, reduction = "umap", raster = F, group.by = "cellchat.id")

labels_seurat <- as.character(data_seurat@meta.data[["cellchat.id"]])
#Lymphoid
labels_seurat[which(labels_seurat %in% c("GDT"))] = "T_cells"
labels_seurat[which(labels_seurat %in% c("B.cells.NR4A1.2.3", "B.cells.Plasma.cells", "B.cells.Unknown", "naive.B.cells", "Plasma.cells", "class-switched.B.cells"))] = "B_cells"
labels_seurat[which(labels_seurat %in% c("NK", "NK.PBMC", "NKT.CD8.EM", "NKT.RA.PBMC", "NKT.PBMC"))] = "NK_cells"
#Myeloid
#labels_seurat[which(strsplit(labels_seurat[i], split = ".",fixed = T)[[1]][1] %in% "Mac")] = "Macrophage"
labels_seurat[which(labels_seurat %in% c("Mac.LST1", "Mac.marker-low", "Mono.CD16+", "Mono.Mac.S100A8.9", "Monocytes.CD14+", "Monocytes.CD14+.PPBP", "Monocytes.CD14+.S100A8.9"))] = "Monocyte"
#Stromal_cells
labels_seurat[which(labels_seurat %in% c("Fibroblasts.HSP.AP1", "Wound_healing_Fibroblasts"))] = "Fibroblasts"
for (i in 1:128678) {
  if(strsplit(labels_seurat[i], split = ".",fixed = T)[[1]][1] %in% c("CD4", "CD8")){
    labels_seurat[i] = "T_cells"
  }
  if(strsplit(labels_seurat[i], split = ".",fixed = T)[[1]][1] %in% "Mac"){
    labels_seurat[i] = "Macrophage"
  }
  }

data_seurat$cellchat.id.2 <- labels_seurat
plot2 <- DimPlot(data_seurat, group.by = "cellchat.id.2",reduction = "umap", raster = F)


#进行singleR注释
load("D:/data/singler_refer/hpca.se.RData")
hpca.se

data_seurat_for_SingleR <- GetAssayData(data_seurat, slot="data") ##获取标准化矩阵
data_seurat.hesc <- SingleR(test = data_seurat_for_SingleR, ref = hpca.se, labels = hpca.se$label.main)
data_seurat.hesc

table(data_seurat.hesc$labels, data_seurat@meta.data$seurat_clusters)

#绘制umap/tsne图
data_seurat@meta.data$labels <-data_seurat.hesc$labels

print(DimPlot(data_seurat, group.by = c("seurat_clusters", "labels"),reduction = "umap", raster = F))


##直接修改类群的值
labels_seurat <- as.character(data_seurat@meta.data[["labels"]])

#Lymphoid
labels_seurat[which(labels_seurat %in% c("B_cell", "Pre-B_cell_CD34-", "Pro-B_cell_CD34+"))] = "B_cell"
#Myeloid
labels_seurat[which(labels_seurat %in% c( "MEP", "Neutrophils", "GMP", "DC", "CMP", "Pro-Myelocyte"))] = "Myeloid"
#Stromal_cells
labels_seurat[which(labels_seurat %in% c( "Osteoblasts", "Chondrocytes", "Smooth_muscle_cells", "Embryonic_stem_cells", 
                                          "Keratinocytes", "Tissue_stem_cells", "Neurons", "MSC", 
                                          "Hepatocytes", "HSC_-G-CSF", "HSC_CD34+", "iPS_cells"))] = "Stromal_cells"
data_seurat$labels2 <- labels_seurat
plot3 <- DimPlot(data_seurat, group.by = c("seurat_clusters", "labels2"),reduction = "umap", raster = F)

#其他类型的基质marker表达
genes_to_check <- c("COL1A1", "DCN", "COL1A2", "C1R", "FAP", "THY1", "COL6A1", "LUM", 
                     "ACTA2", "ACTG2", "ACTN2", "MYL2", "MYH2", "TAGLN", "CNN1")
genes_to_check = unique(intersect(rownames(data_seurat),genes_to_check))
genes_to_check

P3 <- DotPlot(data_seurat, features = genes_to_check, assay='RNA'  )  + coord_flip()
P4 <- VlnPlot(object = data_seurat, features =genes_to_check,log =T )
P5 <- FeaturePlot(object = data_seurat, features=genes_to_check )
#簇9，10，20划分为基质,重新划簇

labels_seurat[which(labels_seurat %in% c("Neurons", "MSC", "Keratinocytes", "iPS_cells", "HSC_CD34+",
                                        "Hepatocytes", "Gametocytes", "Embryonic_stem_cells",
                                        "BM", "BM & Prog.", "Astrocyte", "Tissue_stem_cells"))] = "Stromal_cells"


#AddModuleScore打分函数，依据基因的平均表达水平进行分析，最后得到一组基因的score
data_seurat#输入1：Seurat对象
gene_list#输入2：基因list
data_seurat <- AddModuleScore(
  object = data_seurat,
  features = gene,
  ctrl = 100,
  name = 'CD_Features'
)
colnames(data_seurat@meta.data)
colnames(data_seurat@meta.data)[8] <- 'Inflammatory_Score'

VlnPlot(data_seurat, features = 'Inflammatory_Score',pt.size = 0, adjust = 2,group.by = "orig.ident")

mydata<- FetchData(data_seurat,vars = c("UMAP_1","UMAP_2","Inflammatory_Score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = Inflammatory_Score))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('blue',"cyan",'green','yellow','orange','red'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

save(data_seurat,file = "umap_seurat.RData")

###infercnv+识别恶性细胞

load("D:/data/GSE167036-8/umap_seurat.RData")
load("D:/data/GSE167036-8/data_seurat_malignant.RData")
#######堆叠柱状图
#总体细胞簇类型
color.celltype <- c("#FC8D62","#E78AC3",  "#8DA0CB","#66C2A5",  "#A6D854", "#FFD92F")
table(data_seurat$labels)#y
table(data_seurat$Patient)#x
#prop.table(data_seurat$Patient)
table(data_seurat$labels, data_seurat$Patient)
Cellratio <- prop.table(table(data_seurat$labels, data_seurat$Patient), margin = 2)%>%
  as.data.frame()

ggplot(Cellratio) +
  geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = "identity", position = "stack") + # 如果把 "stack" 改成 "dodge"，可以变成分组柱状图
  xlab(NULL) + ylab('Percentage')+
  scale_fill_manual(values = color.celltype)+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x  = element_text(angle=90,# 设置旋转的角
                                    size=10), # 字体的大小
        axis.text.y  = element_text(size=10), 
        legend.title=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),#去掉网格线
        panel.background = element_blank(),#去掉背景
        axis.line = element_line(size = 1,colour = "black")#加坐标轴
  )
ggsave("9celltype_patient.pdf",width=6,height=4)

#恶性细胞簇患者比例
color.patient <- brewer.pal(10, "Set3")

table(data_seurat_malignant$Patient)
table(data_seurat_malignant$seurat_clusters)
prop.table(table(data_seurat_malignant$seurat_clusters))
table(data_seurat_malignant$Patient, data_seurat_malignant$seurat_clusters)
PatCellratio <- prop.table(table(data_seurat_malignant$Patient, data_seurat_malignant$seurat_clusters),margin = 2)
PatCellratio <- as.data.frame(PatCellratio)

PatCellratio$Var1 <- as.character(PatCellratio$Var1)
PatCellratio$Var1[which(PatCellratio$Var1 %in% c("CA1","LN1"))] <- "patient1"
PatCellratio$Var1[which(PatCellratio$Var1 %in% c("CA1","LN1"))] <- "patient1"
PatCellratio$Var1[which(PatCellratio$Var1 %in% c("CA2","LN2"))] <- "patient2"
PatCellratio$Var1[which(PatCellratio$Var1 %in% c("CA3","LN3"))] <- "patient3"
PatCellratio$Var1[which(PatCellratio$Var1 %in% c("CA4","LN4"))] <- "patient4"
PatCellratio$Var1[which(PatCellratio$Var1 %in% c("CA5","LN5"))] <- "patient5"
PatCellratio$Var1[which(PatCellratio$Var1 %in% c("CA6","LN6"))] <- "patient6"
PatCellratio$Var1[which(PatCellratio$Var1 %in% c("CA7","LN7"))] <- "patient7"
PatCellratio$Var1[which(PatCellratio$Var1 %in% c("CA8","LN8"))] <- "patient8"

ggplot(PatCellratio) +
  geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = "identity", position = "stack") + # 如果把 "stack" 改成 "dodge"，可以变成分组柱状图
  xlab(NULL) + ylab('Percentage')+
  scale_fill_manual(values = color.patient)+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x  = element_text(angle=90,# 设置旋转的角
                                    size=10), # 字体的大小
        axis.text.y  = element_text(size=10), 
        legend.title=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),#去掉网格线
        panel.background = element_blank(),#去掉背景
        axis.line = element_line(size = 1,colour = "black")#加坐标轴
  )
ggsave("10Patient_maligcellcluster.pdf",width=6,height=4)

#恶性细胞簇原发转移
color.tumor <- c("#3C5488B2","#00A087B2")

table(data_seurat_malignant$Patient)
table(data_seurat_malignant$seurat_clusters)
prop.table(table(data_seurat_malignant$seurat_clusters))
table(data_seurat_malignant$Patient, data_seurat_malignant$seurat_clusters)
malCellratio <- prop.table(table(data_seurat_malignant$Patient, data_seurat_malignant$seurat_clusters),margin = 2)
malCellratio <- as.data.frame(malCellratio)

malCellratio$Var1 <- as.character(malCellratio$Var1)
malCellratio$Var1[which(malCellratio$Var1 %in% c("CA1","CA2","CA3","CA4","CA5","CA6","CA7","CA8"))] <- "PT"
malCellratio$Var1[which(malCellratio$Var1 %in% c("LN1","LN2","LN3","LN4","LN5","LN6","LN7","LN8"))] <- "LNM"
ggplot(malCellratio) +
  geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = "identity", position = "stack") + # 如果把 "stack" 改成 "dodge"，可以变成分组柱状图
  xlab(NULL) + ylab('Percentage')+
  scale_fill_manual(values = color.patient)+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x  = element_text(angle=90,# 设置旋转的角
                                    size=10), # 字体的大小
        axis.text.y  = element_text(size=10), 
        legend.title=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),#去掉网格线
        panel.background = element_blank(),#去掉背景
        axis.line = element_line(size = 1,colour = "black")#加坐标轴
  )
  ggsave("11LNMPT_maligcellcluster.pdf",width=6,height=4)

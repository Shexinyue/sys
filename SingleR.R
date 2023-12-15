
library(SingleR)

library(celldex)
library(Seurat)
library(pheatmap)
##下载注释数据库
hpca.se <- celldex::HumanPrimaryCellAtlasData()
load("D:/data/singler_refer/hpca.se.RData")
hpca.se


#查看seuret聚类结果
load("pbmc_tutorial.RData")
pbmc

meta=pbmc@meta.data #pbmc的meta文件，包含了seurat的聚类结果
head(meta)

#查看umap和tsne图
plot1 <- DimPlot(pbmc, reduction = "umap", label = TRUE)
plot2<-DimPlot(pbmc, reduction = "tsne",
               label = TRUE)
plot1 + plot2

#进行singleR注释
pbmc_for_SingleR <- GetAssayData(pbmc, slot="data") ##获取标准化矩阵
pbmc.hesc <- SingleR(test = pbmc_for_SingleR, ref = hpca.se, labels = hpca.se$label.main) #
pbmc.hesc

refdata <- get(load("ref_Monaco_114s.RData"))
sce_for_SingleR <- GetAssayData(sce, slot="data")
clusters <- sce@meta.data$seurat_clusters
pred.hesc <- SingleR(test = sce_for_SingleR, ref = refdata, 
                     labels = refdata$label.fine,
                     #因为样本主要为免疫细胞（而不是全部细胞），因此设置为label.fine
                     method = "cluster", clusters = clusters,
                     #这里我们为上一步分的9个cluster注释celltype
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(pred.hesc$labels)


#seurat 和 singleR的table表
#查看新注释的标签与seurat分类的结果的交叠情况
table(pbmc.hesc$labels,meta$seurat_clusters)

#绘制umap/tsne图
pbmc@meta.data$labels <-pbmc.hesc$labels

print(DimPlot(pbmc, group.by = c("seurat_clusters", "labels"),reduction = "umap"))


##########人工注释单细胞群
genes_to_check <- list(Cancer_cells = c("EPCAM", "KRT18", "KRT7", "KRT8", "CD24"),
                       Epithelial_cell = c("CAPS", "EPCAM"),
                       Cancer_stem_cell = c("OLFM4", "SOX2", "LGR5", "CCKBR"),
                       proliferating = c("MKI67"),
                       T_cells = c("CD3D", "CD3E", "CD3G", "CD2", "TRBC2", "TRBC1", "TRAC", "IL7R", "CD3", "CD4", "CD8"),
                       B_cells = c("CD79A", "CD79B", "CD19", "MS4A1", "IGHG3", "IGKC"),
                       NK_cells = c("FGFBP", "IL7R", "NKG7", "NCAM1", "NCR1", "NKG2"),
                       plasma_cell = c("JCHAIN", "MZB1"),
                       Macrophage = c("CD163", "CD68", "CSF1R", "CD14"),
                       Monocyte = c("HLA-DR", "FCGR3A", "CD68", "ANPEP", "ITGAX", "CD14", "ITGAM", "CD33"),
                       Dendritic = c("CCL17", "CD1C", "CD207", "LILRA4"),
                       Granulocyte = c("S100A12", "S100A8", "S100A9", "VCAN", "CD11B", "GR-1", "LY6G", "CD117", "CD43", "CD49B", "CEBPE", "CD123"),
                       Neutrophil = c("CSF3R"),
                       Mast_cell = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"),
                       Fibroblasts = c("COL1A1", "DCN", "COL1A2", "C1R", "FAP", "THY1", "COL6A1", "LUM"),
                       Smooth_muscle_cell = c("ACTA2", "ACTG2", "ACTN2", "MYL2", "MYH2", "TAGLN", "CNN1"),
                       Pericte = c("ACTA2", "RGS5"),
                       Endothelial = c("PECAM1", "VWF", "CLDN5", "CDH5", "FLT1", "RAMP2")
)
#Fibroblasts , Osteoblasts, Chondrocytes, Smooth_muscle_cells
genes_to_check = unique(intersect(rownames(data_seurat),genes_to_check))#mark与矩阵基因的交集
genes_to_check

P3 <- DotPlot(data_seurat, features = genes_to_check, assay='RNA'  )  + coord_flip()
P4 <- VlnPlot(object = data_seurat, features =genes_to_check,log =T )
P5 <- FeaturePlot(object = data_seurat, features=genes_to_check )

##直接修改类群的值
data_seurat_2 <- data_seurat
new.cluster.ids <- c(rep("Bladder Cancer",8),"Endothelial cell",rep("Bladder Cancer",4))

names(new.cluster.ids) <- levels(data_seurat_2)
sce3 <- RenameIdents(data_seurat_2, new.cluster.ids)#修改Idents

DimPlot(data_seurat_2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(pbmc,                                #seurat对象(数据)
        reduction = 'umap',                  #降维类型，umap，tsne，pca
        group.by = 'celltype',               #分组名称
        pt.size = 1,                         #图中点的大小
        #split.by = 'seurat_clusters',       #拆分绘图的名称
        label = T)                           #图中是否添加label信息

#修改簇的cell.type
DimPlot(data_seurat, group.by = c("seurat_clusters", "labels"),reduction = "umap",label = T)

cell.types = as.character(data_seurat@meta.data[["seurat_clusters"]])
cell.types[which(cell.types %in% c("2", "5","7","8","13","14","16","19"))] = "Epithelial_cells"  
cell.types[which(cell.types %in% c("0", "1", "3","9","10","11","15","23"))] = "Lymphoid"   
cell.types[which(cell.types %in% c("6","17","22","24"))] = "Myeloid"  
cell.types[which(cell.types %in% c("4","18","20","21"))] = "Stomal"
cell.types[which(cell.types %in% "12")] = "Endothelial_cells"
table(cell.types)

data_seurat$celltype <- cell.types

##细胞类型
#Epithelial_cells: Epithelial_cells

#Lymphoid: B_cell,Pre-B_cell_CD34-,Pro-B_cell_CD34+, NK_cell, T_cells,

#Myeloid: Macrophage, Monocyte,Neutrophils,Pro-Myelocyte, DC, CMP, GMP,MEP

#Stromal_cells: Fibroblasts , Osteoblasts, Chondrocytes, Smooth_muscle_cells

#others: Keratinocytes(角质细胞), Astrocyte(星形胶质细胞)，BM(骨髓细胞)，BM & Prog,Gametocytes
#(配子母细胞)，Hepatocytes(肝细胞)，Neuroepithelial_celI(神经上皮细胞),Neurons(神经细胞)，
#Embryonic_stem_cells(胚胎干细胞)，MSC(间充质干细胞), iPS_cells(诱导性多能干细胞)，
#Tissue_stem_cells(组织干细胞)，HSC_-G-CSF(造血干细胞)，HSC_CD34+(造血干细胞)
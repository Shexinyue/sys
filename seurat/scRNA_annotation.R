###scRNA-seq数据使用SingleR注释细胞类型
#note： 需要提前安装SingleR包
##输入参数：data_seurat - 聚类后的seurat对象
#           dir.reference - 参考数据库路径
##输出结果：注释后的seurat对象，注释信息存在data_seurat@meta.data$cell_type中
#
scRNA_annotated <-function(data_seurat, dir.reference = "Z:/0.Public/Data/reference/hpca.se.RData"){
  library(SingleR)
  load(dir.reference)
  data_seurat_for_SingleR <- GetAssayData(data_seurat, slot="data") 
  data_seurat.hesc <- SingleR(test = data_seurat_for_SingleR, ref = hpca.se, labels = hpca.se$label.main)
  data_seurat@meta.data$cell_type <-data_seurat.hesc$labels
  return(data_seurat)
}

#测试代码
#load("Z:/Users/ShenXinYue/scRNA/pbmc_annotation.RData")
#pbmc <- scRNA_annotated(pbmc)
#table(pbmc$cell_type, pbmc$seurat_clusters)##查看注释的细胞类型和簇的交叠情况
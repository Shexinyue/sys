###根据inferVNV分析的结果选出恶性细胞
#
# Note: 如果有多个患者，改变infercnv_out.dir运行多次
#
##输入参数：data_seurat - seurat对象要求该seurat对象运行过scRNA_infercnv函数
#           infercnv_out.dir - infercnv函数的输出文件的路径
#           clusters - 根据heatmap1.pdf和CNV level.pdf选出的CNV 得分高的簇
##输出：seurat对象，在meta.data中添加malignant（"malignant"，"Benign"），细胞的良恶性标签。

scRNA_malignantcell <-function(data_seurat, infercnv_out.dir, clusters){
  kmeans_df_s = read.table(paste0(infercnv_out.dir,"/kmeans_df_s.txt"),sep="\t")
  malignant <- rownames(kmeans_df_s[which(kmeans_df_s$kmeans_class %in% clusters & 
                                            kmeans_df_s$class %in% "cancer" ),])
  if((is(pbmc@meta.data$malignant) == "NULL")[1]){
    anno.df=data.frame(
      CB <- colnames(data_seurat),
      class <- c(rep("Benign",length(colnames(data_seurat))))
    )
    anno.df$class[which(anno.df$CB %in% colnames(malignant))] <- "malignant"
    data_seurat <- AddMetaData(data_seurat, metadata = anno.df$class, col.name = "malignant")
  }
  else{
    data_seurat@meta.data$malignant[which(colnames(data_seurat) %in% malignant)] = "malignant"
  }
  return(data_seurat)
}

#测试代码
#load("Z:/0.Public/Code/Analysis/scRNA_analysis/5.scRNA_inferCNV/TestData/pbmc_infercnv.RData")
#pbmc <- scRNA_malignantcell(pbmc, 
#                            infercnv_out.dir = "D:/data/shili",
#                            clusters = c(1,2,3,6,7))
#table(data_seurat@meta.data$malignant)
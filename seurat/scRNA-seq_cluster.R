###scRNA-seq聚类分析
##输入参数：data_seurat - 预处理后的seurat对象
#           resolutions – FindClusters的分辨率,默认为从0.4到1.4间隔为0.2
#           nfeatures - 高变基因数，默认为2000
cluster <-function(data_seurat, resolutions=seq(0.4,1.4,by = 0.2), nfeatures = 2000){
  #标准化
  data_seurat <- NormalizeData(data_seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  #鉴定高变基因
  data_seurat <- FindVariableFeatures(data_seurat, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  #数据缩放mean 0 and variance 1
  all.genes <- rownames(data_seurat)
  data_seurat <- ScaleData(data_seurat, features = all.genes, verbose = FALSE)
  #线性降维
  data_seurat <- RunPCA(data_seurat, features = VariableFeatures(object = data_seurat), verbose = FALSE)
  #设定PC数,方法为量化拐点
  cur_seu <- data_seurat 
  pct <- cur_seu[["pca"]]@stdev / sum(cur_seu[["pca"]]@stdev)*100
  cumu <- cumsum(pct)
  pc.use <- min(which(cumu >90 &pct< 5)[1],sort(which((pct[1:length(pct) - 1] -pct[2:length(pct)]) >0.1),
                                                decreasing = T)[1]+1)
  dims = 1:pc.use
  data_seurat <- FindNeighbors(data_seurat, reduction = "pca", dims = dims, verbose = FALSE)
  #根据轮廓系数选取resolution
  library(bluster)#approxSilhouette
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
  
  pdf("resolution_SC.pdf")
  barplot(silhouette,names.arg=resolutions,xlab="Resolution",ylab="Silhouette Coefficient")
  dev.off()
  
  data_seurat <- FindClusters(data_seurat, random.seed = 0, resolution = resolution_max, verbose = FALSE)
  return(data_seurat)
}

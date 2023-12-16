#################数据预处理#############
## 运行前先设置路径！！！
## 先运行一遍，根据质控图选择最佳参数！！！

## 相关参数:
## raw.data: 单细胞原始矩阵, 行为基因, 列为细胞;
## protein_coding_path: 蛋白编码基因注释文件的路径, 文件基因symbol列名为gene_name; 
## project: 癌症类型;
## min.features: 保留至少检测到多少个基因的细胞, 默认min.features = 200;
## min.cells: 保留至少在多少个细胞中检测到的基因, 默认min.cells = 5;
## max.features: 过滤检测到超过多少个基因的细胞;
## min.counts: 保留UMI数目大于多少的细胞;
## max.counts: 保留UMI数目小于多少的细胞;
## percent.mt: 线粒体基因比例, 默认percent.mt = 20;

## result: 处理好的高质量的seurat object;

Data_preprocessing <- function(raw.data,protein_coding_path,project,min.features = 200,min.cells = 5,max.features,min.counts,max.counts,percent.mt = 20){
  library(Seurat)
  ## 读取蛋白编码基因文件
  protein_coding <- read.table(protein_coding_path, header = T, sep = "\t")
  ## 保留蛋白编码基因
  raw.data <- raw.data[which(rownames(raw.data) %in% protein_coding$gene_name),]

  ## 创建seurat对象
  raw.seurat.obj <- CreateSeuratObject(counts = raw.data, project = project, min.features = min.features, min.cells = min.cells, names.field = 1, names.delim = "-")

  ## 质量控制
  raw.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(raw.seurat.obj, pattern = "^MT-")
  raw.seurat.obj[["percent.HB"]]<- PercentageFeatureSet(raw.seurat.obj,features="HBB")
  pdf("QC.pdf", width = 8, height = 4)
  VlnPlot(raw.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), pt.size = 0.1, ncol = 4)
  dev.off()
  ## 细胞nFeature_RNA和nCount_RNA太高，可能是双细胞
  ## 线粒体基因所占比例, MT-开头对应线粒体, 占比越多, 表明死细胞越多
  raw.seurat.obj <- subset(raw.seurat.obj, subset = nFeature_RNA > min.features & nFeature_RNA < max.features & nCount_RNA > min.counts & nCount_RNA < max.counts & percent.mt < percent.mt)
  
  return(raw.seurat.obj)

}

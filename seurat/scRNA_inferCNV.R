###scRNA-seq数据inferVNV分析
#
# Note: 需要提前安装infercnv包
#       建议提取出每个患者分别运行
#
##输入参数：data_seurat - seurat对象要求该seurat对象有count数据，且“seurat.obj@meta.data$cell_type ”中存储细胞注释的结果，最少要有“恶性”和”非恶性“两种标签。
#           gene.order - geneOrderingFile文件的路径
#           ref.group.name - character向量，确定哪些标签的细胞被指定作为参考的正常细胞集
#           cutoff - 数字，去掉细胞中平均read数<cutoff的基因（default: 0.1）。 
#                   注意：对于Smart-seq2来说cutoff=1 works well ，对于10x Genomics来说 cutoff=0.1 works well
#           out.dir - 结果文件输出路径
#           analysis_mode - 可选（samples，subclusters，cells），选择在患者、亚群（克隆）或细胞水平上预测CNV区域(default: samples)
#           cluster_by_groups - 逻辑值，是否按组单独聚集(例如患者)
#           denoise - 逻辑值，是否过滤正常参考细胞当中的信号噪音(default: TRUE)
#           HMM - 逻辑值，是否运行HMM预测CNV区域 (default:FALSE)
#                 注！经测试，运行HMM与否对最终的结果文件没有任何影响。
#                   ！运行HMM会单独产生HMM开头的文件，HMM中明确了每个基因的拷贝数状态，如果不需要详细分析数据的拷贝数情况，则用去噪后的结果即可。
#                   ！如果只是为了区分肿瘤细胞，不用运行HMM
#           K - kmeans的聚类数（default:7）
##输出结果：infercnv结果文件
#           heatmap1.pdf - cnv得分聚类热图
#           CNV level.pdf - 每个簇的CNV score

scRNA_inferCNV <-function(data_seurat, 
                          gene.order = "Z:/0.Public/Data/reference/geneOrderingFile.txt", 
                          ref.group.name,
                          cutoff = 0.1, 
                          out.dir,
                          analysis_mode = "samples", 
                          cluster_by_groups = F, 
                          denoise = T, 
                          HMM = F,
                          K = 7){
  library(infercnv)
  library(Seurat)

  matrix_counts <- as.matrix(data_seurat[["RNA"]]$counts)
  ##创建infercnv对象
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = matrix_counts,
                                       annotations_file= data.frame(data_seurat@meta.data$cell_type,row.names=colnames(data_seurat)), 
                                       gene_order_file= gene.order, 
                                       ref_group_names= "normal")
  
  ##运行inferCNV
  infercnv_obj <- infercnv::run(infercnv_obj, 
                      cutoff= cutoff, 
                      out_dir= out.dir,
                      analysis_mode= analysis_mode,
                      cluster_by_groups = cluster_by_groups, 
                      denoise= denoise,
                      HMM= HMM)
  #
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(tidyverse)
  
  #####
  infercnv_obj = readRDS(paste0(out.dir,"/run.final.infercnv_obj"))
  expr <- infercnv_obj@expr.data
  normal_loc <- infercnv_obj@reference_grouped_cell_indices
  normal_loc <- as.character(unlist(normal_loc))
  cancer_loc <- infercnv_obj@observation_grouped_cell_indices
  cancer_loc <- as.character(unlist(normal_loc))
  
  anno.df=data.frame(
    CB=c(colnames(expr)[normal_loc],colnames(expr)[cancer_loc]),
    class=c(rep("normal",length(normal_loc)),rep("cancer",length(cancer_loc)))
  )
  
  gn <- rownames(expr)
  geneFile <- read.table(gene.order,header = F,sep = "\t",stringsAsFactors = F)
  rownames(geneFile)=geneFile$V1
  sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
  expr=expr[intersect(gn,geneFile$V1),]
  sub_geneFile$V2 <- paste("chr",sub_geneFile$V2,sep = "")
  
  set.seed(20231129)
  kmeans.result <- kmeans(t(expr), K)
  kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
  kmeans_df$CB=rownames(kmeans_df)
  kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #合并
  kmeans_df_s=arrange(kmeans_df,kmeans_class) #排序
  rownames(kmeans_df_s)=kmeans_df_s$CB
  kmeans_df_s$CB=NULL
  kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:7排序
  
  #定义热图的注释，及配色
  top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
  color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:K] #类别数
  names(color_v)=as.character(1:K)
  left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("cancer"="red","normal" = "blue"),kmeans_class=color_v))
  
  #下面是绘图
  pdf(paste0(out.dir,"/heatmap1.pdf"),width = 15,height = 10)
  ht = Heatmap(t(expr)[rownames(kmeans_df_s),], #绘图数据的CB顺序和注释CB顺序保持一致
               col = colorRamp2(c(0.85,1,1.15), c("#377EB8","#F0F0F0","#E41A1C")), #如果是10x的数据，这里的刻度会有所变化
               cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
               column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), #这一步可以控制染色体顺序，即使你的基因排序文件顺序是错的
               column_gap = unit(2, "mm"),
               
               heatmap_legend_param = list(title = "Modified expression",direction = "vertical",
                                           title_position = "leftcenter-rot",at=c(0.85,1,1.15),
                                           legend_height = unit(3, "cm")),
               
               top_annotation = top_anno,left_annotation = left_anno, #添加注释
               row_title = NULL,column_title = NULL,
               use_raster = T)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  #每一类对应的CB保存在kmeans_df_s数据框中
  write.table(kmeans_df_s, file = paste0(out.dir,"/kmeans_df_s.txt"), quote = FALSE, sep = '\t', row.names = T, col.names = T)
  
  ###infercnv_score
  CNV_score=as.data.frame(colSums((expr-1)*(expr-1)))
  colnames(CNV_score)="CNV_score"
  CNV_score$CB=rownames(CNV_score)
  kmeans_df_s$CB=rownames(kmeans_df_s)
  CNV_score=CNV_score%>%inner_join(kmeans_df_s,by="CB")
  
  CNV_score%>%ggplot(aes(kmeans_class,CNV_score))+geom_violin(aes(fill=kmeans_class),color="NA")+
    scale_fill_manual(values = color_v)+
    theme_bw()
  ggsave(paste0(out.dir,"/CNV level.pdf"),width = 10,height = 6,units = "cm")
}

#测试代码
#load("Z:/0.Public/Code/Analysis/scRNA_analysis/5.scRNA_inferCNV/TestData/pbmc_infercnv.RData") 
#scRNA_inferCNV(data_seurat = pbmc,
#               ref.group.name = c("B_cell", "CMP", "DC", "Monocyte", "NK_cell", 
#                                  "Platelets", "Pre-B_cell_CD34-", "Pro-B_cell_CD34+"), 
#               out.dir = "D:/data/shili")
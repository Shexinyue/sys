#安装包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")
library(infercnv)
library(Seurat)
library(ggplot2)

setwd("D:/data/inferCNV")
#查看数据
load("D:/data/GSE167036-8/umap_seurat.RData")
seurat_object <- data_seurat
levels(seurat_object)

#筛选分析需要的细胞类型
seurat_object <- subset(seurat_object, labels == "Epithelial_cells"|
                        labels == "Myeloid"|
                          labels == "Lymphoid")
#提取患者
seurat_object_CA1 <- seurat_object[,which(seurat_object@meta.data$Patient %in% "CA1")]

#提取counts                                                                    
counts <- as.matrix(seurat_object[["RNA"]]@counts)
saveRDS(as.matrix(seurat_object_CA1[["RNA"]]@counts), "sc.10x.counts.CA1.matrix")
#细胞类型注释文件
anno= seurat_object$labels
write.table(seurat_object_CA1$labels, "D:/data/inferCNV/seurat.obj.CA1.label.txt", sep = "\t", quote = F, col.names = F)
#基因注释文件
gene_order <- "D:/data/inferCNV/geneOrderingFile.txt"

matrix_counts <- readRDS("sc.10x.counts.CA1.matrix")
ref <- read.delim("seurat.obj.label.txt", header = F)
#创建inferCNV对象
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = matrix_counts,
                                    annotations_file = "seurat.obj.CA1.label.txt",
                                    delim="\t",
                                    gene_order_file = "geneOrderingFile.txt",
                                    ref_group_names=c("Lymphoid", "Myeloid"))
#运行inferCNV
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,
                             out_dir = "D:/data/inferCNV/CA1",
                             analysis_mode="samples",
                             cluster_by_groups = F,       
                             HMM = FALSE,
                             denoise = TRUE,
                             )    
#用inferCNV的官方的画图函数
infercnv::plot_cnv(infercnv_obj, #上两步得到的infercnv对象
                   plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
                   output_filename = "better_plot",output_format = "pdf", #保存为pdf文件
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2))) #改颜色

##结果读取
inferCNV_dir = "D:/data/inferCNV/CA1/"
infercnv_result <- list()
re_files <- list.files(inferCNV_dir)

cnv_signal	<- readRDS(paste0(inferCNV_dir,re_files[grep("run.final.infercnv_obj",re_files)]))

# 提取observation细胞的CNV signal
infercnv_result[[1]] <- cnv_signal@expr.data[,unlist(cnv_signal@observation_grouped_cell_indices)]

# 提取reference细胞的CNV signal
infercnv_result[[2]] <- cnv_signal@expr.data[,unlist(cnv_signal@reference_grouped_cell_indices)]

# 提取基因坐标信息
infercnv_result[[3]] <- cnv_signal@gene_order

names(infercnv_result)<-c("CNV_signal","CNV_signal_RefNormal","Gene_coordinate")

##获得结果文件infercnv_result

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)

# 基因坐标信息
chr <- factor(infercnv_result$Gene_coordinate$chr)
 
# CNVt图谱/obs
if(type=="obs"){expr <- infercnv_result$CNV_signal}
if(type=="ref"){expr <- infercnv_result$CNV_signal_RefNormal}

# 染色体配色
get_group_color_palette <- function () {
  return(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3")))
}  
color <- get_group_color_palette()(length(unique(chr)))

top_color <- HeatmapAnnotation(
  block = anno_block(gp = gpar(fill = color), # 设置填充色
                     labels = levels(chr), 
                     labels_gp = gpar(cex = 0.9, col = "black"))) 

# 画图
key_range<-c(0.7,1,1.3)
ht_tumor <- Heatmap(as.matrix(t(expr)),
                    cluster_rows = T,
                    clustering_method_rows="ward.D",
                    cluster_columns = F,
                    show_column_names = F,
                    show_row_names = F,
                    col=colorRamp2(key_range, c("darkblue", "white", "darkred")),
                    column_split = chr,
                    heatmap_legend_param = list(
                      title = "Modified Expression",
                      title_position = "leftcenter-rot", # 图例标题位置
                      at=c(key_range[1],key_range[3]), #图例范围
                      legend_height = unit(3, "cm") #图例长度
                    ),
                    top_annotation = top_color,
                    row_title = "Observations (Cells)",
                    row_title_side = c("right"),
                    column_title = "Genomic Region",
                    column_title_side = c("bottom"),
                    use_raster = T)
draw(ht_tumor, heatmap_legend_side = "left") # 图例位置


#####
infercnv_obj = readRDS("./CA1/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <- c(normal_loc$Lymphoid,normal_loc$Myeloid)
cancer_loc <- infercnv_obj@observation_grouped_cell_indices
cancer_loc <- cancer_loc$Epithelial_cells

anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[cancer_loc]),
  class=c(rep("normal",length(normal_loc)),rep("cancer",length(cancer_loc)))
)

head(anno.df)
table(anno.df$class)

gn <- rownames(expr)
geneFile <- read.table("geneOrderingFile.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
expr[1:4,1:4]
sub_geneFile$V2 <- paste("chr",sub_geneFile$V2,sep = "")

set.seed(20231003)
kmeans.result <- kmeans(t(expr), 7)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #合并
kmeans_df_s=arrange(kmeans_df,kmeans_class) #排序
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:7排序
head(kmeans_df_s)

#定义热图的注释，及配色
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:7] #类别数
names(color_v)=as.character(1:7)
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("cancer"="red","normal" = "blue"),kmeans_class=color_v))

#下面是绘图
pdf("heatmap1.pdf",width = 15,height = 10)
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
write.table(kmeans_df_s, file = "kmeans_df_s.txt", quote = FALSE, sep = '\t', row.names = T, col.names = T)

###infercnv_score
CNV_score=as.data.frame(colSums((expr-1)*(expr-1)))
colnames(CNV_score)="CNV_score"
CNV_score$CB=rownames(CNV_score)
kmeans_df_s$CB=rownames(kmeans_df_s)
CNV_score=CNV_score%>%inner_join(kmeans_df_s,by="CB")

CNV_score%>%ggplot(aes(kmeans_class,CNV_score))+geom_violin(aes(fill=kmeans_class),color="NA")+
  scale_fill_manual(values = color_v)+
  theme_bw()
ggsave("CNV level.pdf",width = 10,height = 6,units = "cm")
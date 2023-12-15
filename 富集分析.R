library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(pheatmap)

load("D:/data/GSE167036-8/data_seurat_malignant.RData")
table(data_seurat_malignant$seurat_clusters)

data_seurat_malignant.markers <- FindAllMarkers(data_seurat_malignant, only.pos = TRUE,
                                                min.pct = 0.1, logfc.threshold = 1, verbose = FALSE)
all_markers <- subset(data_seurat_malignant.markers, p_val_adj<0.05) #所有差异基因

data_seurat_malignant <- ScaleData(object = data_seurat_malignant,features = rownames(data_seurat_malignant))

pdf(file = "12difgene_heatmap.pdf", width = 13,height = 15)
DoHeatmap(data_seurat_malignant, features = all_markers$gene, disp.min=-2.5, disp.max=2.5)
dev.off()

deg.list <- list(C0 = all_markers$gene[which(all_markers$cluster == "0")],
                 C1 = all_markers$gene[which(all_markers$cluster == "1")],
                 C2 = all_markers$gene[which(all_markers$cluster == "2")],
                 C3 = all_markers$gene[which(all_markers$cluster == "3")],
                 C4 = all_markers$gene[which(all_markers$cluster == "4")],
                 C5 = all_markers$gene[which(all_markers$cluster == "5")],
                 C6 = all_markers$gene[which(all_markers$cluster == "6")],
                 C7 = all_markers$gene[which(all_markers$cluster == "7")],
                 C8 = all_markers$gene[which(all_markers$cluster == "8")],
                 C9 = all_markers$gene[which(all_markers$cluster == "9")],
                 C10 = all_markers$gene[which(all_markers$cluster == "10")],
                 C11 = all_markers$gene[which(all_markers$cluster == "11")],
                 C12 = all_markers$gene[which(all_markers$cluster == "12")],
                 C13 = all_markers$gene[which(all_markers$cluster == "13")])

#sig_dge.up <- subset(all_markers, p_val_adj<0.05&avg_log2FC>0.15)
#sig_dge.up <- sig_dge.up[order(sig_dge.up$avg_log2FC,decreasing = T),]
#sig_dge.up_TOP30 <- rownames(sig_dge.up[1:30,])
#diffall <-c(sig_dge.up_TOP30, sig_dge.down_TOP30) 

#matrix <- AverageExpression(object = data_seurat_malignant,assays = 'RNA',slot = "scale.data")[[1]]
#matrix <- matrix[rownames(matrix)%in%diffall,]
#matrix[matrix>2]=2; matrix[matrix< -2]= -2
#p=pheatmap( matrix ,show_colnames =T,
#            show_rownames = T,
#            cluster_cols = T, cluster_row = T,
#            border_color = NA,
#            color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

##BP, CC和MF三种通路一起富集
ego_ALL <- enrichGO(gene = rownames(all_markers),
                    OrgDb = 'org.Hs.eg.db',
                    keyType = 'SYMBOL',
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)#设置为ALL时BP, CC, MF都计算
ego_ALL <- data.frame(ego_ALL)
write.csv(ego_ALL,"enrichGO_all.csv")

ego_BP_C0 <- enrichGO(gene = deg.list[[1]],
                    OrgDb = 'org.Hs.eg.db',
                    keyType = 'SYMBOL',
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

#画图
pdf(file = "./13go富集分析/C0.pdf",height = 6,width = 8)
barplot(ego_BP_C0, showCategory = 10) + ggtitle("barplot for Biological process")
#dotplot(C0_enrich, showCategory = 10, size = NULL,font.size = 10, title = "C0_EnrichmentGO")#气泡图
dev.off()
#
ego_BP <- ego_BP[order(ego_BP$p.adjust),]
ego_BP_top30 <- ego_BP[1 : 30,]
ggplot(data=ego_BP_top30, aes(x=Description,y=Count)) +
  geom_bar(stat="identity",width=0.8,fill='salmon1') +
  coord_flip() + xlab("GO term") + ylab("Num of Genes") +
  theme_bw()


###msigdbr富集分析
# GSEA 要求输入的是一个将基因按照log2FC进行排序排好序的列表
Markers_genelist <- all_markers$avg_log2FC
names(Markers_genelist)= rownames(all_markers)
Markers_genelist <- sort(Markers_genelist, decreasing = T)



#获取人类的hall mark 基因集
library(msigdbr)
h.human <- msigdbr(species = "Homo sapiens", category = "H")
h.names <- unique(h.human$gs_name)
h.sets <- vector("list", length = length(h.names))
names(h.sets) <- h.names
for (i in names(h.sets)) {
  h.sets[[i]] <- subset(h.human,gs_name == i, "gene_symbol") %>% unlist(.)
}

m_df = msigdbr(species = 'Homo sapiens' , category = "H")
mf_df= m_df %>% dplyr::select(gs_name,gene_symbol)
#
gsea.results <- GSEA(Markers_genelist, TERM2GENE = mf_df)
gsea.results

gseaplot(gsea.results, gsea.results@result$ID)

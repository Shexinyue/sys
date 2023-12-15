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

rm(list=ls())
foldChange = 1  #Fold change
padj = 0.05  #adjust p-value
setwd("E:/LZ/24080")  #modify work path
if(! dir.exists("02_DEG")){
  file.create("02_DEG")
}
setwd("02_DEG")
library(edgeR)
library(tidyverse)
###-----------train------------
rt=read.table("../01_TCGA/files/mRNA_count.txt", sep="\t", header=T, check.names=F)  #change to your own file name_
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rm(list=c("exp","rt"))
data=avereps(data)
data=data[rowMeans(data)>1, ]
# 标准化
data <- log2(data+1)
# 获取normal和case
# 转置,得到df
df <- as.data.frame(t(data))
group <- df %>% 
  mutate(group = ifelse(grepl("-11A-", rownames(df)),"normal","tumor"))
group <- data.frame(sample = rownames(group),group = group$group)
group$sample <- substring(group$sample, 1, 16)
write.table(group,file="../00_rawdata/group.txt",sep="\t",quote=F,col.names=T,row.names = F)

group <- group %>% 
  pull(group)
table(group)
# normal  tumor 
#   50     374 
# group=c(rep("normal", 47), rep("tumor", 382))  #modify the number of cancer and normal samples
rm(list=c("df"))

design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("normal","tumor"))

topTags(et)
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

write.table(diff,file="edgerOut.xls",sep="\t",quote=F)
# DEG
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
length(rownames(diffSig))
# 2234----0.585
# 1134----1
# 240---2
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)   #???????л???У?????ı???ֵ??normalizeExp.txt??

diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)         #????????????У?????ı???ֵ??diffmRNAExp.txt??









#火山图--------------------------------------
allDiff$change = as.factor(
  ifelse(allDiff$FDR <0.05 & abs(allDiff$logFC) > foldChange,
         ifelse(allDiff$logFC > foldChange ,'UP','DOWN'),'NOT')
)

dat_rep<-allDiff[rownames(allDiff)%in%
                   rownames(rbind(head(diffSig[order(diffSig$logFC,decreasing = T),],10),
                                  head(diffSig[order(diffSig$logFC,decreasing = F),],10))),]

#dat_rep<-DEG[rownames(DEG)%in%DEERG$symbol,]


library(ggplot2)
library(ggthemes)
library(RColorBrewer)
# library(Ipaper)
library(scales)
library(ggrepel)



volcano_plot<- ggplot(data = allDiff,
                      aes(x = logFC,
                          y = -log10(FDR),
                          color =change)) +
  scale_color_manual(values = c("#69bcce", "darkgray","#EF5B5B")) +
  scale_x_continuous(
    breaks = seq(-6, 4, by = 2))+  # 根据需要调整刻度
  #limits = c(-10, 10))  # 设置对称的横坐标范围 +
  # scale_y_continuous(trans = "log1p",
  #                    breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 2.4, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-foldChange,foldChange),
             lty = 4,
             col = "gray40",
             lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "gray40",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 4,
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log2(Fold Change)",
       y = "-log10 (adj.P.Val)")
volcano_plot
pdf('volcano.pdf',  w=8,h=6,family='Times')
volcano_plot
dev.off()
png('volcano.png',  w=8,h=6,family='Times', units='in',res=600)
volcano_plot
dev.off()







#-----------------------------热图-----------------------------------------
rm(list= c("normalizeExp"))
library(ComplexHeatmap)
library(tidyverse)
dat_rep1<-allDiff[rownames(allDiff)%in%
                    rownames(rbind(head(diffSig[order(diffSig$logFC,decreasing = T),],10),
                                   head(diffSig[order(diffSig$logFC,decreasing = F),],10))),]



x <- read.table("./normalizeExp.txt",sep="\t",header = T,row.names = 1,check.names = F)
x<- log2(x+1)
mat <- t(scale(t(x)))#归一化
# mat <- x
# mat <- 2^mat-1
# mat <- log2(mat+0.01)

mat[mat < (-2)] <- (-2)
mat[mat > 2] <- 2

# dat_expr1 <- t(dat_expr) %>% as.data.frame()
# dat_expr1$sample <- rownames(dat_expr1)
# dat_expr1 <- dat_expr1[,c('sample','OLR1')]
#group_list <- read.csv(file = "../00_rawdata/02GSE51588/02normalize/GSE51588.group.txt",
#sep="\t")
group1 <- read.csv("../00_rawdata/group.txt",sep="\t")

colnames(group1) <- c("sample","group")

# group1 <- merge(group1, dat_expr1,by='sample')
# groupT <- group1[group1$type=='PD',]
# groupT <- group1[order(group1$OLR1,decreasing = T),]
# 
# group1 <- group1[order(group1$OLR1,decreasing = T),]
group1 <- group1[order(group1$group,decreasing = T),]

#group1 <- rbind(group1[1:13, ], group1[33, ],group1[14:32, ])
mat <- as.data.frame(mat)
colnames(mat) <- substring(colnames(mat),1,16)
mat <- as.matrix(mat)
dat_rep2 <- dat_rep1[order(dat_rep1$logFC,decreasing = T),]
mat <- mat[,group1$sample]
mat <- mat[rownames(dat_rep2),]

mat <- mat %>% as.data.frame()


png('03.heatmap.png',w=6,h=8,units='in',res=600,family='Times')
# 将数据框转换为矩阵
mat <- as.matrix(mat)
heatmap_plot <- 
  HeatmapAnnotation(Group = group1$group, col = list(Group = c("tumor" = "#f94141", "normal" = "#589ff3"))) %v%
  Heatmap(mat,
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression",
          ###cluster_cols = F,
          #cluster_rows = F,
          height = unit(12, "cm"),
          cluster_columns = F,
          cluster_rows = T,
          col = colorRampPalette(c("#6b98c4", "white","#c72228"))(100))
heatmap_plot
dev.off()

pdf('03.heatmap.pdf',w=6,h=8)
# 将数据框转换为矩阵
mat <- as.matrix(mat)
heatmap_plot <- 
  HeatmapAnnotation(Group = group1$group, col = list(Group = c("tumor" = "#f94141", "normal" = "#589ff3"))) %v%
  Heatmap(mat,
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression",
          ###cluster_cols = F,
          #cluster_rows = F,
          height = unit(12, "cm"),
          cluster_columns = F,
          cluster_rows = T,
          col = colorRampPalette(c("#6b98c4", "white","#c72228"))(100))
heatmap_plot
dev.off()




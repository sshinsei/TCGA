rm(list=ls())
setwd("G:/LZ/24080")
if(!dir.exists("05_clustDEG")){dir.create("05_clustDEG")}
setwd("./05_clustDEG")

library(limma)
library(ggplot2)
library(dplyr)
library(ggrepel)


# 加载表达矩阵数据
exp <- read.csv(file ="../02_DEG/normalizeExp.txt",
                sep="\t",row.names=1)
colnames(exp) <- gsub("\\.","-",colnames(exp))
colnames(exp) <- substring(colnames(exp),1,16)
# 加载分组变量
group_list <- read.csv(file = "../03_clust/cluster.txt",
                       sep="\t",header = F)
colnames(group_list) <- c("sample","group")
group_list$group <- paste0("clust",group_list$group)
table(group_list$group)
# clust1 clust2 
#    220    154

# 取样本子集
exp <- exp[,colnames(exp) %in% group_list$sample]

# 转置表达矩阵
dat_expr <- exp
dat_expr <- t(dat_expr)
# rownames(dat_expr) <- gsub(".*_","",rownames(dat_expr))



# group_list匹配表达矩阵的顺序
group <- group_list[match(rownames(dat_expr), group_list$sample), ]
group <- group %>% 
 pull("group")
# group

# 分组矩阵
design <- model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))
rownames(design) = rownames(dat_expr)
design

# 比较矩阵
contrast.matrix <- makeContrasts(paste0(unique(group),collapse = "-"),
                                 levels = design)
contrast.matrix

# 调整差异比较的方向，clust1为-1 :clust2相对于clust1上下调
contrast.mat = makeContrasts(contrasts="clust2-clust1", levels=design) #调整差异比较的方向，Normal为-1
contrast.mat



# 差异分析

# 转置
dat_expr <- t(dat_expr)

fit = lmFit(dat_expr, design)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")


DEG=na.omit(fit)

# logFC设置为1
logFC_cutoff <- 1
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val <0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
#[1] 17491     7 

dim(sig_diff)
#[1]646    7 


summary(sig_diff$change)
# DOWN  NOT   UP 
#  283   0    363



# 保存数据--------------------------#
DEG_write <- DEG[,-c(2,3,6)]
write.csv(DEG_write,file = "DEG_all.csv")
sig_diff_write <- sig_diff[,-c(2,3,6)]
write.csv(sig_diff_write, file = "DEG_sig.csv")


#火山图--------------------------------------


dat_rep<-DEG[rownames(DEG)%in%
               rownames(rbind(head(sig_diff[order(sig_diff$logFC,decreasing = T),],10),
                              head(sig_diff[order(sig_diff$logFC,decreasing = F),],10))),]

#dat_rep<-DEG[rownames(DEG)%in%DEERG$symbol,]


library(ggplot2)
library(ggthemes)
library(RColorBrewer)
# library(Ipaper)
library(scales)
library(ggrepel)

volcano_plot<- ggplot(data = DEG,
                      aes(x = logFC,
                          y = -log10(adj.P.Val),
                          color =change)) +
  scale_color_manual(values = c("#69bcce", "darkgray","#fed53b")) +
  scale_x_continuous(
    breaks = seq(-6, 4, by = 2))+  # 根据需要调整刻度
  #limits = c(-10, 10))  # 设置对称的横坐标范围 +
  # scale_y_continuous(trans = "log1p",
  #                    breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 2.4, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff),
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

library(ComplexHeatmap)
library(tidyverse)
dat_rep1<-DEG[rownames(DEG)%in%
                rownames(rbind(head(sig_diff[order(sig_diff$logFC,decreasing = T),],10),
                               head(sig_diff[order(sig_diff$logFC,decreasing = F),],10))),]



x<-dat_expr
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
group1 <- group_list
# group1 <- merge(group1, dat_expr1,by='sample')
# groupT <- group1[group1$type=='PD',]
# groupT <- group1[order(group1$OLR1,decreasing = T),]
# 
# group1 <- group1[order(group1$OLR1,decreasing = T),]
group1 <- group1[order(group1$group,decreasing = T),]

#group1 <- rbind(group1[1:13, ], group1[33, ],group1[14:32, ])

dat_rep2 <- dat_rep1[order(dat_rep1$logFC,decreasing = T),]
mat <- mat[,group1$sample]
mat <- mat[rownames(dat_rep2),]

mat <- mat %>% as.data.frame()


png('03.heatmap_clust.png',w=6,h=8,units='in',res=600,family='Times')
# 将数据框转换为矩阵
mat <- as.matrix(mat)
heatmap_plot <- 
  HeatmapAnnotation(Group = group1$group, col = list(Group = c("clust1" = "#FF0000", "clust2" = "#0000FF"))) %v%
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
          col = colorRampPalette(c("#69bcce", "white","#fed53b"))(100))
heatmap_plot
dev.off()

pdf('03.heatmap_clust.pdf',w=6,h=8)
# 将数据框转换为矩阵
mat <- as.matrix(mat)
heatmap_plot <- 
  HeatmapAnnotation(Group = group1$group, col = list(Group = c("clust1" = "#FF0000", "clust2" = "#0000FF"))) %v%
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
          col = colorRampPalette(c("#69bcce", "white","#fed53b"))(100))
heatmap_plot
dev.off()







library(ConsensusClusterPlus)
setwd("E:/LZ/24080") 
if(!dir.exists("03_clust")){dir.create("03_clust")}
setwd("03_clust")
rm(list=ls())
library(tidyverse)
# 表达矩阵
norExp = read.csv(file ="../02_DEG/normalizeExp.txt",
                  sep="\t")
colnames(norExp)[2:ncol(norExp)] <- gsub("\\.","-",colnames(norExp)[2:ncol(norExp)])
colnames(norExp)[2:ncol(norExp)] <- substring(colnames(norExp)[2:ncol(norExp)],1,16)
# 分组信息
group <- read.csv(file = "../00_rawdata/group.txt",
                  sep="\t")
#delete normal sample
rt <- norExp %>% column_to_rownames("id")
mask <- group[group$group == "tumor",]$sample
rt <- rt[ , colnames(rt) %in% mask]  

# hub_gene
genes <- read.csv("../06_multicox/multiCox.csv", header = T)
genes_1 = genes[ , 1]

rt <- rt %>% rownames_to_column("id")
rows = na.omit(match(genes_1, rt$id)) # PA:少一个基因
filter <- rt[rows, ]
unique(filter$id == genes[ , 1])  #check
rownames(filter) <- NULL
filter <- filter[ , -1]
#filter <- scale(filter)
filter = as.matrix(filter)
results = ConsensusClusterPlus(filter,
                               maxK=10,
                               reps=1000,
                               pItem=0.8,
                               pFeature=1,
                               title= "./clusterResult",
                               clusterAlg="km",
                               distance="euclidean",
                               seed=1234,
                               plot="pdf")

clusterNum = 2
cluster = results[[clusterNum]][["consensusClass"]]
write.table(cluster, file="cluster.txt", sep="\t", quote=F, col.names=F)



################################################################################
#################################  clust PCA  #########################################
################################################################################

rm(list = ls())
setwd("E:/LZ/24080") 
if(!dir.exists("03_clust")){dir.create("03_clust")}
setwd("03_clust")


hub_gene <- read.csv("../06_multicox/multiCox.csv", header = T)
exp <- read.csv(file ="../02_DEG/normalizeExp.txt",
                sep="\t",row.names = 1)
colnames(exp) <- gsub("\\.","-",colnames(exp))
colnames(exp) <- substring(colnames(exp),1,16)

group <- read.csv(file = "./cluster.txt",
                  sep="\t",header = F)
colnames(group) <- c('sample', 'group')

group$group <- paste0("clust",group$group)
# group$group <- factor(group$group, levels = c("clust2", "clust1"))

hubgene <-  hub_gene$id
expr <- exp

# 取样本子集
expr <- expr[,colnames(expr) %in% group$sample]


library(dplyr)
group <- group %>%dplyr::select(sample, group)
condition <- group


hubgene_expr<-expr[hubgene,]%>%t(.)%>%data.frame(.)

train_Group<-data.frame(sample=condition$sample,Group=condition$group) 


#Sample grouping
library(FactoMineR)
library(factoextra)
library(tinyarray)
library(tidyverse)



type = factor(x = rownames(hubgene_expr),labels = train_Group$Group)
summary(type)
#clust1 clust2 
# 220    154 

pdf(file = "PCA.pdf",width=5.5,height = 4.5)
draw_pca(t(hubgene_expr),type) +
  ggplot2::theme_bw()+
  ggplot2::theme(legend.position='top',
                 panel.grid = element_blank(),
                 axis.title.x = element_text(size=18),
                 axis.title.y = element_text(size=18),
                 axis.text.x = element_text(size = 18,color ='black'),
                 axis.text.y = element_text(size=18,color ='black'),#坐标轴数字黑色加粗
                 legend.text = element_text(size = 15),
                 #axis.line = element_line(colour = "black",linewidth =1.0),#坐标轴1.0磅，黑色
                 axis.ticks.length.x.bottom = unit (0.2, "cm"),#修改坐标轴刻度线长度
                 axis.ticks.length.y.left = unit(0.2,'cm'),
                 #axis.ticks = element_line(colour = "black",linewidth =1.0)
  )
dev.off()


png(file = "PCA.png",width=5.5,height = 4.5,units = "in",
    bg = "white", res = 500)
draw_pca(t(hubgene_expr),type) +
  ggplot2::theme_bw()+
  ggplot2::theme(legend.position='top',
                 panel.grid = element_blank(),
                 axis.title.x = element_text(size=18),
                 axis.title.y = element_text(size=18),
                 axis.text.x = element_text(size = 18,color ='black'),
                 axis.text.y = element_text(size=18,color ='black'),#坐标轴数字黑色加粗
                 legend.text = element_text(size = 15),
                 #axis.line = element_line(colour = "black",linewidth =1.0),#坐标轴1.0磅，黑色
                 axis.ticks.length.x.bottom = unit (0.2, "cm"),#修改坐标轴刻度线长度
                 axis.ticks.length.y.left = unit(0.2,'cm'),
                 #axis.ticks = element_line(colour = "black",linewidth =1.0)
  )
dev.off()


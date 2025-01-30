library(ConsensusClusterPlus)
setwd("E:/LZ/24105") 
if(!dir.exists("02_clust")){dir.create("02_clust")}
setwd("02_clust")
rm(list=ls())
library(tidyverse)
# 表达矩阵
norExp = read.csv(file ="../00_rawdata/GSE14520/GSE14520.txt",
                  sep="\t")
colnames(norExp)[2:ncol(norExp)] <- gsub("\\.","-",colnames(norExp)[2:ncol(norExp)])
colnames(norExp)[2:ncol(norExp)] <- substring(colnames(norExp)[2:ncol(norExp)],1,16)
# 分组信息
group <- read.csv(file = "../00_rawdata/GSE14520/GSE14520.group2.txt",
                  sep="\t")

#delete normal sample
rt <- norExp # %>% column_to_rownames("id")
mask <- group[group$group == "Case",]$sample # 247
rt <- rt[ , colnames(rt) %in% mask]  

# hub_gene
genes <- read.table(file = "../03_upset/feature_genes.txt", 
                    sep = "\t", header = F, check.names = F)
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

clusterNum = 3
cluster = results[[clusterNum]][["consensusClass"]]
write.table(cluster, file="cluster.txt", sep="\t", quote=F, col.names=F)


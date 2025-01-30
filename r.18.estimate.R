rm(list=ls())
library(limma)
library(estimate)
setwd('E:/LZ/24080')
if (! dir.exists('./18_estimate')){
  dir.create('./18_estimate')
}
setwd('./18_estimate')
rt=read.csv(file = "../02_DEG/normalizeExp.txt", sep="\t", 
            header=T, check.names=F,row.names = 1) 
#colnames(rt) = rt[1, ]
#rt = rt[-c(1:3), ]
#如果一个基因占了多行，取均值ֵ
#rt = rt[ , -c(2:52)]  #delete normal sample
#delete normal sample
colnames(rt) <- substring(colnames(rt),1,16)
rt <- rt[,!grepl("11A",colnames(rt))]

rt = as.matrix(rt)
exp = rt
rm(list = c("rt"))
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)
rm(list = c("exp"))

#Remove normal samples and keep only tumor samples
#group = sapply(strsplit(colnames(data), "\\-"), "[", 4)
#group = sapply(strsplit(group,""), "[", 1)
#group = gsub("2", "1", group)
#data = data[ , group == 0]
out = data[rowMeans(data) > 0, ]
out = rbind(ID = colnames(out), out)

#data_1 = data[ , -c(1:8)]  #remove normal sample
#out = data_1[rowMeans(data_1)>0, ]
#out = rbind(ID = colnames(out), out)
#???????????ľ????ļ?
write.table(out, file = "uniq.symbol.txt", sep = "\t", quote = F, col.names = F)

#????estimate??
filterCommonGenes(input.f = "uniq.symbol.txt",
                  output.f = "commonGenes.gct", id = "GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds = "estimateScore.gct")

#????ÿ????Ʒ?Ĵ???
scores = read.table("estimateScore.gct", skip = 2, header = T)
rownames(scores) = scores[ , 1]
scores = t(scores[ , 3:ncol(scores)])
rownames(scores) = gsub("\\.", "\\-", rownames(scores))
out = rbind(ID = colnames(scores), scores)
write.table(out, file = "scores.txt", sep = "\t", quote = F, col.names = F)


library(tidyverse) 
library(dplyr)
df <- read.csv("../00_rawdata/risk.txt",sep="\t")
cluster <- df %>% 
  dplyr::select("id","risk") %>% 
  mutate(clust = ifelse(risk == "low", "Cluster1","Cluster2"),
         risk = ifelse(risk == "low", "risk_low","risk_high")) %>% 
  arrange(clust)
cluster <- cluster[,c("id","clust","risk")]
write.table(cluster, file = "cluster.txt", sep = "\t", 
            quote = F, col.names =F,row.names = F)

################################## violinplot ####################################
rm(list=ls())
library(ggpubr)
Type = read.table("./cluster.txt", 
                  sep = "\t", check.names = F,
                  header = F) # 369
#Type=Type[order(Type[,2]),]

score = read.table("./scores.txt", 
                   sep = "\t", check.names = F,
                   header = T) # 382

#rows = na.omit(match(rownames(Type), rownames(score)))
#score = score[rows, ]
#unique(row.names(Type) == row.names(score))  #check
#colnames(Type) = c("cluster", "Subtype")
#rows = match(rownames(score), rownames(Type))
#cluster = cbind(Type[rows, ], score)
cluster = merge.data.frame(Type, score, by.x = "V1", by.y = "ID") # 326
#cluster = cluster[ , -1]
cluster = cluster[order(cluster$V3, decreasing = T), ]
colnames(cluster)[1:3] = c("id", "cluster", "Subtype")
#cluster$Subtype = factor(cluster$Subtype, levels = c("Subtype_1", "Subtype_2", "Subtype_3"))
#my_comparisons=list(c("Subtype_1","Subtype_2"),c("Subtype_2","Subtype_3"),
#                    c("Subtype_3","Subtype_1"))
#library(ggsci)
#pdf(file="TumorPurity_vioplot.pdf", width=6, height=5, )
#p <- ggviolin(cluster, x="Subtype", y="ESTIMATEScore", fill = "Subtype", 
#              palette = c("#E64B35FF", "#0072B5FF", "#00A087FF"), 
#              add = "boxplot", add.params = list(fill="white"))+
#  stat_compare_means(comparisons = my_comparisons,
#                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
#                                      symbols = c("***", "**", "*", "ns")),
#                     label = "p.signif")
#p
#ggsave(filename = "ESTIMATEScore_vioplot.tiff", plot = p, dpi = 200, width = 5.5,
#       height = 4.5)

#cluster$Subtype <- factor(cluster$Subtype, levels = c("Low", "High"))
#setwd("D:\\jxhe\\2021_11\\liver\\03ssGSEA\\ssGSEA_FPKM\\07estimateVioplot")

#--------------------------------ESTIMATEScore.tiff-----------------------------------
p <- ggviolin(cluster, x = "Subtype", y = "ESTIMATEScore", fill = "Subtype",
              palette = c("#4DBBD5FF", "#DC0000FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("ESTIMATEScore")+
  stat_compare_means(label = "p.format",label.x = 1.4)  #????ͼƬ?޸?
p
ggsave(filename = "ESTIMATEScore.tiff", plot = p, dpi = 400, width = 6, height = 5)
ggsave(filename = "ESTIMATEScore.pdf", plot = p, dpi = 500, width = 6, height = 5)


#--------------------------------StromalScore-----------------------------------
p <- ggviolin(cluster, x = "Subtype", y = "StromalScore", fill = "Subtype",
              palette = c("#4DBBD5FF", "#DC0000FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("StromalScore")+
  stat_compare_means(label = "p.format",label.x = 1.4)  #????ͼƬ?޸?
p
ggsave(filename = "StromalScore.tiff", plot = p, dpi = 400, width = 6, height = 5)
ggsave(filename = "StromalScore.pdf", plot = p, dpi = 500, width = 6, height = 5)


#--------------------------------ImmuneScore-----------------------------------
p <- ggviolin(cluster, x = "Subtype", y = "ImmuneScore", fill = "Subtype",
              palette = c("#4DBBD5FF", "#DC0000FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("ImmuneScore")+
  stat_compare_means(label = "p.format",label.x = 1.4)  #????ͼƬ?޸?
p
ggsave(filename = "ImmuneScore.tiff", plot = p, dpi = 400, width = 6, height = 5)
ggsave(filename = "ImmuneScore.pdf", plot = p, dpi = 500, width = 6, height = 5)


#--------------------------------TumorPurity-----------------------------------
p <- ggviolin(cluster, x = "Subtype", y = "TumorPurity", fill = "Subtype",
              palette = c("#4DBBD5FF", "#DC0000FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("TumorPurity")+
  stat_compare_means(label = "p.format",label.x = 1.4)  #????ͼƬ?޸?
p
ggsave(filename = "TumorPurity.tiff", plot = p, dpi = 400, width = 6, height = 5)
ggsave(filename = "TumorPurity.pdf", plot = p, dpi = 500, width = 6, height = 5)










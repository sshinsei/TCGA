rm(list=ls())
setwd("E:/LZ/24105/") 
if(! dir.exists("13_immu")) dir.create("13_immu")
setwd("13_immu")


#################################################################################
#########################    file prepare   #####################################
#################################################################################
library(GSVA)
library(limma)
library(GSEABase)
rt=read.csv(file = "../00_rawdata/normalizeExp.txt", sep="\t", 
            header=T, check.names=F)
# 去掉正常样本：
# colnames(rt)[2:ncol(rt)] <- substring(colnames(rt)[2:ncol(rt)],1,16)
rt <- rt[,!grepl("11A",colnames(rt))]

#rt = rt[ , -c(2:52)]
rt = as.matrix(rt)
#rt = t(rt)
rownames(rt) = rt[ , 1]
rt = rt[-1, ]
#rownames(rt)=rt[,1]
exp=rt[ , 2:ncol(rt)]
rm(list=c("rt"))
dimnames=list(rownames(exp), colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0, ]
gmtFile="immune.gmt"    #GMT file
geneSet=getGmt(gmtFile, geneIdType = SymbolIdentifier())

# 可能会报matrixStats 中 useNames = NA的错->解决：安装旧版本的matrixStats
ssgsea_par <- ssgseaParam(mat, geneSet)  # all other values are default values
ssgseaScore <- gsva(ssgsea_par)

# ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

# 标准化ssGSEA score
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

ssgseaOut = normalize(ssgseaScore)
ssgseaOut = rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut, file = "ssgseaOut.txt", sep = "\t", quote = F, col.names = F)



# rm(list=ls())
library(limma)
library(estimate)
# setwd("E:/01_ccRCCs/11risk_immu/02estimate")  

data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)
# rm(list = c("exp"))

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

# make clust
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




#################################################################################
#########################    plot heatmap   #####################################
#################################################################################
rm(list=ls())
library(pheatmap)  
setwd("E:/LZ/24105/13_immu/")   
rt=read.table("./ssgseaOut.txt",
              sep="\t",header=T,row.names=1,check.names=F) # 382

#colnames(rt) = substring(colnames(rt), 1, 12)

Type=read.table("cluster.txt",sep="\t",check.names=F,header=F)
#cols = na.omit(match(row.names(Type), colnames(rt)))
#rt = rt[ , cols]

score=read.table("./scores.txt",sep="\t",check.names=F, header=T)
#rownames(score) = substring(rownames(score), 1, 12)

#rows = na.omit(match(row.names(Type), rownames(score)))
#score = score[rows, ]
#colnames(Type) = c("cluster", "Subtype")
#rownames(score) == rownames(Type) #check
#rows = na.omit(match(row.names(score), rownames(Type)))
cluster = merge.data.frame(Type, score, by.x = "V1", by.y = "ID")
cluster = cluster[order(cluster$V2), ]

if(F){
  # -----------rt的列名-------------------------- #
  colnames(rt) <- gsub("11A.*$", "11A", colnames(rt))
  colnames(rt) <- gsub("01A.*$", "01A", colnames(rt))
  colnames(rt) <- gsub("01B.*$", "01B", colnames(rt))
  rt = rt[ , cluster[ , 1]]
  # 找出有重复的列名，重复的列名会加.1
  temp <- colnames(rt[,cluster$V1 != colnames(rt)])
  # 删除重复的列
  rt <- rt[,! colnames(rt) %in% temp]
  # --------------------------------------------- #
}
colnames(rt) <- gsub("\\.","-",colnames(rt))
rt = rt[ , cluster[ , 1]]
unique(cluster$V1 == colnames(rt))  #check
rownames(cluster) = cluster$V1
cluster = cluster[ , -c(1, 2)]
colnames(cluster)[1] = "Subtype"
#??????ͼ
library(ggsci)
ann_color <- list(Subtype = c(risk_low = "#4DBBD5FF", risk_high = "#DC0000FF"))
tiff(file="estimateHM.tiff",width = 8, height = 5, units = "in", bg = "white",
     res = 400)
pheatmap(as.matrix(rt), annotation_col = cluster, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3,
         annotation_colors = ann_color)
dev.off()

pdf("estimateHM.pdf", height = 5, width = 8)
#annotation_col = data.frame(Subtype=factor(c(rep("Subtype_1",33),
#rep("Subtype_2",14),rep("Subtype_3",2))))
#row.names(annotation_col) = rownames(cluster)
pheatmap(as.matrix(rt), annotation_col=cluster, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3,
         annotation_colors = ann_color)
dev.off()


#################################################################################
#########################    plot violin plot   #################################
#################################################################################
rm(list=ls())
setwd("E:/LZ/24105/13_immu/")   
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


#################################################################################
#########################    plot checkpoint    #################################
#################################################################################
rm(list=ls())
setwd("E:/LZ/24105/13_immu/")
imm_checkp = read.table("E:/LZ/24080/20_TIDE/check_point_gene.txt",sep="\t",header=F,check.names=F)
rt = read.csv(file = "../00_rawdata/normalizeExp.txt", sep="\t", 
              header=T, check.names=F, row.names = 1)


#rt = rt[ , -c(1:51)]  #delete normal sample

# 去掉正常样本：
colnames(rt) <- substring(colnames(rt),1,16)
rt <- rt[,!grepl("11A",colnames(rt)[2:ncol(rt)])]

rows = na.omit(match(imm_checkp$V1, rownames(rt)))
rt_imm_checkp = rt[rows, ]

write.table(rt_imm_checkp, file="Immune_checkpoint_exp.txt", sep="\t",
            quote=F, col.names=T, row.names=T)


library(ggpubr)
library(tidyverse)

rt=rt_imm_checkp
#colnames(rt) = substring(colnames(rt), 1, 12)
Type=read.table("./cluster.txt", 
                sep = "\t", check.names = F, row.names = 1,
                header = F)
#Type=Type[order(Type[ , 2]), ]
# 处理rt的样本名
# getImmune中已经处理过样本，所以不用处理
if(F){
  # ----------------rt的列名-------------------------- #
  rt <- as.data.frame(t(rt))
  rt <- rt %>% 
    rownames_to_column("id")
  
  rt$id <- gsub("11A.*$", "11A", rt$id)
  rt$id <- gsub("01A.*$", "01A", rt$id)
  rt$id <- gsub("01B.*$", "01B", rt$id)
  
  # 去重
  rt <- rt %>% 
    distinct(id, .keep_all = TRUE) %>% 
    column_to_rownames("id")
  
  # 判断列名是否都唯一
  !any(duplicated(colnames(rt)))
  
  # 转置
  rt <- as.data.frame(t(rt))
  # --------------------------------------------- #
}


# 先取type的子集
Type <- Type[rownames(Type) %in% colnames(rt),]

cols = na.omit(match(rownames(Type), colnames(rt)))
rt=t(rt[ , cols])
unique(rownames(rt) == rownames(Type)) #check

#׼??????ͼ???????ļ?
data = data.frame()
rows = na.omit(match(rownames(rt), rownames(Type)))
for(i in colnames(rt)){
  a = as.data.frame(cbind(expression=rt[ , i], gene=i, Subtype=as.vector(Type[rows, 2])))
  a[ , 1] = as.numeric(a[ , 1])
  test <- wilcox.test(expression ~ Subtype, data = a)
  pValue = test$p.value
  if(!is.na(pValue) & pValue < 0.01){
    data = rbind(data, a)
    print(pValue)
  }
  #data=rbind(data,cbind(expression=rt[ , i], gene=i, Subtype=as.vector(Type[,2])))
}
write.table(data, file="data.txt",sep="\t",row.names=F,quote=F)

#????????ͼ
data=read.table("data.txt",sep="\t",header=T,check.names=F)  
data$expression = log2(data$expression + 0.01)
data$Subtype = factor(data$Subtype, levels = c("risk_low", "risk_high"))
p = ggboxplot(data, x="gene", y="expression", fill = "Subtype",
              ylab = "Gene expression",
              xlab = "", palette = c("#4DBBD5FF", "#DC0000FF"))+
  theme(text = element_text(size = 21))
p = p+rotate_x_text(60)
p = p+stat_compare_means(aes(group=Subtype),
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns")),
                         label = "p.signif")
p
ggsave(filename = "checkpoint.tiff", plot = p, dpi = 400, width = 17, height = 8)
ggsave(filename = "checkpoint.pdf", plot = p, dpi = 500, width = 17, height = 8)

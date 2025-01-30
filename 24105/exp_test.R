rm(list = ls())
setwd('E:/LZ/24105')
if(!dir.exists("./exp_test")){dir.create("./exp_test")}
setwd("./exp_test")
library(tidyverse)
##-------------------------------------train-----------------------------------------------
# 表达矩阵
norExp <- read.csv(file ="../00_rawdata/normalizeExp.txt",
                   sep="\t",row.names=1)
colnames(norExp) <- gsub("\\.","-",colnames(norExp))
# 分组矩阵
group <- read.csv(file = "../00_rawdata/group.txt",
                  sep="\t")
group <- group %>% mutate(
  group = if_else(group == "tumor","Case","Normal")
)
# 特征基因：待修改
#genes = c('TYRO3','SLC2A1','ARG2', 'SPHK1','EPO','GULP1')
#genes <- as.data.frame(genes)
genes <- read.csv("../03_upset/feature_genes.txt",header = F,sep="\t")
#rows = match(group$sample, rownames(norExp))
#data = norExp[rows, ]
#rownames(data) == genes$V1
data = data.frame()
for (i in 1:nrow(genes)) {
  gene = genes[i, 1]
  da = as.data.frame(t(norExp[gene, ]))
  da <- da %>% rownames_to_column("sample")
  # 构建分组
  da <- da %>% 
    left_join(.,group,by="sample") 
  # print(table(da$Group))
  da$gene = gene
  rownames(da) = NULL
  da <- da[,2:4]
  colnames(da) = c("expression", "Group", "gene")
  #da[ , 1] = log2(da[ , 1] + 0.01)
  #data <- data[,2:4]
  data = rbind(data, da)
}

#boxplot
library(ggpubr)
library(ggsci)
data$Group=factor(data$Group, levels = c("Normal", "Case"))
# data$expression = log2(data$expression + 1)

p=ggboxplot(data, x="gene", y="expression", color = "Group",
            ylab="Gene expression",
            xlab="",
            palette = c("#00A087B2", "#E64B35FF") )
p=p+rotate_x_text(60)
#pdf(file="boxplot_1.pdf", width=7.5, height=5.5)                          #????ͼƬ?ļ?
#p+stat_compare_means(aes(group=Group),
#                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
#                                      symbols = c("***", "**", "*", "ns")),
#                     label = "p.signif")
#dev.off()

p=p+stat_compare_means(aes(group=Group),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")
p + theme(text = element_text(size = 18))
ggsave(filename = "train.png", plot = p, dpi = 400, width = 10,
       height = 5)
ggsave(filename = "train.pdf", plot = p, dpi = 500, width = 10,
       height = 5)


##-------------------------------------test 1-----------------------------------------------
rm(list = ls())
setwd('E:/LZ/24105')
if(!dir.exists("./exp_test")){dir.create("./exp_test")}
setwd("./exp_test")
library(tidyverse)
# 表达矩阵
norExp = read.csv(file ="../00_rawdata/GSE14520/GSE14520.txt",
                       sep="\t",row.names=1)
# 分组矩阵
group <- read.csv(file = "../00_rawdata/GSE14520/GSE14520.group2.txt",
                  sep="\t")
#genes = c('TYRO3','SLC2A1','ARG2', 'SPHK1','EPO','GULP1')
#genes <- as.data.frame(genes)
genes <- read.csv("../03_upset/feature_genes.txt",header = F,sep="\t")
#genes <- read.csv("../00_rawdata/ERGs_2.txt",header = F,sep="\t")
#rows = match(group$sample, rownames(norExp))
#data = norExp[rows, ]
#rownames(data) == genes$V1
data = data.frame()
for (i in 1:nrow(genes)) {
  gene = genes[i, 1]
  da = as.data.frame(t(norExp[gene, ]))
  da <- da %>% rownames_to_column("sample")
  # 构建分组
  da <- da %>% 
    left_join(.,group,by="sample") 
  # print(table(da$Group))
  da$gene = gene
  rownames(da) = NULL
  da <- da[,2:4]
  colnames(da) = c("expression", "Group", "gene")
  #da[ , 1] = log2(da[ , 1] + 0.01)
  #data <- data[,2:4]
  data = rbind(data, da)
}

#boxplot
library(ggpubr)
library(ggsci)
data$Group=factor(data$Group, levels = c("Normal", "Case"))
#data$expression = log2(data$expression + 1)

p=ggboxplot(data, x="gene", y="expression", color = "Group",
            ylab="Gene expression",
            xlab="",
            palette = c("#00A087B2", "#E64B35FF") )
p=p+rotate_x_text(60)
#pdf(file="boxplot_1.pdf", width=7.5, height=5.5)                          #????ͼƬ?ļ?
#p+stat_compare_means(aes(group=Group),
#                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
#                                      symbols = c("***", "**", "*", "ns")),
#                     label = "p.signif")
#dev.off()

p=p+stat_compare_means(aes(group=Group),
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),
                       label = "p.signif")
p + theme(text = element_text(size = 18))
ggsave(filename = "GSE14520.png", plot = p, dpi = 400, width = 10,
       height = 5)
ggsave(filename = "GSE14520.pdf", plot = p, dpi = 500, width = 10,
       height = 5)


##-------------------------------------test 2-----------------------------------------------
rm(list = ls())
setwd('E:/LZ/24105')
if(!dir.exists("./exp_test")){dir.create("./exp_test")}
setwd("./exp_test")
library(tidyverse)
# 表达矩阵
norExp = read.csv(file ="../00_rawdata/GSE104580/02normalize/GSE104580.txt",
                  sep="\t",row.names=1)
# 分组矩阵
group <- read.csv(file = "../00_rawdata/GSE104580/02normalize/GSE104580.group.txt",
                  sep="\t")
#genes = c('TYRO3','SLC2A1','ARG2', 'SPHK1','EPO','GULP1')
#genes <- as.data.frame(genes)
#genes <- read.csv("../00_rawdata/genes.txt",header = F,sep="\t")
genes <- read.csv("../03_upset/feature_genes.txt",header = F,sep="\t")
#rows = match(group$sample, rownames(norExp))
#data = norExp[rows, ]
#rownames(data) == genes$V1
data = data.frame()
for (i in 1:nrow(genes)) {
  gene = genes[i, 1]
  da = as.data.frame(t(norExp[gene, ]))
  da <- da %>% rownames_to_column("sample")
  # 构建分组
  da <- da %>% 
    left_join(.,group,by="sample") 
  # print(table(da$Group))
  da$gene = gene
  rownames(da) = NULL
  da <- da[,2:4]
  colnames(da) = c("expression", "Group", "gene")
  #da[ , 1] = log2(da[ , 1] + 0.01)
  #data <- data[,2:4]
  data = rbind(data, da)
}

#boxplot
library(ggpubr)
library(ggsci)
data$Group=factor(data$Group, levels = c("Normal", "Case"))
# data$expression = log2(data$expression + 1)

p=ggboxplot(data, x="gene", y="expression", color = "Group",
            ylab="Gene expression",
            xlab="",
            palette = c("#00A087B2", "#E64B35FF") )
p=p+rotate_x_text(60)
#pdf(file="boxplot_1.pdf", width=7.5, height=5.5)                          #????ͼƬ?ļ?
#p+stat_compare_means(aes(group=Group),
#                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
#                                      symbols = c("***", "**", "*", "ns")),
#                     label = "p.signif")
#dev.off()

p=p+stat_compare_means(aes(group=Group),
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),
                       label = "p.signif")
p + theme(text = element_text(size = 18))
ggsave(filename = "GSE104580.png", plot = p, dpi = 400, width = 10,
       height = 5)
ggsave(filename = "GSE104580.pdf", plot = p, dpi = 500, width = 10,
       height = 5)



##-------------------------------------test 3-----------------------------------------------
rm(list = ls())
setwd('E:/LZ/24105')
if(!dir.exists("./exp_test")){dir.create("./exp_test")}
setwd("./exp_test")
library(tidyverse)
# 表达矩阵
norExp = read.csv(file ="../00_rawdata/GSE76427/GSE76427.txt",
                  sep="\t",row.names=1)
# 分组矩阵
group <- read.csv(file = "../00_rawdata/GSE76427/GSE76427.group.txt",
                  sep="\t")
# 特征基因：待修改
genes = c('TYRO3','SLC2A1','ARG2', 'SPHK1','EPO','GULP1')
genes <- as.data.frame(genes)
#genes <- read.csv("../04_venn/feature_genes.txt",header = F,sep="\t")
#rows = match(group$sample, rownames(norExp))
#data = norExp[rows, ]
#rownames(data) == genes$V1
data = data.frame()
for (i in 1:nrow(genes)) {
  gene = genes[i, 1]
  da = as.data.frame(t(norExp[gene, ]))
  da <- da %>% rownames_to_column("sample")
  # 构建分组
  da <- da %>% 
    left_join(.,group,by="sample") 
  # print(table(da$Group))
  da$gene = gene
  rownames(da) = NULL
  da <- da[,2:4]
  colnames(da) = c("expression", "Group", "gene")
  #da[ , 1] = log2(da[ , 1] + 0.01)
  #data <- data[,2:4]
  data = rbind(data, da)
}

#boxplot
library(ggpubr)
library(ggsci)
data$Group=factor(data$Group, levels = c("Normal", "Case"))
# data$expression = log2(data$expression + 1)

p=ggboxplot(data, x="gene", y="expression", color = "Group",
            ylab="Gene expression",
            xlab="",
            palette = c("#00A087B2", "#E64B35FF") )
p=p+rotate_x_text(60)
#pdf(file="boxplot_1.pdf", width=7.5, height=5.5)                          #????ͼƬ?ļ?
#p+stat_compare_means(aes(group=Group),
#                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
#                                      symbols = c("***", "**", "*", "ns")),
#                     label = "p.signif")
#dev.off()

p=p+stat_compare_means(aes(group=Group),
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),
                       label = "p.signif")
p + theme(text = element_text(size = 18))
ggsave(filename = "test3.png", plot = p, dpi = 400, width = 10,
       height = 5)
ggsave(filename = "test3.pdf", plot = p, dpi = 500, width = 10,
       height = 5)

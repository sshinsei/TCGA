
###########--NORMAL-CLUST1-CLUST2的表达量验证--#########################
rm(list = ls())
setwd('E:/LZ/24080')
if(!dir.exists("10_clustexp")){
  dir.create("10_clustexp")}
setwd("10_clustexp")


# 是否画图
rm(list=ls())
plot <- T #T/F

hub_gene <- read.csv('../06_multicox/multiCox.csv', header = T)

library(tidyverse)
train_data <- read.csv(file ="../02_DEG/normalizeExp.txt",
                       sep="\t",row.names = 1)
colnames(train_data) <- gsub("\\.","-",colnames(train_data))
colnames(train_data) <- substring(colnames(train_data),1,16)

group1 <- read.csv(file = "../03_clust/cluster.txt",
                        sep="\t",header = F)
colnames(group1) <- c('sample', 'group')
# train_group的group上paste clust
group1$group <- paste0("clust",group1$group)
group1$group <- factor(group1$group, levels = c("clust1", "clust2"))

group2 <- read.csv(file = "../00_rawdata/group.txt",
                   sep="\t")
train_group <- left_join(group2,group1,by="sample")
train_group <- train_group %>% 
  mutate(
    group.x = case_when(
      group.y == "clust1" ~ "clust1",
      group.y == "clust2" ~ "clust2",
      TRUE ~ group.x   # 其他情况保持原值
    )) %>% 
  select(sample,group.x)
colnames(train_group) <- c("sample","group")

hubgene <-  hub_gene$id
expr <- train_data
# 取样本子集
expr <- expr[,colnames(expr) %in% train_group$sample]

train_group <- train_group
library(dplyr)
train_group <- train_group %>%dplyr::select(sample, group)
condition <- train_group

#hubgene<- c("FANCI","TRIP13","UBE2T","FEN1","CENPN")
hubgene_expr<-expr[hubgene,]%>%t(.)%>%data.frame(.)
# hubgene_expr<-hubgene_expr[,order(colnames(hubgene_expr),decreasing = F)]

train_Group<-data.frame(sample=condition$sample,Group=condition$group) 

##按照字母进行排序
dat <-hubgene_expr %>% as.data.frame %>% tibble::rownames_to_column(var = "sample")
identical(dat$sample,train_Group$sample)
dat <- merge(dat, train_Group, by = "sample")
dat2 <- tidyr::gather(dat, Gene, expression, -c(sample, Group))
dat2 <- na.omit(dat2)
library(rstatix)
stat_res <- dat2 %>% 
  group_by(Gene) %>% 
  wilcox_test(expression ~ Group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr,none
  add_significance("p")
stat_res  

# 筛选基因：只保留clust1-clust2之间显著的
plot_genes <- stat_res %>% 
  filter(group1 == "clust1" & group2 == "clust2" & p.adj.signif != "ns") %>% 
  pull(Gene)

# 画图数据准备
# hubgene_expr取子集
hubgene_expr <- hubgene_expr[,colnames(hubgene_expr) %in% plot_genes]
boxplot_dat <-data.frame(cbind(hubgene_expr,condition$group))
colnames(boxplot_dat)[ncol(boxplot_dat)]<-'Group'
head(boxplot_dat)
boxplot_dat$id <- rownames(boxplot_dat)
library(reshape2)
boxplot_dat <- melt(boxplot_dat,id = c("id","Group"))
boxplot_dat$value <-as.numeric(boxplot_dat$value)
boxplot_dat$Group<- factor(boxplot_dat$Group)
colnames(boxplot_dat)<-c('id','Group','Gene','value')

boxplot_dat<-boxplot_dat[order(boxplot_dat$Gene,decreasing = T),]

boxplot_dat$value <- log2(boxplot_dat$value + 1)

# 画图
library(ggsci)
library(ggplot2)
library(ggpubr)

##########-------------------PA-related基因的表达量--------------------------------------------
# 表达矩阵
norExp <- expr
# 分组矩阵
group <- condition

# 差异基因
temp <- read.table("../02_DEG/diffmRNAExp.txt", header = T, sep = "\t")
# 与DEG取交集e.tabel(
a <- list(temp$id,plot_genes)
genes <- Reduce(intersect,a) # 3 
genes <- as.data.frame(genes)
write.table(genes,"feature_genes.txt",col.names = F,row.names = F, quote= F,sep="\t")
#rows = match(group$sample, rownames(norExp))
#data = norExp[rows, ]
#rownames(data) == genes$V1
data = data.frame()
for (i in 1:nrow(genes)) {
  gene = genes[i, 1]
  da = as.data.frame(t(norExp[gene, ])) # norEXP是表达矩阵
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
data$Group=factor(data$Group)
# 标准化
data$expression = log2(data$expression + 1)

p=ggboxplot(data, x="gene", y="expression", color = "Group",
            ylab="Gene expression",
            xlab="",
            palette = c("#1f77b4", "#ff7f0e", "#2ca02c") )
p=p+rotate_x_text(60)

p=p+stat_compare_means(aes(group=Group),
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),
                       label = "p.signif")
p + theme(text = element_text(size = 18))
ggsave(filename = "PA_genes.png", plot = p, dpi = 400, width = 10,
       height = 5)
ggsave(filename = "PA_genes.pdf", plot = p, dpi = 500, width = 10,
       height = 5)
  



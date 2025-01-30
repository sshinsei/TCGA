rm(list = ls())
setwd('E:/LZ/24105')
if(!dir.exists("04_roc")){
  dir.create("04_roc")}
setwd("04_roc")

#  ROC曲线
#ROC----train---------

library(tidyverse)
hub_gene <- read.csv('../03_model/lasso_svm_xgb_gene.txt', header = F)

train_data <- read.csv(file ="../00_rawdata/normalizeExp.txt",
                       sep="\t",row.names=1)
colnames(train_data) <- gsub("\\.","-",colnames(train_data))
train_group <- read.csv(file = "../00_rawdata/group.txt",
                        sep="\t")
train_group <- train_group %>% mutate(
  group = if_else(group == "tumor","Case","Normal")
)
train_group$group <- factor(train_group$group, levels = c("Case", "Normal"))


##候选基因
#gene <- read.csv("../06_violin/DE-GSE165082-wilcox_result.csv",header = T,row.names = 1)
#train_data <- t(train_data)
data_candi <- train_data[rownames(train_data)%in%hub_gene$V1,train_group$sample]


exp <-data_candi

library(pROC)
library(ggplot2)
exp2 <- t(exp)
exp2 <- cbind(exp2,train_group)
## 绘制ROC曲线
library(pROC)
# corlor <- c("#FF2E63","#87CEFA")
corlor <- c("#FF2E63","#FF4500","firebrick","blue","#3CB371","gray50","#E7B800","#87CEFA","#CD5C5C","#F0E68C","#3CB371")
len <- length(colnames(exp2))-2
# png(file = 'train-ROC.all.png',width = 400,height = 400)
tiff(file="train-ROC.tiff",width = 5, height = 5, units = "in",
     bg = "white", res = 500)
for (i in c(1:len)) {
  roc<-roc(exp2$group,exp2[,i],levels=c("Case", "Normal"))
  # roc = smooth(roc)
  if (i ==1) {
    plot(roc,
         print.auc=T,
         print.auc.x=0.5,print.auc.y=0.5,
         print.auc.pattern=paste(colnames(exp2)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main="GSE51588 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }else{
    plot(roc,
         print.auc=T,
         add=T,
         print.auc.x=0.5,print.auc.y=(0.5-i*0.05+0.05),
         print.auc.pattern=paste(colnames(exp2)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main="GSE51588 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }
  i<-i+1
}
dev.off()


pdf(file = 'train-ROC.all.pdf',w=5,h=5)
for (i in c(1:len)) {
  roc<-roc(exp2$group,exp2[,i],levels=c("Case", "Normal"))
  # roc = smooth(roc)
  if (i ==1) {
    plot(roc,
         print.auc=T,
         print.auc.x=0.5,print.auc.y=0.5,
         print.auc.pattern=paste(colnames(exp2)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main="GSE51588 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }else{
    plot(roc,
         print.auc=T,
         add=T,
         print.auc.x=0.5,print.auc.y=(0.5-i*0.05+0.05),
         print.auc.pattern=paste(colnames(exp2)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main="train ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }
  i<-i+1
}
dev.off()



##------test 1-------------
rm(list = ls())
setwd('E:/LZ/24105')
if(!dir.exists("04_roc")){
  dir.create("04_roc")}
setwd("04_roc")
library(tidyverse)
#data <- read.csv('../00_Rowdata/dat(GSE14905).csv',row.names = 1)
#group <- read.csv('../00_Rowdata/group(GSE14905).csv',row.names = 1)
test_data <- read.csv(file ="../00_rawdata/GSE14520/GSE14520.txt",
                      row.names=1, sep="\t")
test_group <- read.csv(file = "../00_rawdata/GSE14520/GSE14520.group.txt",
                       sep="\t")
test_group$group <- factor(test_group$group, levels = c("Case", "Normal"))

##候选基因
#hub_gene <- read.csv('../05_model/13.svm_lasso.txt', header = F)
hub_gene <- read.csv('../03_model/lasso_svm_xgb_gene.txt', header = F)
#hub_gene <- c("IL1R2", "IRAK3", "SLA")
#hub_gene <- c("NPTX2","LOXL1","COL1A2","C3AR1") #LOXL1,COL1A2无
#hub_gene <- as.data.frame(hub_gene)
#colnames(hub_gene) <- c("V1") # P2RY8 缺少
#gene <- read.csv("../06_violin/DE-GSE165082-wilcox_result.csv",header = T,row.names = 1)
#data <- t(data)
data_candi_test <- test_data[rownames(test_data)%in%hub_gene$V1,test_group$sample]


exp_test <- data_candi_test

#library(pROC)
#library(ggplot2)
exp2_test <- t(exp_test)
exp2_test <- cbind(exp2_test,test_group)
## 绘制ROC曲线
library(pROC)
# corlor <- c("#FF2E63","#87CEFA")
corlor <- c("#FF2E63","#FF4500","firebrick","blue","#3CB371","gray50","#E7B800","#87CEFA","#CD5C5C","#F0E68C","#3CB371")
len <- length(colnames(exp2_test))-2
# png(file = 'test-ROC.all.png',width = 400,height = 400)
tiff(file="test1-ROC.tiff",width = 5, height = 5, units = "in",
     bg = "white", res = 500)
for (i in c(1:len)) {
  roc<-roc(exp2_test$group,exp2_test[,i],levels=c("Case","Normal"))
  # 报错Error in smooth_roc_binormal(roc, n) : ROC curve not smoothable (not enough points).则把roc = smooth(roc)注释#
  # roc = smooth(roc)
  if (i ==1) {
    plot(roc,
         print.auc=T,
         print.auc.x=0.5,print.auc.y=0.5,
         print.auc.pattern=paste(colnames(exp2_test)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main=" GSE12021 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }else{
    plot(roc,
         print.auc=T,
         add=T,
         print.auc.x=0.5,print.auc.y=(0.5-i*0.05+0.05),
         print.auc.pattern=paste(colnames(exp2_test)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main="GSE12021 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }
  i<-i+1
}
dev.off()


pdf(file = 'test1-ROC.all.pdf',w=5,h=5)
for (i in c(1:len)) {
  roc<-roc(exp2_test$group,exp2_test[,i],levels=c("Case","Normal"))
  # roc = smooth(roc)
  if (i ==1) {
    plot(roc,
         print.auc=T,
         print.auc.x=0.5,print.auc.y=0.5,
         print.auc.pattern=paste(colnames(exp2_test)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main="GSE12021 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }else{
    plot(roc,
         print.auc=T,
         add=T,
         print.auc.x=0.5,print.auc.y=(0.5-i*0.05+0.05),
         print.auc.pattern=paste(colnames(exp2_test)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main="GSE12021 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }
  i<-i+1
}
dev.off()

##------test 3-------------
rm(list = ls())
setwd('E:/LZ/24105')
if(!dir.exists("04_roc")){
  dir.create("04_roc")}
setwd("04_roc")
library(tidyverse)
#data <- read.csv('../00_Rowdata/dat(GSE14905).csv',row.names = 1)
#group <- read.csv('../00_Rowdata/group(GSE14905).csv',row.names = 1)
test_data <- read.csv(file ="../00_rawdata/GSE104580/val.combat.exp.txt",
                      row.names=1, sep="\t")
test_group <- read.csv(file = "../00_rawdata/06val_combat/val.group.txt",
                       sep="\t")
test_group$group <- factor(test_group$group, levels = c("Case", "Normal"))

##候选基因
hub_gene <- read.csv('../05_model/13.svm_lasso.txt', header = F)
#hub_gene <- read.csv('../05_model/lasso_svm_rf_gene.txt', header = F)
#hub_gene <- c("IL1R2", "IRAK3", "SLA")
#hub_gene <- c("NPTX2","LOXL1","COL1A2","C3AR1") #LOXL1,COL1A2无
#hub_gene <- as.data.frame(hub_gene)
#colnames(hub_gene) <- c("V1") # P2RY8 缺少
#gene <- read.csv("../06_violin/DE-GSE165082-wilcox_result.csv",header = T,row.names = 1)
#data <- t(data)
data_candi_test <- test_data[rownames(test_data)%in%hub_gene$V1,test_group$sample]


exp_test <- data_candi_test

#library(pROC)
#library(ggplot2)
exp2_test <- t(exp_test)
exp2_test <- cbind(exp2_test,test_group)
## 绘制ROC曲线
library(pROC)
# corlor <- c("#FF2E63","#87CEFA")
corlor <- c("#FF2E63","#FF4500","firebrick","blue","#3CB371","gray50","#E7B800","#87CEFA","#CD5C5C","#F0E68C","#3CB371")
len <- length(colnames(exp2_test))-2
# png(file = 'test-ROC.all.png',width = 400,height = 400)
tiff(file="test3-ROC.tiff",width = 5, height = 5, units = "in",
     bg = "white", res = 500)
for (i in c(1:len)) {
  roc<-roc(exp2_test$group,exp2_test[,i],levels=c("Case","Normal"))
  # 报错Error in smooth_roc_binormal(roc, n) : ROC curve not smoothable (not enough points).则把roc = smooth(roc)注释#
  # roc = smooth(roc)
  if (i ==1) {
    plot(roc,
         print.auc=T,
         print.auc.x=0.5,print.auc.y=0.5,
         print.auc.pattern=paste(colnames(exp2_test)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main=" GSE12021+GSE55457 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }else{
    plot(roc,
         print.auc=T,
         add=T,
         print.auc.x=0.5,print.auc.y=(0.5-i*0.05+0.05),
         print.auc.pattern=paste(colnames(exp2_test)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main="GSE12021+GSE55457 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }
  i<-i+1
}
dev.off()


pdf(file = 'test3-ROC.all.pdf',w=5,h=5)
for (i in c(1:len)) {
  roc<-roc(exp2_test$group,exp2_test[,i],levels=c("Case","Normal"))
  # roc = smooth(roc)
  if (i ==1) {
    plot(roc,
         print.auc=T,
         print.auc.x=0.5,print.auc.y=0.5,
         print.auc.pattern=paste(colnames(exp2_test)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main="GSE12021+GSE55457 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }else{
    plot(roc,
         print.auc=T,
         add=T,
         print.auc.x=0.5,print.auc.y=(0.5-i*0.05+0.05),
         print.auc.pattern=paste(colnames(exp2_test)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main="GSE12021+GSE55457 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }
  i<-i+1
}
dev.off()


##------test 2-------------
rm(list = ls())
setwd('E:/02-OA')
if(!dir.exists("06_roc")){
  dir.create("06_roc")}
setwd("06_roc")
library(tidyverse)
#data <- read.csv('../00_Rowdata/dat(GSE14905).csv',row.names = 1)
#group <- read.csv('../00_Rowdata/group(GSE14905).csv',row.names = 1)
test_data <- read.csv(file ="../00_rawdata/05GSE55457/02normalize/GSE55457.txt",
                      sep="\t",row.names = 1)
test_group <- read.csv(file ="../00_rawdata/05GSE55457/02normalize/GSE55457.group.txt",
                       sep="\t")
test_group$group <- factor(test_group$group, levels = c("Case", "Normal"))

##候选基因
#hub_gene <- read.csv('../06_venn/feature_genes.txt', header = F)
hub_gene <- read.csv('../05_model/13.svm_lasso.txt', header = F)
# hub_gene <- c("IL1R2", "IRAK3", "SLA")
# hub_gene <- c("NPTX2","LOXL1","COL1A2","C3AR1")
# hub_gene <- as.data.frame(hub_gene)
# colnames(hub_gene) <- c("V1")
#gene <- read.csv("../06_violin/DE-GSE165082-wilcox_result.csv",header = T,row.names = 1)
#data <- t(data)
data_candi_test <- test_data[rownames(test_data)%in%hub_gene$V1,test_group$sample]


exp_test <- data_candi_test

#library(pROC)
#library(ggplot2)
exp2_test <- t(exp_test)
exp2_test <- cbind(exp2_test,test_group)
# 降序排序
# exp2_test <- exp2_test[order(exp2_test$group,decreasing = T),]

## 绘制ROC曲线
library(pROC)
# corlor <- c("#FF2E63","#87CEFA")
corlor <- c("#FF2E63","#FF4500","firebrick","blue","#3CB371","gray50","#E7B800","#87CEFA","#CD5C5C","#F0E68C","#3CB371")
len <- length(colnames(exp2_test))-2
# png(file = 'test-ROC.all.png',width = 400,height = 400)
tiff(file="test2-ROC.tiff",width = 5, height = 5, units = "in",
     bg = "white", res = 500)
for (i in c(1:len)) {
  roc<-roc(exp2_test$group,exp2_test[,i],levels=c("Case","Normal"))
  # 报错Error in smooth_roc_binormal(roc, n) : ROC curve not smoothable (not enough points).则把roc = smooth(roc)注释#
  # roc = smooth(roc)
  if (i ==1) {
    plot(roc,
         print.auc=T,
         print.auc.x=0.5,print.auc.y=0.5,
         print.auc.pattern=paste(colnames(exp2_test)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main=" GSE55457 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }else{
    plot(roc,
         print.auc=T,
         add=T,
         print.auc.x=0.5,print.auc.y=(0.5-i*0.05+0.05),
         print.auc.pattern=paste(colnames(exp2_test)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main="GSE55457 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }
  i<-i+1
}
dev.off()


pdf(file = 'test2-ROC.all.pdf',w=5,h=5)
for (i in c(1:len)) {
  roc<-roc(exp2_test$group,exp2_test[,i],levels=c("Case","Normal"))
  # roc = smooth(roc)
  if (i ==1) {
    plot(roc,
         print.auc=T,
         print.auc.x=0.5,print.auc.y=0.5,
         print.auc.pattern=paste(colnames(exp2_test)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main="GSE55457 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }else{
    plot(roc,
         print.auc=T,
         add=T,
         print.auc.x=0.5,print.auc.y=(0.5-i*0.05+0.05),
         print.auc.pattern=paste(colnames(exp2_test)[i],"=",'%.3f'),
         #auc.polygon=T,
         #auc.polygon.con="#fff7f7",
         #grid=c(0.5,0.2),
         grid.col=c("black","black"),
         #print.thres=T,
         main="GSE55457 ROC curve",
         col=corlor[i],
         legacy.axes=T)
  }
  i<-i+1
}
dev.off()



#######表达量验证
##############
#-----------------------------训练---------------------------------------
# 是否画图
rm(list=ls())
plot <- T #T/F

hub_gene <- read.csv('../05_model/13.svm_lasso.txt', header = F)

library(tidyverse)
train_data <- read.csv(file ="../00_rawdata/02GSE51588/02normalize/GSE51588.txt",
                       sep="\t",row.names = 1)
train_group <- read.csv(file ="../00_rawdata/02GSE51588/02normalize/GSE51588.group.txt",
                        sep="\t")

hubgene <-  hub_gene$V1
expr <- train_data
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
library(rstatix)
stat_res <- dat2 %>% 
  group_by(Gene) %>% 
  wilcox_test(expression ~ Group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr,none
  add_significance("p")
stat_res  


# 画图数据准备
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



# 画图
library(ggsci)
library(ggplot2)
library(ggpubr)

if(plot == T){
  ggviolin(
    boxplot_dat,
    x = "Group",
    y = "value",
    fill = "Group", 
    add = "boxplot",
    palette = c( "#B72230","#104680"),
    add.params = list(fill = "white"),
    #trim = T
    # add = "jitter",
    # ylim = c(4, 8.8)
  )+stat_pvalue_manual(stat_res,
                       y.position = 3, # 有几个基因就设置几个值
                       size = 6,
                       color='red',
                       family = "Times",
                       label = "p.signif",
                       #parse = T,
                       face = "bold")+
    ggtitle("GSE51588")+
    theme_bw() + xlab("")+ylab("Expression")+
    theme(axis.title.x =element_text(size=20,family = "Times", face = "bold"),
          axis.text.x =element_text(size=18,family = "Times", face = "bold",color='black'),
          axis.title.y =element_text(size=20,family = "Times", face = "bold"),
          axis.text.y=element_text(size=18,family = "Times", face = "bold"),
          strip.text = element_text(size = 20,family = "Times", face = "bold"),
          plot.title=element_text(size=22,family = "Times", face = "bold",hjust=0.5),
          legend.position = "none")+
    facet_wrap(.~Gene,nrow=1
               ,scales = "free_y"
    )+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), 
          legend.title=element_text(size=15) , legend.text=element_text(size=14))
  ggsave(width=10,height=8,'train_boxplot.pdf')
  ggsave(width=10,height=8,'train_boxplot.png')
}


# 验证集 ---------------------------------------------------------------
# ----------------------12021----------------------------------------
rm(list = ls())
setwd('E:/02-OA')
if(!dir.exists("06_roc")){
  dir.create("06_roc")}
setwd("06_roc")

test_data <- read.csv(file ="../00_rawdata/04GSE12021/02normalize/GSE12021.txt",
                      row.names=1, sep="\t")
test_group <- read.csv(file ="../00_rawdata/04GSE12021/02normalize/GSE12021.group.txt",
                       sep="\t")
# hub_gene <- read.csv('../05_model/lasso_svm_rf_gene.txt', header = F)
# hub_gene <- c("IL1R2", "IRAK3", "SLA")
#hub_gene <- c("NPTX2","LOXL1","COL1A2","C3AR1")
#hub_gene <- as.data.frame(hub_gene)
#colnames(hub_gene) <- c("V1")
hub_gene <- read.csv('../05_model/13.svm_lasso.txt', header = F)
hubgene <- hub_gene$V1
expr <- test_data
test_group <- test_group


library(dplyr)
test_group <- test_group %>%dplyr::select(sample, group)
condition <- test_group
condition$group <- factor(condition$group,levels=c("Normal","Case"))
# condition <- condition[order(condition$group,decreasing = F),]

#hubgene<- c("FANCI","TRIP13","UBE2T","FEN1","CENPN")
hubgene_expr<-expr[hubgene ,]%>%t(.)%>%data.frame(.)
# hubgene_expr<-hubgene_expr[,order(colnames(hubgene_expr),decreasing = F)]

train_Group<-data.frame(sample=condition$sample,Group=condition$group) 

##按照字母进行排序
dat <-hubgene_expr %>% as.data.frame %>% tibble::rownames_to_column(var = "sample")
identical(dat$sample,train_Group$sample)
dat <- merge(dat, train_Group, by = "sample")
dat2 <- tidyr::gather(dat, Gene, expression, -c(sample, Group))
library(rstatix)
stat_res_test <- dat2 %>% 
  group_by(Gene) %>% 
  wilcox_test(expression ~ Group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr,none
  add_significance("p")
stat_res_test  

boxplot_dat_test <-data.frame(cbind(hubgene_expr,condition$group))
colnames(boxplot_dat_test)[ncol(boxplot_dat_test)]<-'Group'
head(boxplot_dat_test)
boxplot_dat_test$id <- rownames(boxplot_dat_test)
library(reshape2)
boxplot_dat_test <- melt(boxplot_dat_test,id = c("id","Group"))
boxplot_dat_test$value <-as.numeric(boxplot_dat_test$value)
boxplot_dat_test$Group<- factor(boxplot_dat_test$Group)
colnames(boxplot_dat_test)<-c('id','Group','Gene','value')

boxplot_dat_test<-boxplot_dat_test[order(boxplot_dat_test$Gene,decreasing = T),]
boxplot_dat_test$value <- log2(boxplot_dat_test$value+1)
plot <- T
# 画图
if(plot == T){
  
  # 
  # boxplot_dat1 <- boxplot_dat %>%  
  #   filter(Gene == "ISG15")
  # boxplot_dat2 <- boxplot_dat1 %>%  
  #   filter(Group == "Control")
  # median(boxplot_dat2$value)
  # 
  # boxplot_dat3 <- boxplot_dat1 %>%  
  #   filter(Group == "SLE")
  # median(boxplot_dat3$value)
  # 
  # 
  # 
  # 
  # boxplot_dat1 <- boxplot_dat$value[which(boxplot_dat$Gene == "ISG15")]
  # 
  # mean(boxplot_dat1)
  # 
  # Control_sample <- group$sample[which(group$group == group2)]
  # 
  
  boxplot_dat_test$Group <- factor(boxplot_dat_test$Group,levels = c("Case","Normal"))
  
  library(ggsci)
  library(ggplot2)
  library(ggpubr)
  ggviolin(
    boxplot_dat_test,
    x = "Group",
    y = "value",
    fill = "Group", 
    add = "boxplot",
    palette = c( "#B72230","#104680"),
    add.params = list(fill = "white"),
    #trim = T
    # add = "jitter",
    # ylim = c(4, 8.8)
  )+stat_pvalue_manual(stat_res_test,
                       y.position = c(12,15,12), # 有几个基因就设置几个值
                       size = 6,
                       color='red',
                       family = "Times",
                       label = "p.signif",
                       #parse = T,
                       face = "bold")+
    ggtitle("GSE12021")+
    theme_bw() + xlab("")+ylab("Expression")+
    theme(axis.title.x =element_text(size=20,family = "Times", face = "bold"),
          axis.text.x =element_text(size=18,family = "Times", face = "bold",color='black'),
          axis.title.y =element_text(size=20,family = "Times", face = "bold"),
          axis.text.y=element_text(size=18,family = "Times", face = "bold"),
          strip.text = element_text(size = 20,family = "Times", face = "bold"),
          plot.title=element_text(size=22,family = "Times", face = "bold",hjust=0.5),
          legend.position = "none")+
    facet_wrap(.~Gene,nrow=1
               ,scales = "free_y"
    )+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), legend.title=element_text(size=15) , legend.text=element_text(size=14))
  ggsave(width=10,height=8,'test_boxplot.pdf')
  ggsave(width=10,height=8,'test_boxplot.png')
  
}



# ----------------------55457----------------------------------------
rm(list = ls())
setwd('E:/02-OA')
if(!dir.exists("06_roc")){
  dir.create("06_roc")}
setwd("06_roc")

test_data <- read.csv(file ="../00_rawdata/05GSE55457/02normalize/GSE55457.txt",
                      row.names=1, sep="\t")
test_group <- read.csv(file ="../00_rawdata/05GSE55457/02normalize/GSE55457.group.txt",
                       sep="\t")
# hub_gene <- read.csv('../05_model/lasso_svm_rf_gene.txt', header = F)
# hub_gene <- c("IL1R2", "IRAK3", "SLA")
#hub_gene <- c("NPTX2","LOXL1","COL1A2","C3AR1")
#hub_gene <- as.data.frame(hub_gene)
#colnames(hub_gene) <- c("V1")
#hub_gene <- read.csv('../04_venn/feature_genes.txt', header = F)
hub_gene <- read.csv('../05_model/13.svm_lasso.txt', header = F)
hubgene <- hub_gene$V1
expr <- test_data
test_group <- test_group


library(dplyr)
test_group <- test_group %>%dplyr::select(sample, group)
condition <- test_group
condition$group <- factor(condition$group,levels=c("Normal","Case"))
# condition <- condition[order(condition$group,decreasing = F),]

#hubgene<- c("FANCI","TRIP13","UBE2T","FEN1","CENPN")
hubgene_expr<-expr[hubgene ,]%>%t(.)%>%data.frame(.)
# hubgene_expr<-hubgene_expr[,order(colnames(hubgene_expr),decreasing = F)]

train_Group<-data.frame(sample=condition$sample,Group=condition$group) 

##按照字母进行排序
dat <-hubgene_expr %>% as.data.frame %>% tibble::rownames_to_column(var = "sample")
identical(dat$sample,train_Group$sample)
dat <- merge(dat, train_Group, by = "sample")
dat2 <- tidyr::gather(dat, Gene, expression, -c(sample, Group))
library(rstatix)
stat_res_test <- dat2 %>% 
  group_by(Gene) %>% 
  wilcox_test(expression ~ Group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr,none
  add_significance("p")
stat_res_test  

boxplot_dat_test <-data.frame(cbind(hubgene_expr,condition$group))
colnames(boxplot_dat_test)[ncol(boxplot_dat_test)]<-'Group'
head(boxplot_dat_test)
boxplot_dat_test$id <- rownames(boxplot_dat_test)
library(reshape2)
boxplot_dat_test <- melt(boxplot_dat_test,id = c("id","Group"))
boxplot_dat_test$value <-as.numeric(boxplot_dat_test$value)
boxplot_dat_test$Group<- factor(boxplot_dat_test$Group)
colnames(boxplot_dat_test)<-c('id','Group','Gene','value')

boxplot_dat_test<-boxplot_dat_test[order(boxplot_dat_test$Gene,decreasing = T),]
boxplot_dat_test$value <- log2(boxplot_dat_test$value+1)
plot <- T
# 画图
if(plot == T){
  
  # 
  # boxplot_dat1 <- boxplot_dat %>%  
  #   filter(Gene == "ISG15")
  # boxplot_dat2 <- boxplot_dat1 %>%  
  #   filter(Group == "Control")
  # median(boxplot_dat2$value)
  # 
  # boxplot_dat3 <- boxplot_dat1 %>%  
  #   filter(Group == "SLE")
  # median(boxplot_dat3$value)
  # 
  # 
  # 
  # 
  # boxplot_dat1 <- boxplot_dat$value[which(boxplot_dat$Gene == "ISG15")]
  # 
  # mean(boxplot_dat1)
  # 
  # Control_sample <- group$sample[which(group$group == group2)]
  # 
  
  boxplot_dat_test$Group <- factor(boxplot_dat_test$Group,levels = c("Case","Normal"))
  
  library(ggsci)
  library(ggplot2)
  library(ggpubr)
  ggviolin(
    boxplot_dat_test,
    x = "Group",
    y = "value",
    fill = "Group", 
    add = "boxplot",
    palette = c( "#B72230","#104680"),
    add.params = list(fill = "white"),
    #trim = T
    # add = "jitter",
    # ylim = c(4, 8.8)
  )+stat_pvalue_manual(stat_res_test,
                       y.position = c(12, 15, 12), # 有几个基因就设置几个值
                       size = 6,
                       color='red',
                       family = "Times",
                       label = "p.signif",
                       #parse = T,
                       face = "bold")+
    ggtitle("Test")+
    theme_bw() + xlab("")+ylab("Expression")+
    theme(axis.title.x =element_text(size=20,family = "Times", face = "bold"),
          axis.text.x =element_text(size=18,family = "Times", face = "bold",color='black'),
          axis.title.y =element_text(size=20,family = "Times", face = "bold"),
          axis.text.y=element_text(size=18,family = "Times", face = "bold"),
          strip.text = element_text(size = 20,family = "Times", face = "bold"),
          plot.title=element_text(size=22,family = "Times", face = "bold",hjust=0.5),
          legend.position = "none")+
    facet_wrap(.~Gene,nrow=1
               ,scales = "free_y"
    )+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), legend.title=element_text(size=15) , legend.text=element_text(size=14))
  ggsave(width=10,height=8,'test2_boxplot.pdf')
  ggsave(width=10,height=8,'test2_boxplot.png')
  
}



# ----------------------06val_combat----------------------------------------
rm(list = ls())
setwd('E:/02-OA')
if(!dir.exists("06_roc")){
  dir.create("06_roc")}
setwd("06_roc")

test_data <- read.csv(file ="../00_rawdata/06val_combat/val.combat.exp.txt",
                      row.names=1, sep="\t")
test_group <- read.csv(file ="../00_rawdata/06val_combat/val.group.txt",
                       sep="\t")
# hub_gene <- read.csv('../05_model/lasso_svm_rf_gene.txt', header = F)
# hub_gene <- c("IL1R2", "IRAK3", "SLA")
#hub_gene <- c("NPTX2","LOXL1","COL1A2","C3AR1")
#hub_gene <- as.data.frame(hub_gene)
#colnames(hub_gene) <- c("V1")
#hub_gene <- read.csv('../04_venn/feature_genes.txt', header = F)
hub_gene <- read.csv('../05_model/13.svm_lasso.txt', header = F)
hubgene <- hub_gene$V1
expr <- test_data
test_group <- test_group


library(dplyr)
test_group <- test_group %>%dplyr::select(sample, group)
condition <- test_group
condition$group <- factor(condition$group,levels=c("Normal","Case"))
# condition <- condition[order(condition$group,decreasing = F),]

#hubgene<- c("FANCI","TRIP13","UBE2T","FEN1","CENPN")
hubgene_expr<-expr[hubgene ,]%>%t(.)%>%data.frame(.)
# hubgene_expr<-hubgene_expr[,order(colnames(hubgene_expr),decreasing = F)]

train_Group<-data.frame(sample=condition$sample,Group=condition$group) 

##按照字母进行排序
dat <-hubgene_expr %>% as.data.frame %>% tibble::rownames_to_column(var = "sample")
identical(dat$sample,train_Group$sample)
dat <- merge(dat, train_Group, by = "sample")
dat2 <- tidyr::gather(dat, Gene, expression, -c(sample, Group))
library(rstatix)
stat_res_test <- dat2 %>% 
  group_by(Gene) %>% 
  wilcox_test(expression ~ Group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr,none
  add_significance("p")
stat_res_test  

boxplot_dat_test <-data.frame(cbind(hubgene_expr,condition$group))
colnames(boxplot_dat_test)[ncol(boxplot_dat_test)]<-'Group'
head(boxplot_dat_test)
boxplot_dat_test$id <- rownames(boxplot_dat_test)
library(reshape2)
boxplot_dat_test <- melt(boxplot_dat_test,id = c("id","Group"))
boxplot_dat_test$value <-as.numeric(boxplot_dat_test$value)
boxplot_dat_test$Group<- factor(boxplot_dat_test$Group)
colnames(boxplot_dat_test)<-c('id','Group','Gene','value')

boxplot_dat_test<-boxplot_dat_test[order(boxplot_dat_test$Gene,decreasing = T),]
boxplot_dat_test$value <- log2(boxplot_dat_test$value+1)
plot <- T
# 画图
if(plot == T){
  
  # 
  # boxplot_dat1 <- boxplot_dat %>%  
  #   filter(Gene == "ISG15")
  # boxplot_dat2 <- boxplot_dat1 %>%  
  #   filter(Group == "Control")
  # median(boxplot_dat2$value)
  # 
  # boxplot_dat3 <- boxplot_dat1 %>%  
  #   filter(Group == "SLE")
  # median(boxplot_dat3$value)
  # 
  # 
  # 
  # 
  # boxplot_dat1 <- boxplot_dat$value[which(boxplot_dat$Gene == "ISG15")]
  # 
  # mean(boxplot_dat1)
  # 
  # Control_sample <- group$sample[which(group$group == group2)]
  # 
  
  boxplot_dat_test$Group <- factor(boxplot_dat_test$Group,levels = c("Case","Normal"))
  
  library(ggsci)
  library(ggplot2)
  library(ggpubr)
  ggviolin(
    boxplot_dat_test,
    x = "Group",
    y = "value",
    fill = "Group", 
    add = "boxplot",
    palette = c( "#B72230","#104680"),
    add.params = list(fill = "white"),
    #trim = T
    # add = "jitter",
    # ylim = c(4, 8.8)
  )+stat_pvalue_manual(stat_res_test,
                       y.position = c(3.8,4.5,3.8), # 有几个基因就设置几个值
                       size = 6,
                       color='red',
                       family = "Times",
                       label = "p.signif",
                       #parse = T,
                       face = "bold")+
    ggtitle("GSE12021+GSE55457")+
    theme_bw() + xlab("")+ylab("Expression")+
    theme(axis.title.x =element_text(size=20,family = "Times", face = "bold"),
          axis.text.x =element_text(size=18,family = "Times", face = "bold",color='black'),
          axis.title.y =element_text(size=20,family = "Times", face = "bold"),
          axis.text.y=element_text(size=18,family = "Times", face = "bold"),
          strip.text = element_text(size = 20,family = "Times", face = "bold"),
          plot.title=element_text(size=22,family = "Times", face = "bold",hjust=0.5),
          legend.position = "none")+
    facet_wrap(.~Gene,nrow=1
               ,scales = "free_y"
    )+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), legend.title=element_text(size=15) , legend.text=element_text(size=14))
  ggsave(width=10,height=8,'test3_boxplot.pdf')
  ggsave(width=10,height=8,'test3_boxplot.png')
  
}

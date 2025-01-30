rm(list=ls())
setwd("E:/LZ/24105")  #modify work path
if(! dir.exists("08_roc")){
  dir.create("08_roc")
}
setwd("./08_roc")
library(survival)
library(timeROC)
library(ggsci)
rt <- read.table("../00_rawdata/GSE14520/risk.txt", header = T, sep = "\t", check.names = F, row.names = 1)

ROC_rt <- timeROC(T = rt$futimes, delta = rt$fustate, marker = rt$riskScore,
                  cause = 1, weighting = 'aalen', times = c(1, 3, 5),
                  ROC = TRUE, iid = FALSE)

pdf(file = "ROC.pdf", width = 5, height = 5)
par(oma = c(0.5,1,0,1), font.lab = 1.5, font.axis = 1.5)
plot(ROC_rt, time = 1, title = FALSE, lwd = 2, col = "#DC0000FF")
plot(ROC_rt, time = 3, title = FALSE, lwd = 2, col = "#00A087FF", add = TRUE)
plot(ROC_rt, time = 5, title = FALSE, lwd = 2, col = "#4DBBD5FF", add = TRUE)
#plot(ROC_rt, time = 10, title = FALSE, lwd = 2, col = "#F39B7FFF", add = TRUE)

#plot(ROC_rt,time=8,title=FALSE,lwd=2,col='blue',add=TRUE)
#title( main = "ROC curve")
legend('bottomright',
        c(
        paste0('AUC at 1 years: ', round(ROC_rt$AUC[1], 3)),
        paste0('AUC at 3 years: ', round(ROC_rt$AUC[2], 3)),
        paste0('AUC at 5 years: ', round(ROC_rt$AUC[3], 3))),
        col = c("#DC0000FF", "#00A087FF", "#4DBBD5FF"),  #, "#F39B7FFF"
       lwd = 2, bty = 'n')
dev.off()
  
tiff(file="ROC.tiff",width = 5, height = 5, units = "in",
     bg = "white", res = 500)
par(oma = c(0.5,1,0,1), font.lab = 1.5, font.axis = 1.5)
plot(ROC_rt, time = 1, title = FALSE, lwd = 2, col = "#DC0000FF")
plot(ROC_rt, time = 3, title = FALSE, lwd = 2, col = "#00A087FF", add = TRUE)
plot(ROC_rt, time = 5, title = FALSE, lwd = 2, col = "#4DBBD5FF", add = TRUE)
#plot(ROC_rt, time = 10, title = FALSE, lwd = 2, col = "#F39B7FFF", add = TRUE)

#plot(ROC_rt,time=8,title=FALSE,lwd=2,col='blue',add=TRUE)
#title( main = "ROC curve")
legend('bottomright',
       c(
         paste0('AUC at 1 years: ', round(ROC_rt$AUC[1], 3)),
         paste0('AUC at 3 years: ', round(ROC_rt$AUC[2], 3)),
         paste0('AUC at 5 years: ', round(ROC_rt$AUC[3], 3))),
       col = c("#DC0000FF", "#00A087FF", "#4DBBD5FF"),
       lwd = 2, bty = 'n')
dev.off()


###############################################################################
###############           TCGA            #####################################
###############################################################################
rm(list=ls())
setwd("E:/LZ/24105")  #modify work path
if(! dir.exists("08_roc")){
  dir.create("08_roc")
}
setwd("./08_roc")

plot <- T #T/F

hub_gene <- read.csv('../07_multicox/multiCox.csv')
# hub_gene <- read.csv('../06_lasso/geneCoef.txt',sep="\t",header = T)

library(tidyverse)
train_data <- read.csv("../00_rawdata/normalizeExp.txt", 
                       sep = "\t", row.names=1,header = T)
#colnames(train_data) <- substring(colnames(train_data), 1, 16)
colnames(train_data) <- gsub("\\.","-", colnames(train_data))

train_group <- read.csv(file = "../00_rawdata/group.txt",
                        sep="\t",header = T)
train_group$group <- ifelse(train_group$group == "tumor","Case","Normal")

hubgene <-  hub_gene[,1]
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

boxplot_dat$value[boxplot_dat$Gene %in% c("HOXC8", "HS3ST5")] <- 
  log2(boxplot_dat$value[boxplot_dat$Gene %in% c("HOXC8", "HS3ST5")])

if(plot == T){
  ggviolin(
    boxplot_dat,
    x = "Group",
    y = "value",
    fill = "Group", 
    add = "boxplot",
    palette = c("#ef5675","#7a9ac0"),
    add.params = list(fill = "white"),
    #trim = T
    # add = "jitter",
    # ylim = c(4, 8.8)
  )+stat_pvalue_manual(stat_res,
                       y.position = c(19,14), # 有几个基因就设置几个值
                       size = 6,
                       color='red',
                       family = "Times",
                       label = "p.signif",
                       #parse = T,
                       face = "bold")+
    ggtitle("TCGA")+
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
  ggsave(width=12,height=6,'TCGA.pdf')
  ggsave(width=12,height=6,'TCGA.png')
}


###############################################################################
###############            104580            #####################################
###############################################################################
rm(list=ls())

plot <- T #T/F

hub_gene <- read.csv('../07_multicox/multiCox.csv')
# hub_gene <- read.csv('../06_lasso/geneCoef.txt',sep="\t",header = T)
colnames(hub_gene)[1] <- "V1"

library(tidyverse)
train_data <- read.csv("../00_rawdata/GSE104580/02normalize/GSE104580.txt", 
                       sep = "\t", row.names=1,header = T)
#colnames(train_data) <- substring(colnames(train_data), 1, 16)
#colnames(train_data) <- gsub("\\.","-", colnames(train_data))

train_group <- read.csv(file = "../00_rawdata/GSE104580/02normalize/GSE104580.group.txt",
                        sep="\t",header = T)

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

boxplot_dat$value[boxplot_dat$Gene %in% c("HOXC8", "HS3ST5")] <- 
  log2(boxplot_dat$value[boxplot_dat$Gene %in% c("HOXC8", "HS3ST5")])

boxplot_dat$value <- log2(boxplot_dat$value+1)
#boxplot_dat$Group <- ifelse(boxplot_dat$Group == "Case", "Response", "Not-response")
#boxplot_dat$Group <- factor(boxplot_dat$Group,levels = c("Response", "Not-response"))
if(plot == T){
  ggviolin(
    boxplot_dat,
    x = "Group",
    y = "value",
    fill = "Group", 
    add = "boxplot",
    palette = c("#ef5675","#7a9ac0"),
    add.params = list(fill = "white"),
    #trim = T
    # add = "jitter",
    # ylim = c(4, 8.8)
  )+stat_pvalue_manual(stat_res,
                       y.position = 5, # 有几个基因就设置几个值
                       size = 6,
                       color='red',
                       family = "Times",
                       label = "p.signif",
                       #parse = T,
                       face = "bold")+
    ggtitle("GSE104580")+
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
  ggsave(width=12,height=6,'GSE104580.pdf')
  ggsave(width=12,height=6,'GSE104580.png')
}



###############################################################################
###############            14520            #####################################
###############################################################################
rm(list=ls())

plot <- T #T/F

hub_gene <- read.csv('../07_multicox/multiCox.csv')
# hub_gene <- read.csv('../06_lasso/geneCoef.txt',sep="\t",header = T)
colnames(hub_gene)[1] <- "V1"

library(tidyverse)
train_data <- read.csv("../00_rawdata/GSE14520/GSE14520.txt", 
                       sep = "\t", row.names=1,header = T)
#colnames(train_data) <- substring(colnames(train_data), 1, 16)
#colnames(train_data) <- gsub("\\.","-", colnames(train_data))

train_group <- read.csv(file = "../00_rawdata/GSE14520/GSE14520.group2.txt",
                        sep="\t",header = T)

hubgene <-  hub_gene$V1
expr <- train_data
train_group <- train_group
library(dplyr)
train_group <- train_group %>%dplyr::select(sample, group)
condition <- train_group

# hubgene<- c("GAS6","CD14","RAB17","CYP27A1","S100A9","AQP9","IGFBP2")
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

boxplot_dat$value[boxplot_dat$Gene %in% c("HOXC8", "HS3ST5")] <- 
  log2(boxplot_dat$value[boxplot_dat$Gene %in% c("HOXC8", "HS3ST5")])

if(plot == T){
  ggviolin(
    boxplot_dat,
    x = "Group",
    y = "value",
    fill = "Group", 
    add = "boxplot",
    palette = c("#ef5675","#7a9ac0"),
    add.params = list(fill = "white"),
    #trim = T
    # add = "jitter",
    # ylim = c(4, 8.8)
  )+stat_pvalue_manual(stat_res,
                       y.position = c(14,14), # 有几个基因就设置几个值
                       size = 6,
                       color='red',
                       family = "Times",
                       label = "p.signif",
                       #parse = T,
                       face = "bold")+
    ggtitle("GSE14520")+
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
  ggsave(width=12,height=6,'GSE14520.pdf')
  ggsave(width=12,height=6,'GSE14520.png')
}




###################################################################################
#################################GSE14520 SUBGROUP#############################
###################################################################################
# 是否画图
rm(list=ls())
plot <- T #T/F

hub_gene <- read.csv('../07_multicox/multiCox.csv')

library(tidyverse)
train_data <- read.csv(file ="../00_rawdata/GSE14520/GSE14520.txt",
                       sep="\t",row.names = 1)
# group: TACE与否 
train_group <- read.csv(file = "../00_rawdata/GSE14520/GSE14520.group.txt",
                   sep="\t",header = T)
colnames(train_group) <- c('sample', 'group')
# train_group的group上paste clust
# group1$group <- paste0("clust",group1$group)
train_group$group <- ifelse(train_group$group == "Case","TACE","Resection Only")
train_group$group <- factor(train_group$group, levels = c("TACE","Resection Only"))



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
    palette = c("#1f77b4", "#ff7f0e", "#2ca02c"),
    add.params = list(fill = "white"),
    #trim = T
    # add = "jitter",
    # ylim = c(4, 8.8)
  )+stat_pvalue_manual(stat_res,
                       y.position = c(2,2.5,2.25), # 有几个基因就设置几个值
                       size = 6,
                       color='red',
                       family = "Times",
                       label = "p.signif",
                       #parse = T,
                       face = "bold")+
    ggtitle(" ")+
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
  ggsave(width=10,height=8,'control_clust_boxplot.pdf')
  ggsave(width=10,height=8,'control_clust_boxplot.png')
}







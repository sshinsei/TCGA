rm(list=ls())
setwd("E:/LZ/24080")  #modify work path
if(! dir.exists("04_clustsurvival")){
  dir.create("04_clustsurvival")
}
setwd("./04_clustsurvival")  
library(survival)
library(ggsci)
# 合并亚型文件和生存数据
os <- read.table(file = "../00_rawdata/clinical.csv", sep = ",", header = T)
clust <- read.csv("../03_clust/cluster.txt", sep="\t",header = F)
colnames(clust) <- c('sample', 'group')
clust$group <- paste0("clust",clust$group)
rt <- left_join(clust, os, by="sample")
rt = na.omit(rt)
rt = rt[rt$futimes > 30, ] # 344

# rt = read.table("../00_rawdata/risk.txt", header = T, sep = "\t", row.names = 1)
diff = survdiff(Surv(futimes, fustate) ~ group, data = rt)
pValue = 1 - pchisq(diff$chisq, df = 1)
#pValue=round(pValue,3)
pValue = signif(pValue, 4)
pValue = format(pValue, scientific = TRUE)

fit <- survfit(Surv(futimes, fustate) ~ group, data = rt)
summary(fit)  #See five-year survival rates
library(ggsci)
library(survival)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(survminer)
library(TSHRC)
p <- ggsurvplot(fit,  #创建的拟合对象
                data = rt,  #指定变量数据来源
                conf.int = TRUE,  #显示置信区间
                pval = TRUE,  #添加P值
                risk.table = TRUE,  #绘制累计风险曲线
                surv.median.line = "hv",  #添加中位生存时间线
                add.all = FALSE,  #添加总患者生存曲线
                palette = c("#DC0000FF", "#4DBBD5FF"))  #自定义调色板
p
tiff(file="survival.tiff",width = 6, height = 6.5, units = "in",
     bg = "white", res = 500)
ggsurvplot(fit,  #创建的拟合对象
           data = rt,  #指定变量数据来源
           conf.int = TRUE,  #显示置信区间
           pval = TRUE,  #添加P值
           risk.table = TRUE,  #绘制累计风险曲线
           surv.median.line = "hv",  #添加中位生存时间线
           add.all = FALSE,  #添加总患者生存曲线
           palette = c("#4DBBD5FF", "#DC0000FF"))  #自定义调色板
dev.off()

pdf(file = "survival.pdf", width = 6, height = 6.5)
ggsurvplot(fit,  #创建的拟合对象
           data = rt,  #指定变量数据来源
           conf.int = TRUE,  #显示置信区间
           pval = TRUE,  #添加P值
           risk.table = TRUE,  #绘制累计风险曲线
           surv.median.line = "hv",  #添加中位生存时间线
           add.all = FALSE,  #添加总患者生存曲线
           palette = c("#4DBBD5FF", "#DC0000FF"))  #自定义调色板
dev.off()


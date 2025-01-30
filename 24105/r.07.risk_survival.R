rm(list=ls())
setwd("E:/LZ/24105")  #modify work path
if(! dir.exists("07_survival")){
  dir.create("07_survival")
}
setwd("./07_survival")  
library(survival)
library(ggsci)
rt = read.table("../00_rawdata/GSE14520/risk.txt", header = T, sep = "\t", row.names = 1)
diff = survdiff(Surv(futimes, fustate) ~ risk, data = rt)
pValue = 1 - pchisq(diff$chisq, df = 1)
#pValue=round(pValue,3)
pValue = signif(pValue, 4)
pValue = format(pValue, scientific = TRUE)

fit <- survfit(Surv(futimes, fustate) ~ risk, data = rt)
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
           palette = c("#DC0000FF", "#4DBBD5FF"))  #自定义调色板
dev.off()

pdf(file = "survival.pdf", width = 6, height = 6.5)
ggsurvplot(fit,  #创建的拟合对象
           data = rt,  #指定变量数据来源
           conf.int = TRUE,  #显示置信区间
           pval = TRUE,  #添加P值
           risk.table = TRUE,  #绘制累计风险曲线
           surv.median.line = "hv",  #添加中位生存时间线
           add.all = FALSE,  #添加总患者生存曲线
           palette = c("#DC0000FF", "#4DBBD5FF"))  #自定义调色板
dev.off()

################################################################################
############################     Adjuvant    ###################################
################################################################################

sub <- read.csv("../00_rawdata/GSE14520/GSE14520.subgroup.txt",header = T,sep = "\t")
adj <- rt[rownames(rt) %in% sub[sub$group == "Adjuvant",]$sample,]

diff = survdiff(Surv(futimes, fustate) ~ risk, data = adj)
pValue = 1 - pchisq(diff$chisq, df = 1)
#pValue=round(pValue,3)
pValue = signif(pValue, 4)
pValue = format(pValue, scientific = TRUE)

fit <- survfit(Surv(futimes, fustate) ~ risk, data = adj)
summary(fit)  #See five-year survival rates

p <- ggsurvplot(fit,  #创建的拟合对象
                data = rt,  #指定变量数据来源
                conf.int = TRUE,  #显示置信区间
                pval = TRUE,  #添加P值
                risk.table = TRUE,  #绘制累计风险曲线
                surv.median.line = "hv",  #添加中位生存时间线
                add.all = FALSE,  #添加总患者生存曲线
                palette = c("#DC0000FF", "#4DBBD5FF"))  #自定义调色板
p
tiff(file="adj.tiff",width = 6, height = 6.5, units = "in",
     bg = "white", res = 500)
ggsurvplot(fit,  #创建的拟合对象
           data = adj,  #指定变量数据来源
           conf.int = TRUE,  #显示置信区间
           pval = TRUE,  #添加P值
           risk.table = TRUE,  #绘制累计风险曲线
           surv.median.line = "hv",  #添加中位生存时间线
           add.all = FALSE,  #添加总患者生存曲线
           palette = c("#DC0000FF", "#4DBBD5FF"))  #自定义调色板
dev.off()

pdf(file = "adj.pdf", width = 6, height = 6.5)
ggsurvplot(fit,  #创建的拟合对象
           data = adj,  #指定变量数据来源
           conf.int = TRUE,  #显示置信区间
           pval = TRUE,  #添加P值
           risk.table = TRUE,  #绘制累计风险曲线
           surv.median.line = "hv",  #添加中位生存时间线
           add.all = FALSE,  #添加总患者生存曲线
           palette = c("#DC0000FF", "#4DBBD5FF"))  #自定义调色板
dev.off()



################################################################################
############################      Post       ###################################
################################################################################

sub <- read.csv("../00_rawdata/GSE14520/GSE14520.subgroup.txt",header = T,sep = "\t")
post <- rt[rownames(rt) %in% sub[sub$group == "Post",]$sample,]

diff = survdiff(Surv(futimes, fustate) ~ risk, data = post)
pValue = 1 - pchisq(diff$chisq, df = 1)
#pValue=round(pValue,3)
pValue = signif(pValue, 4)
pValue = format(pValue, scientific = TRUE)

fit <- survfit(Surv(futimes, fustate) ~ risk, data = post)
summary(fit)  #See five-year survival rates

p <- ggsurvplot(fit,  #创建的拟合对象
                data = rt,  #指定变量数据来源
                conf.int = TRUE,  #显示置信区间
                pval = TRUE,  #添加P值
                risk.table = TRUE,  #绘制累计风险曲线
                surv.median.line = "hv",  #添加中位生存时间线
                add.all = FALSE,  #添加总患者生存曲线
                palette = c("#DC0000FF", "#4DBBD5FF"))  #自定义调色板
p
tiff(file="post.tiff",width = 6, height = 6.5, units = "in",
     bg = "white", res = 500)
ggsurvplot(fit,  #创建的拟合对象
           data = post,  #指定变量数据来源
           conf.int = TRUE,  #显示置信区间
           pval = TRUE,  #添加P值
           risk.table = TRUE,  #绘制累计风险曲线
           surv.median.line = "hv",  #添加中位生存时间线
           add.all = FALSE,  #添加总患者生存曲线
           palette = c("#DC0000FF", "#4DBBD5FF"))  #自定义调色板
dev.off()

pdf(file = "post.pdf", width = 6, height = 6.5)
ggsurvplot(fit,  #创建的拟合对象
           data = post,  #指定变量数据来源
           conf.int = TRUE,  #显示置信区间
           pval = TRUE,  #添加P值
           risk.table = TRUE,  #绘制累计风险曲线
           surv.median.line = "hv",  #添加中位生存时间线
           add.all = FALSE,  #添加总患者生存曲线
           palette = c("#DC0000FF", "#4DBBD5FF"))  #自定义调色板
dev.off()


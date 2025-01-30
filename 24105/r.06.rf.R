rm(list=ls())
library(glmnet)
library(survival)
library(tidyverse)
setwd("E:/LZ/24105")  #modify work path
if(! dir.exists("06_rf")){
  dir.create("06_rf")
}
setwd("./06_rf")
set.seed(12348)

if (!require(randomForestSRC)) install.packages("randomForestSRC")
if (!require(survival)) install.packages("survival")

library(randomForestSRC)
library(survival)

rt = read.table("../00_rawdata/tumor_clinical.txt",
                header = T, sep = "\t", row.names = 1,
                check.names = FALSE) # 246
rt = na.omit(rt) # 241
rt = rt[rt$futimes > 30, ] # 241
#rt[ , 3:ncol(rt)] = log2(rt[ , 3:ncol(rt)] + 0.001)
rt[ , "futimes"] = rt[ , "futimes"]/365 # 344
#rt = na.omit(rt)
# p<0.05的特征基因
gene = read.table("../05_unicox/univariateCox_P0.001.csv", header = T, sep=",") # 93
#cols = match(gene$V1, colnames(rt))
rt = rt[ , c("futimes", "fustate", as.vector(gene[ , 1]))]

v.obj <- rfsrc(Surv(futimes, fustate) ~ ., data = rt, ntree = 100, nsplit = 10,
               na.action = "na.impute", tree.err = TRUE, importance = TRUE, block.size = 1)


# 获取特征重要性数值
importance_values <- data.frame(gene = names(v.obj$importance),importance = v.obj$importance)
importance_values <- importance_values %>% arrange(desc(importance))

# 设定阈值筛选重要基因
# threshold <- 0
# important_genes <- importance_values[importance_values$importance > threshold,]
important_genes <- importance_values[1:10,]
print(important_genes)
write.csv(important_genes,"rf_genes.csv")


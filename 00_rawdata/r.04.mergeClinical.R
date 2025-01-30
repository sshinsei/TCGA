rm(list=ls())
setwd("E:/LZ/24080")  #modify work path
if(! dir.exists("04_unicox")){
  dir.create("04_unicox")
}
setwd("./04_unicox")
library(tidyverse)
##----------要删除正常样本----------------##
norExp = read.table(file = "../02_DEG/normalizeExp.txt", 
                    sep = "\t", header = T, check.names = F) # 430,有一列是id
cli_data = read.table(file = "../00_rawdata/clinical.csv", sep = ",", header = T)
fea_genes = read.table(file = "../03_venn/feature_genes.txt", 
                       sep = "\t", header = F, check.names = F)
rows = na.omit(match(fea_genes$V1, norExp$id))
norExp_1 = norExp[rows, ]
unique(norExp_1$id == fea_genes$V1)  #check
rownames(norExp_1) = norExp_1$id
norExp_1 = norExp_1[ , -1]
norExp_1 = t(norExp_1)
norExp_1 = as.data.frame(norExp_1)
norExp_1$id = rownames(norExp_1)
#rows = na.omit(match(rownames(norExp_1), cli_data$Id))
#norExp_2 = cbind(norExp_1, cli_data[rows, ])
#norExp_2 = norExp_2[ , c((ncol(norExp_2)-2):ncol(norExp_2), 1:(ncol(norExp_2)-3))]
#unique(rownames(norExp_2) == norExp_2$Id)  #check

norExp_1$id <- substring(norExp_1$id, 1, 16) # 424


norExp_2 = merge(cli_data, norExp_1, by.x = "sample", by.y = "id") # 418
# 删除正常样本
norExp_2 = norExp_2[!grepl("11A",norExp_2$sample),] # 367
# 去重
norExp_2 <- norExp_2 %>% 
  distinct(sample, .keep_all = TRUE) # 367



write.table(norExp_2, file = "../00_rawdata/tumor_clinical.txt", sep = "\t", row.names = F, quote = F)

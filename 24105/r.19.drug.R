rm(list=ls())
setwd("E:/LZ/24105")
if (!dir.exists("19_drug")){dir.create("19_drug")}
setwd("19_drug")
library(pRRophetic)
library(edgeR)
library(tidyverse)
norExp = read.csv(file = "../02_DEG/normalizeExp.txt", sep="\t", 
                  header=T, check.names=F, row.names = 1)

colnames(norExp) <- substring(colnames(norExp),1,16)

# 删除正常样本
mask <- !grepl("11A", colnames(norExp))
norExp <- norExp[,mask]
write.table(norExp,"geneExp.txt",sep="\t",quote = F,
            row.names = T, col.names = T)
#norExp = norExp[ , -c(2:52)]  #delete normal sample
norExp = as.matrix(norExp)
# rownames(norExp) = norExp[ , 1]
# exp = norExp[ , 2:ncol(norExp)]
exp = norExp

dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp),
              dimnames = dimnames)
data <- as.matrix(data)
#data = avereps(data)
#data = data[rowMeans(data) > 1, ]
#data = log2(data + 1)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
#data(cgp2016ExprRma)
possibleDrugs2016 <- unique(drugData2016$Drug.name)
#possibleDrugs2016
head(possibleDrugs2016)

drug_sen = data.frame()
for (i in possibleDrugs2016[1:251]) {
  temp = try(predictedPtype <- pRRopheticPredict(testMatrix = data, 
                                                 drug = i, tissueType = "all", 
                                                 batchCorrect = "eb", selection = 1,
                                                 dataset = "cgp2016"), silent = TRUE)
  print(temp)
  if("try-error" %in% class(temp)){
    next
  } else {
    predictedPtype = t(as.data.frame(predictedPtype))
    rownames(predictedPtype) = i
    drug_sen = rbind(drug_sen, predictedPtype) 
  }
}
write.table(drug_sen, file = "drug_predict.txt", sep = "\t", quote = F, col.names = T)


save.image(file = "./IC50.RData")

## --------------------------------运行一次后直接load-----------------------------------------------------
rm(list=ls())
setwd("E:/LZ/24105")
if (!dir.exists("19_drug")){dir.create("19_drug")}
setwd("19_drug")
load("E:/LZ/24080/16_drug/IC50.RData")
#plot
library(dplyr)
library(stringr)
#drug_sen = read.table("drug_predict.txt", header = T, sep = "\t",
#check.names = F, row.names = 1) %>% t()
#rownames(drug_sen) = substring(rownames(drug_sen), 1, 12)
drug_sen = t(drug_sen)
drug_sen = as.data.frame(drug_sen)
drug_sen$id = rownames(drug_sen)
drug_sen$id = str_replace_all(drug_sen$id, "\\.", "-")
#colnames(drug_sen) = drug_sen[1, ]
#drug_sen = drug_sen[-1, ]
library(tidyverse)
#drug_sen$id = rownames(drug_sen)
#drug_sen$id = gsub(drug_sen$id, ".", "-")
risk = read.table("../00_rawdata/risk.txt", header = T, sep = "\t", check.names = F)
#risk = risk[order(risk$risk, decreasing = T), ]
#inte = intersect(risk$id, rownames(drug_sen))
#rownames(risk) = risk$id
#risk = risk[inte, ]
#rows = na.omit(match(rownames(risk), rownames(drug_sen)))
#rt1 = cbind(risk[, c(1,5)], drug_sen[rows, ])
rt1 = merge(risk[ , c(1, 5)], drug_sen, by.x = "id", by.y = "id")
rt1 = rt1[order(rt1$risk), ]
rt1$risk = factor(rt1$risk, levels = c("high", "low"))
rt1 = rt1[ , -1]
library(ggplot2)
library(ggsci)
library(ggpubr)
for (i in 2:ncol(rt1)) {
  da = na.omit(rt1[ , c(1, i)])
  #colnames(da)[2] = gsub("-", "_", colnames(da)[2])
  colname = colnames(da)[2]
  colnames(da)[2] = "expression"
  # 计算表达量
  expValue <- da %>% 
    group_by(risk) %>% 
    summarise(expValue = median(expression))
  # 提取表达量
  highval <- expValue$expValue[1]
  lowval <- expValue$expValue[2]
  test = wilcox.test(expression ~ risk, data = da)
  pValue = test$p.value
  if(!is.na(pValue) & pValue < 0.05){
    p <- ggviolin(da, x = "risk", y = "expression", fill = "risk",
                  palette = c("#DC0000FF", "#4DBBD5FF"), add = "boxplot",
                  add.params = list(fill="white"))+
      ylab(paste(colname, "sensitivity (IC50)", sep = " "))+
      stat_compare_means(label = "p.format",label.x = 1.25, size = 5)+
      theme(text = element_text(size = 20))
    #p
    ggsave(filename = paste(colname, "tiff", sep = "."), plot = p, dpi = 500,
           width = 6, height = 5)
    ggsave(filename = paste(colname, ".pdf", sep = ""), plot = p, dpi = 400,
           width = 9, height = 8)
  }
}

# t <- as.data.frame(colnames(drug_sen))
# write.table(t, file = "drug_list.txt", sep = "\t", quote = F, col.names = T)




#################################################################################
##########################  SAVE LIST  ########################################
#################################################################################

rm(list=ls())
setwd("E:/LZ/24105")
if (!dir.exists("19_drug")){dir.create("19_drug")}
setwd("19_drug")
# 设置文件夹路径
folder_path <- "E:/LZ/24105/19_drug"
# 获取文件夹下所有 .tiff 文件的完整路径
tiff_files <- list.files(path = folder_path, pattern = "\\.tiff$", full.names = TRUE)
# 提取文件名（不包含路径）
file_names <- basename(tiff_files)
# 如果需要去掉文件扩展名
file_names_no_ext <- tools::file_path_sans_ext(file_names)
# 输出文件名
print(file_names_no_ext)
write.table(file_names_no_ext,"drugnames.txt",sep="\t",quote = F,row.names = F,col.names = F)

#################################################################################
################################# cellminer #####################################
#################################################################################
#Drug Sensitivity Analysis
rm(list=ls())
setwd("E:/LZ/24105")
if (!dir.exists("19_drug")){dir.create("19_drug")}
setwd("19_drug")
library(impute)
library(limma)

rt <- read.table("../00_rawdata/drug.txt", sep = "\t", header = T, check.names = F)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
drug <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(drug),colnames(drug))
data <- matrix(as.numeric(as.matrix(drug)), nrow = nrow(drug), dimnames = dimnames)
mat<- impute.knn(data) # 填补缺失值
drug <- mat$data
drug <- avereps(drug)

#Read expression input file
exp <- read.table("geneExp.txt", sep= "\t", header = T, row.names = 1,
                  check.names = F) #固定的文件
dim(exp)
# [1] 23808    60
exp[1:4, 1:4]

#Extract specific gene expression
gene <- read.csv("../07_multicox/multiCox.csv") #your genes
genelist <- as.vector(gene[, 1])
genelist
genelist <- gsub(" ", "", genelist)
genelist <- intersect(genelist, row.names(exp))
exp <- exp[genelist, ]

#Drug Sensitivity Calculation
#drug <- as.data.frame(drug)
outTab <- data.frame()
for (Gene in row.names(exp)) {
  x <- as.numeric(exp[Gene, ])
  for (Drug in row.names(drug)) {
    y <- as.numeric(drug[Drug, ])
    corT <- cor.test(x, y, method = "spearman")
    cor <- corT$estimate
    pvalue <- corT$p.value
    if(!is.na(pvalue) & pvalue < 0.05) {
      outVector <- cbind(Gene, Drug, cor,pvalue)
      outTab <- rbind(outTab, outVector)
    }
  }
}
outTab <- outTab[order(as.numeric(as.vector(outTab$pvalue))), ]
write.table(outTab, file= "drugCor.txt", sep= "\t", row.names=F, quote=F)

#outTab_1 <- outTab[outTab$cor>0, ]
#visualization
library(ggplot2)
library(ggpubr)
#setwd("E:\\LiZhen\\CESC\\14drug\\02CellMiner")
# DGIDB:Vandetanib
outTab <- read.table("drugCor.txt", sep = "\t",header = T,check.names = F)
outTab1 <- outTab %>% filter(Drug %in% c('Sorafenib', 'Lenalidomide', 'Regorafenib', 
                                         'Apatinib', 'Cabozantinib',"Vandetanib",
                                         "Bortezomib","XL-184","Sunitinib"))
                             # & Gene %in% genelist)
# lenvatinib:30,243
# dgidb:搜索获批药物
# CDX2: 1.TEGAFUR:80
# NEK2: 1.PAZOPANIB:188 2.PALBOCICLIB: 146
outTab <- outTab1
#1.Scatter plot
plotList_1 <- list()
corPlotNum <- 5
if(nrow(outTab) < corPlotNum){
  corPlotNum = nrow(outTab)
}
for (i in 1:corPlotNum) {
  Gene <- outTab[i, 1]
  Drug <- outTab[i, 2]
  
  x <- as.numeric(exp[Gene, ])
  y <- as.numeric(drug[Drug, ])
  cor <- sprintf( "%.03f", as.numeric(outTab[i, 3]))
  pvalue = 0
  if(as.numeric(outTab[i, 4]) < 0.001){
    pvalue = "p<0.001"
  } else {
    pvalue = paste0("p=", sprintf( "%.03f", as.numeric(outTab[i, 4])))
  }
  if(cor != 0){
    df1 <- as.data.frame(cbind(x, y))
    p1 = ggplot(data = df1, aes(x = x, y = y))+
      geom_point(size = 2)+
      stat_smooth(method = "lm", se = FALSE, formula = y~x)+
      labs(x = "Expression", y = "IC50", title = paste0(Gene, ", ", Drug),
           subtitle = paste0("Cor=", cor, ", ", pvalue))+
      theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20))
    plotList_1[[i]] = p1
  }
}

#2.boxplot
library(ggsci)
plotList_2 <- list()
corPlotNum <- 5
if(nrow(outTab) < corPlotNum){
  corPlotNum = nrow(outTab)
}
for (i in 1:corPlotNum) {
  Gene <- outTab[i, 1]
  Drug <- outTab[i, 2]
  x <- as.numeric(exp[Gene, ])
  y <- as.numeric(drug[Drug, ])
  df1 <- as.data.frame(cbind(x,y))
  colnames(df1)[2] <- "IC50"
  df1$group <- ifelse(df1$x > median(df1$x), "high", "low")
  compaired <- list(c("low", "high"))
  p1 = ggboxplot(df1, x = "group", y = "IC50", fill = "group",
                 palette = c( "#E64B35FF", "#0072B5FF"), add = "jitter", size = 0.5,
                 xlab = paste0("The_expression_of_", Gene),
                 ylab = paste0( "IC50_of_", Drug))+
    stat_compare_means(comparisons = compaired, method = "wilcox.test",
                       symnum.args= list(cutpoints = c( 0, 0.001, 0.01, 0.05, 1),
                                         symbols = c( "***", "**", "*", "ns")))+
    theme(text = element_text(size = 15))
  plotList_2[[i]] = p1
}

nrow <- ceiling(sqrt(corPlotNum))
ncol <- ceiling(corPlotNum/nrow)
p1 <- ggarrange(plotlist = plotList_1, nrow = nrow, ncol = ncol)
p2 <- ggarrange(plotlist = plotList_2, nrow = nrow, ncol = ncol)
ggsave(filename = "plotList_1.tiff", plot = p1, dpi = 500, width = 8, height = 12)
ggsave(filename = "plotList_1.pdf", plot = p1, dpi = 500, width = 8, height = 12)
ggsave(filename = "plotList_2.tiff", plot = p2, dpi = 500, width = 18, height = 18)
ggsave(filename = "plotList_2.pdf", plot = p2, dpi = 500, width = 18, height = 18)




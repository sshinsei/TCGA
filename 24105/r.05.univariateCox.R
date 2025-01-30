rm(list=ls())
setwd("E:/LZ/24105")  #modify work path
if(! dir.exists("05_unicox")){
  dir.create("05_unicox")
}
setwd("./05_unicox")
########--------mergedata--------######
library(tidyverse)
##----------要删除正常样本----------------##
norExp = read.table(file = "../00_rawdata/GSE14520/GSE14520.txt", 
                    sep = "\t", header = T, check.names = F) # 430,有一列是id
cli_data = read.table(file = "../00_rawdata/GSE14520/survival.csv", 
                      sep = ",", header = T)
cli_data = cli_data[,c(1,3,2)]
# fea_genes = read.table(file = "../03_upset/feature_genes.txt", sep = "\t", header = F, check.names = F) # 11
fea_genes = read.table(file = "../04_clustDEG/feature_genes.txt", sep = "\t", header = F, check.names = F) # 11
fea_genes$V1 <- unique(fea_genes$V1)
rows = na.omit(match(fea_genes$V1, rownames(norExp))) #69
norExp_1 = norExp[rows, ]
unique(rownames(norExp_1) == fea_genes$V1)  #check
#rownames(norExp_1) = rownames(norExp)
# norExp_1 = norExp_1[ , -1]
norExp_1 = t(norExp_1)
norExp_1 = as.data.frame(norExp_1)
norExp_1$id = rownames(norExp_1)
#rows = na.omit(match(rownames(norExp_1), cli_data$Id))
#norExp_2 = cbind(norExp_1, cli_data[rows, ])
#norExp_2 = norExp_2[ , c((ncol(norExp_2)-2):ncol(norExp_2), 1:(ncol(norExp_2)-3))]
#unique(rownames(norExp_2) == norExp_2$Id)  #check

norExp_1$id <- substring(norExp_1$id, 1, 16) # 488


norExp_2 = merge(cli_data, norExp_1, by.x = "sample", by.y = "id") # 488
norExp_2 <- na.omit(norExp_2) # 242
# 删除正常样本
group <- read.table(file = "../00_rawdata/GSE14520/GSE14520.group2.txt", 
                    sep = "\t", header = T, check.names = F) 
#case <- group[group$group == "Case",]$sample
#!norExp_2$sample %in% case

norExp_2 = norExp_2[norExp_2$sample %in% group[group$group == "Case",]$sample,] # 242
# 去重
norExp_2 <- norExp_2 %>% 
  distinct(sample, .keep_all = TRUE) # 242
write.table(norExp_2, file = "../00_rawdata/GSE14520/tumor_clinical.txt", 
            sep = "\t", row.names = F, quote = F)

library(survival)
rt = read.table("../00_rawdata/GSE14520/tumor_clinical.txt", header = T, sep = "\t",
                check.names = F,row.names = 1)
rt = rt[rt$futimes > 1, ] # 344
rt[ , "futimes"] = rt[ , "futimes"]/12
#3/7groups
#train_sub <- sample(nrow(rt), floor(7/10 * nrow(rt)))
#train_data = rt[sort(train_sub), ]  #训练集
#test_data = rt[-train_sub, ]  #测试集
#write.table(train_data, file = "train.txt", sep = "\t", col.names = T, quote = F)
#write.table(test_data, file = "test.txt", sep = "\t", col.names = T, quote = F)

#rt = read.table("train.txt", header = T, sep = "\t", row.names = 1, check.names = F)
#rt[ , 3:ncol(rt)] = log2(rt[ , 3:ncol(rt)] + 0.001)
outTab=data.frame()
for(i in colnames(rt[ ,3:ncol(rt)])){
  cox <- coxph(Surv(futimes, fustate) ~ rt[ ,i], data = rt)
  coxSummary = summary(cox)
  outTab = rbind(outTab, cbind(gene = i, HR = coxSummary$coefficients[ , "exp(coef)"],
                               HR.95L = coxSummary$conf.int[ , "lower .95"],
                               HR.95H = coxSummary$conf.int[ , "upper .95"],
                               z = coxSummary$coefficients[ , "z"],
                               pvalue = coxSummary$coefficients[ , "Pr(>|z|)"]))
}
write.table(outTab, file = "univariateCox.csv", sep = ",", row.names = F,
            quote = F)
write.table(outTab[outTab$pvalue < 0.05, ], file = "univariateCox_P0.05.csv",
            sep = ",", row.names = F, quote = F)
write.table(outTab[outTab$pvalue < 0.01, ], file = "univariateCox_P0.01.csv",
            sep = ",", row.names = F, quote = F)
write.table(outTab[outTab$pvalue < 0.001, ], file = "univariateCox_P0.001.csv",
            sep = ",", row.names = F, quote = F)

######绘制森林图######
rm(list=ls())
setwd("E:/LZ/24105")  #modify work path
if(! dir.exists("05_unicox")){
  dir.create("05_unicox")
}
setwd("./05_unicox")
########--------mergedata--------######
library(tidyverse)
#读取输入文件
# 0.05=10 0.01=8 0.001=3
# 0.05=8 0.01=7 0.001=2
# 0.05=42 0.01=28 0.001=3
rt <- read.table("univariateCox_P0.01.csv",header=T,sep=",",row.names=1,check.names=F)
#rt = rt[rt$pvalue < 0.0001, ]
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#输出图形
pdf(file="forest.pdf", width = 7, height = 8)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1, nRow)
layout(matrix(c(1, 2), nc=2),width=c(3, 2.7))

#绘制森林图左边的临床信息
xlim = c(0, 3)
par(mar=c(4, 2.5, 2, 1))
plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab="", ylab="")
text.cex=0.9
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',
                                                   cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',
                                                 cex=text.cex,font=2,adj=1,)

#绘制森林图
par(mar=c(4,1,2,1), mgp=c(2,0.5,0))
xlim = c(-0.1, max(as.numeric(hrLow), as.numeric(hrHigh)))
plot(1, xlim=c(0.5,2.0), ylim=ylim, type="n", axes=F, ylab="", xaxs="i",
     xlab = "Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,
       length=0.05,col="#3948ABFF",lwd=1.5)
abline(v=1,col="black",lty=2,lwd=1.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED0000FF", "#42B540FF")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(side = 1)
dev.off()

tiff(file="forest.tiff",width=6,height=7,units ="in",compression="lzw",
     bg="white",res=400)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3, 2.7))

#绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',
                                                   cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',
                                                 cex=text.cex,font=2,adj=1,)

#绘制森林图(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(-0.1,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=c(0.5,2.0),ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,
       length=0.05,col="#3948ABFF",lwd=1.5)
abline(v=1,col="black",lty=2,lwd=1.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED0000FF", "#42B540FF")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()





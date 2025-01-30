#############################-----------GEO---------------------#################
rm(list=ls())
setwd('E:/LZ/24105/10_GEO/')
if(! dir.exists("GSE76427")){dir.create("GSE76427")}
setwd("GSE76427")

###------------merge data---------------
cli = read.table("../../00_rawdata/GSE76427/GSE76427.clinic.txt", header = T, sep = "\t", check.names = FALSE)
nor = read.csv("../../00_rawdata/GSE76427/GSE76427.txt", header = T, sep = "\t", check.names = FALSE,
               row.names = 1)
# 多因素cox后的结果的geneCoef
gene = read.table("../../07_multicox/multiCox.csv", header = T, sep = ",", check.names = FALSE)
#gene <- c("LY6H", "KCNE1L", "RNASEH2A")
#gene <- data.frame(id=gene)
rows = na.omit(match(gene$id, rownames(nor)))
nor_1 = t(nor[rows, ])
#colnames(nor_1) = nor_1[1, ]
#nor_1 = nor_1[-1, ]
nor_1 = as.data.frame(nor_1)
nor_1$id = rownames(nor_1)
data = merge(cli, nor_1, by.x = "sample", by.y = "id")
write.table(data, file = "multiInput.txt", sep = "\t", quote = F, row.names = F)

rt = read.csv("./multiInput.txt", 
              header = T, sep = "\t", check.names = F,
              row.names = 1)
# rt[ , "futimes"] = rt[ , "futimes"]/365 #futime本身是日：除以365
# rt = rt[rt$futimes > 30, ] # ?
rt <- na.omit(rt)
#rt[ , 3:ncol(rt)] = log2(rt[ , 3:ncol(rt)] + 0.01)
#rt = rt[ , c("futime", "fustat", "ID1", "TNFAIP2", "CXCR4", "SERPINE1", "ANXA5")]
cox <- coxph(Surv(futimes, fustate) ~ ., data = rt)
#cox = step(cox,direction = "both")
#rt_1 = read.table("multiInput_1.txt", header = T, sep = "\t", check.names = F,
#                row.names = 1)
#rt_1$futime = (rt_1$futime)/12
riskScore = predict(cox, type = "risk", newdata = rt)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
write.table(cbind(id = rownames(cbind(rt[ , 1:2], riskScore, risk)),
                  cbind(rt[ , 1:2], riskScore, risk)), file = "risk.txt",
            sep = "\t", quote = F, row.names = F)
cox

pdf(file="forest.pdf",
    width = 8,             #ͼƬ?Ŀ???
    height = 5,            #ͼƬ?ĸ߶?
)
ggforest(cox,
         main = "Hazard ratio",
         cpositions = c(0.0,1.0, 1.5), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()

#Output model parameters
multiCoxSum = summary(cox)
outTab = data.frame()
outTab = cbind(
  coef = multiCoxSum$coefficients[ , "coef"],
  HR = multiCoxSum$conf.int[ , "exp(coef)"],
  HR.95L = multiCoxSum$conf.int[ , "lower .95"],
  HR.95H = multiCoxSum$conf.int[ , "upper .95"],
  pvalue = multiCoxSum$coefficients[ , "Pr(>|z|)"])
outTab = cbind(id = row.names(outTab), outTab)
outTab = gsub("`", "", outTab)
write.table(outTab, file = "multiCox.xls", sep = "\t", row.names = F, quote = F)

#############################--------ROC----------###############################
rm(list=ls())
rt <- read.csv("E:/LZ/24105/10_GEO/GSE76427/risk.txt", 
               header = T, sep = "\t", check.names = F, row.names = 1)
ROC_rt <- timeROC(T = rt$futimes, delta = rt$fustate, marker = rt$riskScore,
                  cause = 1, weighting = 'aalen', times = c(1, 3, 5),
                  ROC = TRUE, iid = FALSE)

pdf(file = "ROC_test.pdf", width = 5, height = 5)
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

tiff(file="ROC_test.tiff",width = 5, height = 5, units = "in",
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
       col = c("#DC0000FF", "#00A087FF", "#4DBBD5FF"),  #, "#F39B7FFF"
       lwd = 2, bty = 'n')
dev.off()


###############################----------KM------------------------##############################
rm(list=ls())
library(survival)
library(ggsci)
rt = read.table("./risk.txt", 
                header = T, sep = "\t", row.names = 1)
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




###############################################################################
###############        riskcurve          #####################################
###############################################################################
rm(list=ls())
setwd('E:/LZ/24105/10_GEO/GSE76427/')
data = read.table("./multiInput.txt", 
                  sep="\t", header=T, check.names=F)
risk = read.table("./risk.txt", 
                  sep="\t", header=T, check.names=F) 
#gene = read.table("geneCoef.txt", sep="\t", header=T, check.names=F) 
#cols = match(gene$id, colnames(data))
#data = data[ , c(1, cols)]
rt = merge(data, risk[ , c(1,4,5)], by.x = "sample", by.y = "id")
# rt$futime = rt$futime/12
write.table(rt, file = "risk_1.txt",sep="\t", row.names=F, quote=F)

rt = rt[order(rt$riskScore, decreasing = F), ]
#???Ʒ???ͼ
riskClass = rt[ , "risk"]
lowLength = length(riskClass[riskClass == "low"])
highLength = length(riskClass[riskClass == "high"])
line=rt[,"riskScore"]
line[line>10]=10

#????״̬ͼ
color = as.vector(rt$fustat)
color[color == 1] = "red"
color[color == 0] = "blue"

#???Ʒ???ͼ&????״̬ͼ
tiff(file="riskScore_survStat.tiff", width = 7, height = 3.5, units = "in", bg = "white",
     res = 500)
layout(matrix(c(1, 2), 2, 1), width = 5, heights = c(1.6, 1.6), respect = T)
par(mar=c(0, 4, 1, 1))
plot(rt$futimes,
     pch=19,
     xlab="",
     ylab="Survival time (years)",
     col = color, cex.axis = 0.75, cex.lab = 0.75, cex=0.5, xaxt="n")
abline(v=lowLength,lty=2)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","blue"),cex=0.75)

par(mar=c(2, 4, 0.25, 1))
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("#4DBBD5FF",lowLength),
           rep("#DC0000FF",highLength)),
     cex.axis = 0.75, cex.lab = 0.75)
trainMedianScore=median(rt$riskScore)
abline(h=trainMedianScore,v=lowLength,lty=2)
legend("topleft", c("high risk", "low risk"),bty="n",pch=19,col=c("#DC0000FF","#4DBBD5FF"),
       cex=0.75)
dev.off()

pdf(file="riskScore_survStat.pdf", width = 7, height = 3.5)
layout(matrix(c(1, 2), 2, 1), width = 5, heights = c(2, 2), respect = T)
par(mar=c(0, 4, 1, 1))
plot(rt$futimes,
     pch=19,
     xlab="",
     ylab="Survival time (years)",
     col = color, cex.axis = 0.75, cex.lab = 0.75, cex=0.5, xaxt="n")
abline(v=lowLength,lty=2)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","blue"),cex=0.75)

par(mar=c(2, 4, 0.25, 1))
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("#4DBBD5FF",lowLength),
           rep("#DC0000FF",highLength)),
     cex.axis = 0.75, cex.lab = 0.75)
trainMedianScore=median(rt$riskScore)
abline(h=trainMedianScore,v=lowLength,lty=2)
legend("topleft", c("high risk", "low risk"),bty="n",pch=19,col=c("#DC0000FF","#4DBBD5FF"),
       cex=0.75)
dev.off()

#??????????ͼ
library(pheatmap)
#rt=rt[order(rt$riskScore), ]
table(rt$risk)
# high  low 
#  57   58 
#rt <- rt[,c(2,3,5,6,8,10,11,1,4,7,9)]
type = c(rep("low", 57), rep("high", 58))
#gene = read.table("gene.txt",sep = "\t",header = T)
#rt1 = rt[row.names(rt) %in% gene$id,] 
rt1=rt[ , c(4:(ncol(rt)-2))] #只留下gene的表达矩阵
rt1=t(rt1)
#rt1=log2(rt1+0.01)
cluster <- data.frame(type=factor(type))
rownames(cluster) <- colnames(rt1)
ann_color <- list(type = c(low = "#4DBBD5FF", high = "#DC0000FF"))
pdf("heatmap.pdf", height = 2.5, width = 4.5)
#annotation_col = data.frame(Subtype=factor(c(rep("Subtype_1",33),
#rep("Subtype_2",14),rep("Subtype_3",2))))
#row.names(annotation_col) = rownames(cluster)
pheatmap(rt1, annotation=cluster, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=6,
         scale="row",
         show_colnames=F,
         fontsize_col=3,
         annotation_colors = ann_color)
dev.off()

tiff(file="heatmap.tiff", height = 2.5, width = 4.5, units = "in",
     bg = "white", res = 500)
pheatmap(rt1, annotation=cluster, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         #cluster_rows = F,
         fontsize=8,
         fontsize_row=6,
         scale="row",
         show_colnames=F,
         fontsize_col=3,
         annotation_colors = ann_color)
dev.off()














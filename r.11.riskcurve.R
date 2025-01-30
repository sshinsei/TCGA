rm(list=ls())
setwd("E:/LZ/24080")  #modify work path
if(! dir.exists("11_riskcurve")){
  dir.create("11_riskcurve")
}
setwd("./11_riskcurve")
data = read.table("../00_rawdata/tumor_clinical.txt", 
                  sep="\t", header=T, check.names=F)
rt=read.table("../00_rawdata/risk.txt", 
              sep="\t", header=T, check.names=F) 
gene = read.csv("../06_multicox/multiCox.csv", header=T, check.names=F) 
cols = match(gene$id, colnames(data))
data = data[ , c(1, cols)]
rt1 = merge(data, rt, by.x = "sample", by.y = "id")
# write.table(rt1[ , c(1, 12, 13, 2:11, 14, 15)], file = "risk_1.txt",sep="\t",
# row.names=F, quote=F)
write.table(rt1, file = "risk_1.txt",sep="\t",
            row.names=F, quote=F)

########################-------plot----------------------#########################
rm(list = ls())
library(pheatmap)
rt = read.table("risk_1.txt",sep = "\t",header = T,row.names = 1)
rt=rt[order(rt$riskScore), ]
table(rt$risk)
# high  low 
# 172  172
#rt <- rt[,c(2,3,5,6,8,10,11,1,4,7,9)]
type = c(rep("low", 172), rep("high", 172))
#gene = read.table("gene.txt",sep = "\t",header = T)
#rt1 = rt[row.names(rt) %in% gene$id,] 
# rt1=rt[ , c(3:(ncol(rt)-2))]
rt1=rt[ , 1:7] # 只保留基因的表达矩阵
rt1=t(rt1)
#rt1=log2(rt1+0.01)

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
plot(rt$futime,
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
plot(rt$futime,
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


#####################---单基因KM曲线--------######################################
rm(list=ls())
setwd("E:/LZ/24080")  #modify work path
if(! dir.exists("11_riskcurve")){
  dir.create("11_riskcurve")
}
setwd("./11_riskcurve")
outTab=data.frame()

picDir="picture"
dir.create(picDir)

library(survival)
library(qvalue)
rt=read.table("./risk_1.txt",
              header=T,sep="\t",row.names=1,check.names=F)
rt = rt[ , -c(ncol(rt) - 1, ncol(rt))]
#rt[,"futime"]=rt[,"futime"]/365

# 要先新建一个picture的文件夹才能运行
if(! dir.exists("picture")){
  dir.create("picture")
}
# p>0.05的不画图
for(i in colnames(rt[,1:7])){
  cox <- coxph(Surv(futimes, fustate) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  med=median(rt[,i])
  if(med!= 0){
    a=rt[,i]>med
    rt1=rt[a,]
    b=setdiff(rownames(rt),rownames(rt1))
    rt2=rt[b,]
    n1=nrow(rt1)
    n2=nrow(rt2)
    surTab1=summary(survfit(Surv(futimes, fustate) ~ 1, data = rt1))
    surTab2=summary(survfit(Surv(futimes, fustate) ~ 1, data = rt2))
    medianTab1=surTab1$table
    medianTab2=surTab2$table
    diff=survdiff(Surv(futimes, fustate) ~a,data = rt)
    fit <- survfit(Surv(futimes, fustate) ~ a, data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    outTab=rbind(outTab,cbind(gene=i,coxSummary$coefficients,coxSummary$conf.int,KM=pValue,
                              H_med=medianTab1["median"],H_0.95LCL=medianTab1["0.95LCL"],H_0.95UCL=medianTab1["0.95UCL"],
                              L_med=medianTab2["median"],L_0.95LCL=medianTab2["0.95LCL"],L_0.95UCL=medianTab2["0.95UCL"]))
    pval=0
    if(pValue<0.05){
      pval=signif(pValue,4)
      pval=format(pval, scientific = TRUE)
    }else{
      pval=round(pValue,3)
    }
    if(TRUE){
      geneName=unlist(strsplit(i,"\\|",))[1]
      tiffFile=paste(geneName,".survival.tiff",sep="")
      outTiff=paste(picDir,tiffFile,sep="\\")
      tiff(file=outTiff,width = 15,height = 15,units ="cm",compression="lzw",bg="white",res=600)
      plot(fit, col=c("#4DBBD5FF", "#DC0000FF"), xlab="Time (years)", ylab="Overall survival",
           main=paste(geneName,"(p=",pval, ")", sep=""),mark.time=T,ylim=c(0,1.1),
           lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
      legend("topright", c(paste("Low expression"), 
                           paste("High expression")), 
             col=c("#4DBBD5FF", "#DC0000FF"), bty="n", lwd = 2, cex=0.8)
      dev.off()
    }
  }
}





###-----------------pdf-------------------------
for(i in colnames(rt[,1:7])){
  cox <- coxph(Surv(futimes, fustate) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  med=median(rt[,i])
  if(med!=0){
    a=rt[,i]>med
    rt1=rt[a,]
    b=setdiff(rownames(rt),rownames(rt1))
    rt2=rt[b,]
    n1=nrow(rt1)
    n2=nrow(rt2)
    surTab1=summary(survfit(Surv(futimes, fustate) ~ 1, data = rt1))
    surTab2=summary(survfit(Surv(futimes, fustate) ~ 1, data = rt2))
    medianTab1=surTab1$table
    medianTab2=surTab2$table
    diff=survdiff(Surv(futimes, fustate) ~a,data = rt)
    fit <- survfit(Surv(futimes, fustate) ~ a, data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    outTab=rbind(outTab,cbind(gene=i,coxSummary$coefficients,coxSummary$conf.int,KM=pValue,
                              H_med=medianTab1["median"],H_0.95LCL=medianTab1["0.95LCL"],H_0.95UCL=medianTab1["0.95UCL"],
                              L_med=medianTab2["median"],L_0.95LCL=medianTab2["0.95LCL"],L_0.95UCL=medianTab2["0.95UCL"]))
    pval=0
    if(pValue<0.05){
      pval=signif(pValue,4)
      pval=format(pval, scientific = TRUE)
    }else{
      pval=round(pValue,3)
    }
    if(TRUE){
      geneName=unlist(strsplit(i,"\\|",))[1]
      tiffFile=paste(geneName,".survival.pdf",sep="")
      outTiff=paste(picDir,tiffFile,sep="\\")
      pdf(file=outTiff,width = 8,height = 8)
      plot(fit, col=c("#4DBBD5FF", "#DC0000FF"), xlab="Time (years)", ylab="Overall survival",
           main=paste(geneName,"(p=",pval, ")", sep=""),mark.time=T,ylim=c(0,1.1),
           lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
      legend("topright", c(paste("Low expression"), 
                           paste("High expression")), 
             col=c("#4DBBD5FF", "#DC0000FF"), bty="n", lwd = 2, cex=0.8)
      dev.off()
    }
  }
}



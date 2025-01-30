rm(list=ls())
library(glmnet)
library(survival)
setwd("E:/LZ/24105")  #modify work path
if(! dir.exists("06_lasso")){
  dir.create("06_lasso")
}
setwd("./06_lasso")
set.seed(123)
rt = read.table("../00_rawdata/GSE14520/tumor_clinical.txt",
                header = T, sep = "\t", row.names = 1,
                check.names = FALSE)
rt = na.omit(rt)
rt = rt[rt$futimes > 1, ] # 344
#rt[ , 3:ncol(rt)] = log2(rt[ , 3:ncol(rt)] + 0.001)
rt[ , "futimes"] = rt[ , "futimes"]/12 # 344
#rt = na.omit(rt)
# p<0.05的特征基因
gene = read.table("../05_unicox/univariateCox_P0.01.csv", header = T, sep=",") # 47
#cols = match(gene$V1, colnames(rt))
rt = rt[ , c("futimes", "fustate", as.vector(gene[ , 1]))]

x = as.matrix(rt[ , c(3:ncol(rt))])
y = data.matrix(Surv(rt$futimes, rt$fustate))

fit <- glmnet(x, y, family = "cox", maxit = 2000)
#setwd("D:\\jxhe\\2021_11\\liver\\05lasso\\03lasso")
#pdf("lambda.pdf")

tiff("lambda.tiff", width = 18, height =18, units = "cm", compression = "lzw",
     bg = "white", res = 400)
plot(fit, xvar = "lambda", label = TRUE)
abline(v = c(fit$lambda.min, fit$lambda.1se), lty = "dashed")
dev.off()


pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
abline(v = c(fit$lambda.min, fit$lambda.1se), lty = "dashed")
dev.off()

cvfit <- cv.glmnet(x, y, family = "cox", maxit = 2000)  #Cross-validation for glmnet


tiff("cvfit.tiff", width = 18, height = 18, units = "cm", compression = "lzw",
     bg = "white", res = 400)
plot(cvfit)
abline(v = log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty = "dashed")
dev.off()


pdf("cvfit.pdf")
plot(cvfit)
abline(v = log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty = "dashed")
dev.off()


coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene = row.names(coef)[index]
geneCoef = cbind(Gene = lassoGene, Coef = actCoef)

write.table(geneCoef, file = "geneCoef.txt", sep = "\t", quote = F, row.names = F)
# 2
#save.image("lasso.RData")

if(F){
  riskScore=predict(cvfit, newx = x, s = "lambda.min",type="response")
  outCol=c("futime","fustat",lassoGene)
  risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
  outTab=cbind(rt[,1:2],riskScore=as.vector(riskScore),risk)
  write.table(cbind(id=rownames(outTab),outTab),
              file="../00_rawdata/risk.txt",
              sep="\t",
              quote=F,
              row.names=F)
}



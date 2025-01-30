rm(list=ls())
setwd("E:/LZ/24105")  #modify work path
if(! dir.exists("12_nomogram")){
  dir.create("12_nomogram")
}
setwd("./12_nomogram")

cli = read.csv("../00_rawdata/clinical.csv", header=T)
colnames(cli)[1] <- "sample"
colnames(cli)[9] <- "id"
cli <- cli[,-1]
risk = read.table("../00_rawdata/risk.txt", header=T, sep="\t", check.names=F)
#risk = risk[ , -c(2:3)]
# data <- left_join(risk,cli,by="id")
data = merge(cli, risk, by.x = 'id', by.y = 'id')
#indep_1 = indep[ , c(1:8, 12)]
#data[which(data$fustat == "Alive"), 2] = 0
#data[which(data$fustat == "Dead"), 2] = 1
#data$age = as.numeric(data$age)
#data[ , 4] = round(data[ , 4]/365)
data = na.omit(data) # 271
#data = data[ , -5]
write.table(data, file = "indepInput.txt", sep = "\t", quote = F, col.names = T)


###############################------uni-------##################################
rm(list=ls())
library(survival)
library(tidyverse)
rt=read.csv("./indepInput.txt", header=T,sep="\t",check.names=F,
            row.names = 1)
#rt = na.omit(rt)
###------------------每个变量都要处理为2分类变量------------------###
data <- rt
# age:分为>=65和<65
data$age <- ifelse(data$age >= 65, ">=65","<65")
data$age <- as.factor(data$age)

data <- data %>% 
  mutate(stage = case_when(grepl("III|IV", data$stage) ~ "Stage III+IV",
                           grepl("I|II", data$stage) ~ "Stage I+II",
                           TRUE ~ NA),
         `T` = ifelse(grepl("T1|T2",data$T), 
                      "T1+T2","T3+T4"),
         M = case_when(M=="M0"  ~ "M0",
                       M=="M1"  ~ "M1",
                       TRUE ~ NA),
         N = case_when(N=="N0"  ~ "N0",
                       N=="N1"  ~ "N1",
                       TRUE ~ NA)) %>% 
  drop_na()

data$stage <- as.factor(data$stage)
data$`T` <- as.factor(data$`T`)
data$M <- as.factor(data$M)
data$N <- as.factor(data$N)
data$gender <- as.factor(data$gender)

# select(-age) # age不显著，不影响，要看的只是列线图比单个临床因素的效果好或坏
write.table(data,file="indepInput_clean.txt",sep="\t",row.names=F,quote=F,
            fileEncoding = "UTF-8")

outTab = data.frame()
for(i in colnames(data[ , c(2, 6, 7, 8, 11)])){ # stage M N T riscScore
  da = data[ , c('futimes', 'fustate', i)] %>% na.omit()
  cox <- coxph(Surv(futimes, fustate) ~ da[ , i], data = da)
  coxSummary = summary(cox)
  coxP = coxSummary$coefficients[ , "Pr(>|z|)"]
  outTab = rbind(outTab,
                 cbind(id = i,
                       HR = coxSummary$conf.int[ , "exp(coef)"],
                       HR.95L = coxSummary$conf.int[ , "lower .95"],
                       HR.95H = coxSummary$conf.int[ , "upper .95"],
                       pvalue = coxSummary$coefficients[ , "Pr(>|z|)"]))
}

write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)

######????ɭ??ͼ######
#??ȡ?????ļ?
rt <- read.table("uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#????ͼ??
pdf(file="uni_forest.pdf", width = 7,height = 4)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

#????ɭ??ͼ???ߵ??ٴ???Ϣ
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#????ɭ??ͼ
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,
       length=0.05,col="#3948ABFF",lwd=1.5)
abline(v=1,col="black",lty=2,lwd=1.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED0000FF", "#42B540FF")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

#tiff
tiff(file="uni_forest.tiff", width = 7,height = 4, units = "in", bg = "white",
     res = 500)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

#????ɭ??ͼ???ߵ??ٴ???Ϣ
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#????ɭ??ͼ
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,
       length=0.05,col="#3948ABFF",lwd=1.5)
abline(v=1,col="black",lty=2,lwd=1.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED0000FF", "#42B540FF")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

###########################-----------multi-----------------####################
rm(list=ls())
setwd("E:/LZ/24105")  #modify work path
if(! dir.exists("12_nomogram")){
  dir.create("12_nomogram")
}
setwd("./12_nomogram")
rm(list=ls())
rt = read.table("./indepInput_clean.txt",
                header=T,sep="\t",check.names=F,row.names=1)
rt = na.omit(rt)
rt = rt[ , -c(2,3,4,11)] # 删除grade M N和risk
# 转为factor
rt[, 1:4] <- lapply(rt[, 1:4], as.factor) # age T stage 转为FACTOR

multiCox=coxph(Surv(futimes, fustate) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
rownames(outTab) <- c("age","T","stage","gender","riskScore")
# outTab=cbind(id=c("age","T","stage","riskScore"),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=T,quote=F)

######????ɭ??ͼ######
#??ȡ?????ļ?
rt <- read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#????ͼ??
pdf(file="multi_forest.pdf", width = 7,height = 4)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

#????ɭ??ͼ???ߵ??ٴ???Ϣ
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#????ɭ??ͼ
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,
       length=0.05,col="#3948ABFF",lwd=1.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED0000FF", "#42B540FF")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

#tiff
tiff(file="multi_forest.tiff", width = 7,height = 4, units = "in", bg = "white",
     res = 500)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

#????ɭ??ͼ???ߵ??ٴ???Ϣ
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#????ɭ??ͼ
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,
       length=0.05,col="#3948ABFF",lwd=1.5)
abline(v=1,col="black",lty=2,lwd=1.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED0000FF", "#42B540FF")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()


######################------------NOMOGRAM----------------######################
rm(list=ls())
library(rms)
library(dplyr)
rt=read.table("./indepInput_clean.txt",
              sep="\t",header=T,row.names=1,check.names=F) %>% na.omit()
#rt$futime <- rt$futime/365
#rt$gender <- ifelse(rt$gender == "0","F","M")
rt = rt[ , -c(2,3,4,11)] # 删除grade和risk
# 转为factor
rt[, 1:4] <- lapply(rt[, 1:4], as.factor) # age T stage 转为FACTOR


dd <- datadist(rt)
options(datadist="dd")
# rt$futime=as.integer(rt$futime)
# rt$age[which(rt$age == "??65")] <- "<=65"
# rt$age[which(rt$age == ">65")] <- ">65"
# rt$gender <- factor(rt$gender,labels = c("F", "M"))
rt$age<-factor(rt$age,labels=c(">=65","<65"))
rt$T<-factor(rt$T,labels=c("T1+T2","T3+T4"))
# rt$N <- factor(rt$N ,labels=c("N0","N1"))
# rt$M <- factor(rt$M ,labels=c("M0","M1"))
rt$stage <- factor(rt$stage,labels = c("Stage I+II","Stage III+IV"))
rt$gender <- factor(rt$gender, labels = c("female","male"))
#???ɺ???
f <- cph(Surv(futimes, fustate) ~ age+stage+`T`+gender+riskScore, x=TRUE, y=TRUE,
         surv=TRUE, data=rt)
surv <- Survival(f)
# 
# #??��nomogram
# nom <- nomogram(f, fun=list(function(x)surv(3,x),function(x)surv(5,x) ),
#                 lp=F, funlabel=c("3-year survival","5-year survival"),
#                 maxscale=100,
#                 fun.at=c(0.99, 0.8, 0.6, 0.4, 0.2, 0.05))
if(F){
  library(regplot)
  nom <- regplot(f,
                 # observation=pbc[1,],
                 obscol = "#326db1",
                 failtime = c(1,3,5), 
                 plots = c("violin","bars"),
                 droplines = T, # 是否画竖???
                 points = T,
                 title = "Nomogram", # 更换标题
                 # odds = T, # 是否显示OR???
                 showP = T, # 是否显示变量的显著性标记（默认：T???
                 rank = "sd", # 根据sd给变量排???
                 #interval="confidence", # 展示可信区间
                 clickable = F, # 是否可以交互（默认：F???
                 prfail = T)
  tiff(filename = "nomogram.tiff", width = 13, height = 10, units = "in",
       res = 500)
  regplot(f,
          # observation=pbc[1,],
          obscol = "#326db1",
          failtime = c(1,3,5), 
          plots = c("violin","bars"),
          droplines = T, # 是否画竖???
          points = T,
          title = "Nomogram", # 更换标题
          # odds = T, # 是否显示OR???
          showP = T, # 是否显示变量的显著性标记（默认：T???
          rank = "sd", # 根据sd给变量排???
          #interval="confidence", # 展示可信区间
          clickable = F, # 是否可以交互（默认：F???
          prfail = T)
  dev.off()
  
  pdf(file = "nomogram.pdf", width = 13, height = 10)
  plot(nom)
  dev.off()
}
library(ggplot2)
#ggsave(nom, plot = "test.pdf")
#regplot(f,
#failtime = c(3, 5), prfail = T, droplines=T)
nom <- nomogram(f,fun = list(function(x)surv(1,x),function(x)surv(3,x),
                             function(x)surv(5,x)), lp=T,
                funlabel = c('1-year survival Probability',
                             '3-year survival Probability',
                             '5-year survival Probability'), maxscale = 100,
                fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
pdf("nomogram.pdf", width = 13, height = 10)
plot(nom, lplabel="Linear Predictor", xfrac=.35,varname.label=TRUE,
     varname.label.sep="=", ia.space=.2, tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,
     lwd=5, label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()

tiff(filename = "nomogram.tiff", width = 13, height = 10, units = "in",
     res = 400)
plot(nom, lplabel="Linear Predictor", xfrac=.35,varname.label=TRUE,
     varname.label.sep="=", ia.space=.2, tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,
     lwd=5, label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()
# pdf(file="nomogram.pdf", height = 8, width = 11)
# nom
# dev.off()

############DCA曲线
# plot_decision_curve(f,
#                     col = c('red'),
#                     confidence.intervals =FALSE)

# devtools::install_github('yikeshu0611/ggDCA')

####--------------------DCA---------------------
library(ggDCA)
library(survival)
dca <- dca(f, times=c(1*1, 3*1, 5*1))
pdf(file="DCA.pdf", width = 7, height = 6) 
ggplot(dca)
dev.off()

tiff(filename = "DCA.tiff", width = 7, height = 6, units = "in", res = 600)
ggplot(dca)
dev.off()
#################################
t <- c(1,3,5)
d <- dca(f, times = t)
ggplot(d)

library(ggsci)
myfunction <- function(string) {
  # print(var)
  print(string)
  result <- paste(as.character(string),' Year', sep=" ")
  return(result)
}
d$time <- myfunction(d$time)
pdf("dca_test.pdf",height = 4,width = 10)
ggplot(dca, linetype = F) + scale_color_jama(name="Model Type") + 
  theme_bw()+
  facet_wrap(~time)
dev.off()

#2024.5.8-----------校准曲线-----------------------
library(rmda)
# form <- as.formula()
# dca <- decision_curve()
################################
#У׼ͼ
#rt$futime = rt$futime/365
cox1 <- cph(Surv(futimes, fustate) ~ age+`T`+stage+riskScore, surv=T, x=T,
            y=T, time.inc = 1, data = rt) 
#B means samples number, m=B/3
cal1 <- calibrate(cox1, cmethod="KM", method="boot", u=1, m=nrow(rt)/3, 02513)
cox2 <- cph(Surv(futimes, fustate) ~ age+`T`+stage+riskScore, surv=T, x=T,
            y=T, time.inc = 3, data = rt) 
cal2 <- calibrate(cox2, cmethod="KM", method="boot", u=3, m=nrow(rt)/3, 02513)  #B means samples number, m=B/3
cox3 <- cph(Surv(futimes, fustate) ~ age+`T`+stage+riskScore, surv=T, x=T,
            y=T,time.inc = 5,data=rt) 
cal3 <- calibrate(cox3, cmethod="KM", method="boot", u=5, m=nrow(rt)/3, 02513)
# cox4 <- cph(Surv(futime, fustat) ~ `T`+N+M+stage+riskScore, surv=T, x=T,
# y=T,time.inc = 10, data=rt) 
# cal4 <- calibrate(cox4, cmethod="KM", method="boot", u=10, m=nrow(rt)/3, 02513)
# nomgram = predict(multiCox,type="risk",newdata=rt)
pdf("calibrate.pdf", width = 11,  height = 8)
par(mar = c(10,5,3,2),cex = 1.0)
plot(cal1,lwd=3,lty=2,errbar.col="#58aa5a",xlim = c(0,1),ylim = c(0,1),
     xlab ="Nomogram-Predicted Probability of Survival", ylab="Actual Survival")
lines(cal1[,c('mean.predicted','KM')],type = 'b',lwd = 3,col ="#58aa5a" ,pch = 16)
# box(lwd = 1)
# abline(0,1,lty = 3,lwd = 3,col = "red")

plot(cal2,lwd=3,lty=2,errbar.col="#f7903c",xlim = c(0,1),ylim = c(0,1),add = T)
lines(cal2[,c('mean.predicted','KM')],type = 'b',lwd = 3,col ="#f7903c" ,pch = 16)
# box(lwd = 1)
# abline(0,1,lty = 3,lwd = 3,col = "blue")

plot(cal3,lwd=3,lty=2,errbar.col="#4e84bd",xlim = c(0,1),ylim = c(0,1),add = T)
lines(cal3[,c('mean.predicted','KM')],type = 'b',lwd = 3,col ="#4e84bd" ,pch = 16)

# plot(cal4,lwd=3,lty=2,errbar.col="red",xlim = c(0,1),ylim = c(0,1),add = T)
# lines(cal4[ , c('mean.predicted','KM')], type = 'b', lwd = 3, col ="red" , pch = 16)

box(lwd = 1)
abline(0,1,lty = 3,lwd = 3,col = "black")
legend("topleft", #图例的位???
       legend = c("1-Years","3-Years","5-Years"), #图例文字
       col =c("#58aa5a","#f7903c","#4e84bd"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗???
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边???
dev.off()

tiff("calibrate.tiff", width = 11,  height = 8, units = "in", bg = "white",
     res = 500)
par(mar = c(10,5,3,2),cex = 1.0)
plot(cal1,lwd=3,lty=2,errbar.col="#58aa5a",xlim = c(0,1),ylim = c(0,1),
     xlab ="Nomogram-Predicted Probability of Survival", ylab="Actual Survival")
lines(cal1[,c('mean.predicted','KM')],type = 'b',lwd = 3,col ="#58aa5a" ,pch = 16)
# box(lwd = 1)
# abline(0,1,lty = 3,lwd = 3,col = "red")

plot(cal2,lwd=3,lty=2,errbar.col="#f7903c",xlim = c(0,1),ylim = c(0,1),add = T)
lines(cal2[,c('mean.predicted','KM')],type = 'b',lwd = 3,col ="#f7903c" ,pch = 16)
# box(lwd = 1)
# abline(0,1,lty = 3,lwd = 3,col = "blue")

plot(cal3,lwd=3,lty=2,errbar.col="#4e84bd",xlim = c(0,1),ylim = c(0,1),add = T)
lines(cal3[,c('mean.predicted','KM')],type = 'b',lwd = 3,col ="#4e84bd" ,pch = 16)

#plot(cal4,lwd=3,lty=2,errbar.col="red",xlim = c(0,1),ylim = c(0,1),add = T)
#lines(cal4[ , c('mean.predicted','KM')], type = 'b', lwd = 3, col ="red" , pch = 16)

box(lwd = 1)
abline(0,1,lty = 3,lwd = 3,col = "black")
legend("topleft", #图例的位???
       legend = c("1-Years","3-Years","5-Years"), #图例文字
       col =c("#58aa5a","#f7903c","#4e84bd"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗???
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边???
dev.off()


if(F){
  #???ɺ???
  rt$futime = rt$futime/365
  f <- cph(Surv(futime, fustat) ~ age+T+N+M+stage+riskScore,
           x=T, y=T, surv=T, data=rt, time.inc=1)
  surv <- Survival(f)
  nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(3, x),function(x) surv(5, x)),
                  lp=F, funlabel=c("1-year survival", "3-year survival", "5-year survival"),
                  maxscale=100,
                  fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))  
  
  #nomogram???ӻ?
  pdf(file="nomogram.pdf", height = 8, width = 9)
  plot(nom)
  dev.off()
  
}


###################################################################################
################################# cli_plot ########################################
###################################################################################
rm(list=ls())
setwd("E:/LZ/24105")  #modify work path
if(! dir.exists("12_nomogram")){
  dir.create("12_nomogram")
}
setwd("./12_nomogram")
data = read.table("./indepInput_clean.txt", 
                  header = T, check.names = F, sep = "\t", )
#cli = cli[ , -c(2, 3, 10)]
#risk = read.table("risk.txt", header = T, check.names = F, sep = "\t")
#data = merge(cli, risk, by.x = "id", by.y = "id")

library(ggpubr)
library(ggsci)

data$age = factor(data$age, levels = c(">=65", "<65"))
data$stage = factor(data$stage, levels = c("Stage I+II", "Stage III+IV"))
data$T = factor(data$T, levels = c("T1+T2", "T3+T4"))
data$gender = factor(data$gender, levels = c("female", "male"))
data$M = factor(data$M, levels = c("M0", "M1"))
data$N = factor(data$N, levels = c("N0", "N1"))

#-------age---------------------
if(T){
  data1 = data[ , c("age", "riskScore")] %>% na.omit()
  p <- ggviolin(data1, x = "age", y = "riskScore", fill = "age",
                palette = c("#4DBBD5FF", "#DC0000FF"), add = "boxplot",
                add.params = list(fill="white"))+
    ylab("riskScore")+
    stat_compare_means(label = "p.format",label.x = 1.5, size = 4)+
    theme(text = element_text(size = 15))
  p
  ggsave(filename = "age.tiff", plot = p, dpi = 500, width = 6, height = 5)
  ggsave(filename = "age.pdf", plot = p, dpi = 500, width = 6, height = 5)
}



#-------stage---------------------
data1 = data[ , c("stage", "riskScore")] %>% na.omit()
p <- ggviolin(data1, x = "stage", y = "riskScore", fill = "stage",
              palette = c("#4DBBD5FF", "#DC0000FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("riskScore")+
  stat_compare_means(label = "p.format",label.x = 1.5, size = 4)+
  theme(text = element_text(size = 15))
p
ggsave(filename = "stage.tiff", plot = p, dpi = 500, width = 6, height = 5)
ggsave(filename = "stage.pdf", plot = p, dpi = 500, width = 6, height = 5)


#------------------T----------------------
data1 = data[ , c("T", "riskScore")] %>% na.omit()
p <- ggviolin(data1, x = "T", y = "riskScore", fill = "T",
              palette = c("#4DBBD5FF", "#DC0000FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("riskScore")+
  stat_compare_means(label = "p.format",label.x = 1.5, size = 4)+
  theme(text = element_text(size = 15))
p
ggsave(filename = "T.tiff", plot = p, dpi = 500, width = 6, height = 5)
ggsave(filename = "T.pdf", plot = p, dpi = 500, width = 6, height = 5)


#------------------M----------------------
data1 = data[ , c("M", "riskScore")] %>% na.omit()
p <- ggviolin(data1, x = "M", y = "riskScore", fill = "M",
              palette = c("#4DBBD5FF", "#DC0000FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("riskScore")+
  stat_compare_means(label = "p.format",label.x = 1.5, size = 4)+
  theme(text = element_text(size = 15))
p
ggsave(filename = "M.tiff", plot = p, dpi = 500, width = 6, height = 5)
ggsave(filename = "M.pdf", plot = p, dpi = 500, width = 6, height = 5)


#------------------N----------------------
data1 = data[ , c("N", "riskScore")] %>% na.omit()
p <- ggviolin(data1, x = "N", y = "riskScore", fill = "N",
              palette = c("#4DBBD5FF", "#DC0000FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("riskScore")+
  stat_compare_means(label = "p.format",label.x = 1.5, size = 4)+
  theme(text = element_text(size = 15))
p
ggsave(filename = "N.tiff", plot = p, dpi = 500, width = 6, height = 5)
ggsave(filename = "N.pdf", plot = p, dpi = 500, width = 6, height = 5)


#------------------gender----------------------
data1 = data[ , c("gender", "riskScore")] %>% na.omit()
p <- ggviolin(data1, x = "gender", y = "riskScore", fill = "gender",
              palette = c("#4DBBD5FF", "#DC0000FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("riskScore")+
  stat_compare_means(label = "p.format",label.x = 1.5, size = 4)+
  theme(text = element_text(size = 15))
p
ggsave(filename = "gender.tiff", plot = p, dpi = 500, width = 6, height = 5)
ggsave(filename = "gender.pdf", plot = p, dpi = 500, width = 6, height = 5)




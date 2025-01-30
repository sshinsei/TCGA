#wgcna-------
rm(list = ls())
setwd('E:/LZ/24080')
if (! dir.exists('./03_WGCNA_GEO')){
  dir.create('./03_WGCNA_GEO')
}
setwd('./03_WGCNA_GEO')

# group <- read.csv('../00_rawdata/01.group.csv')
# colnames(group) <- c('sample', 'Type')
# group <- group[group$Type=='Tumor',]
# 
# dat_expr <- read.csv('../00_rawdata/01.fpkmlog.csv',row.names = 1,check.names = F) %>% lance::lc.tableToNum()
# expr <- dat_expr[,group$sample]

# 载入数据

eset <- read.csv(file ="../00_rawdata/01_GEO/GSE76427.txt",
                 sep="\t",row.names=1)
group <- read.csv(file = "../00_rawdata/01_GEO/GSE76427.group.txt",
                  sep="\t")
# 加载表达矩阵数据
# eset <- read.csv(file ="../00_rawdata/16combat/combat.exp.txt",
#                sep="\t",row.names = 1)
# 加载分组变量
# group <- read.csv(file = "../00_rawdata/16combat/group.txt",
#                  sep="\t")

# 转置eset
eset <- t(eset)

colnames(group) <- c('sample', 'group')
datExprOri=as.data.frame(eset)
nGenes0 = ncol(datExprOri)
nSamples0 = nrow(datExprOri)
nGenes0 ;nSamples0
# [1] 17491
# [1] 424
datExpr1 <- datExprOri

#去除低表达的基因80%
#https://blog.csdn.net/sakoko_/article/details/108907543
# datExpr1 <- t(datExpr1) %>% as.data.frame()
# countExpr<-data.frame(apply(datExpr1,1,function(x) {re=sum(x==0)}))
# colnames(countExpr)<-'count Expr==0'
# rmORnot<-data.frame(apply(countExpr,1,function(x) {if(x>ncol(datExpr1)*0.2){re=0}else{re=1}}))#1保留，0去掉
# colnames(rmORnot)<-'remove gene or not'
# datExpr1<-datExpr1[which(rmORnot==1),]
# datExpr1 <- t(datExpr1) %>% as.data.frame()
# nGenes1<-ncol(datExpr1)      #基因数目
# nSamples1 <- nrow(datExpr1) #样品名称变量
# nGenes1 ;nSamples1
# # [1] 19445
# # [1] 16



#筛选方差前75%的基因
#https://zhuanlan.zhihu.com/p/498684272
# datExpr1 <- t(datExpr1) %>% as.data.frame()
# m.mad <- apply(datExpr1,1,mad)
# datExpr1 <- datExpr1[which(m.mad >
#                              max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.001)),]
# datExpr1 <- t(datExpr1) %>% as.data.frame()
# nGenes1<-ncol(datExpr1)      #基因数目
# nSamples1 <- nrow(datExpr1) #样品名称变量
# nGenes1 ;nSamples1
# [1] 8827
# [1] 112


#过滤基因，取绝对中位差MAD top 5000 的基因 PMID：36818670
#https://zhuanlan.zhihu.com/p/401089111
# datExpr1 <- t(datExpr1) %>% as.data.frame()
# datExpr1 <- datExpr1[order(apply(datExpr1, 1, mad), decreasing = T)[1:5000],]
# datExpr1 <- t(datExpr1) %>% as.data.frame()
# nGenes1<-ncol(datExpr1)      #基因数目
# nSamples1 <- nrow(datExpr1) #样品名称变量
# nGenes1 ;nSamples1
# [1] 5000
# [1] 365

## #检测缺失值
library(WGCNA)
##检验数据质量
gsg = goodSamplesGenes(datExpr1 , verbose = 3)   #计算数据集中是否含有方差为0的基因
gsg$allOK   ##为TRUE不用删除基因
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(datExpr1)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
  
}      ##输出的datExpr1已经过滤基因
datExpr1<-datExpr1
nGenes1<-ncol(datExpr1)      #基因数目
nSamples1 <- nrow(datExpr1) #样品名称变量
nGenes1 ;nSamples1
# [1] 16739
# [1] 189




##根据层次聚类绘制样本
# 观察是否有离群样品需要剔除
tree=hclust(dist(datExpr1),method ='complete')   ##观察是否有离群样本  c("single", "complete", "average", "mcquitty", "ward.D", "centroid", "median", "ward.D2")
#hclust()函数是stats包中的函数,可以根据距离矩阵实现层次聚类。
pdf(file='01.sampleClustering.pdf',w=10,h=7)
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plot(tree,xlab="", sub="", main="Sample Clustering",
     labels=F,
     cex=1.0,  ##label大小
     font=2,
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2)
#abline(h =300, col = 'red')
dev.off()

png(file='01.sampleClustering.png',w=10,h=7,units='in',res=600,bg='white')
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plot(tree,xlab="", sub="", main="Sample Clustering",
     labels=F,
     cex=1.0,  ##label大小
     font=2,
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2)
#abline(h =300 ,col = 'red')
dev.off()

## 从上图可以看出，疑似异样，考虑去除:cutHeight=300
# clust = cutreeStatic(tree, cutHeight = 300, minSize = 100)
# table(clust)  ###无离群样本
## 去除异常样品后的表达矩阵
# datExpr = datExpr1[clust == 1, ]
# nGenes = ncol(datExpr)      #基因数目
# nSample =nrow(datExpr) #样品名称变量
# nGenes ;nSample
# SampleName<-rownames(datExpr)

if(F){
  # ## 表型数据聚类树
  # group <- group[group$sample %in% SampleName,]
  # condition <- group
  # rownames(condition) <- condition$sample
  # condition$Tumor<-ifelse(condition$Type=='Tumor',1,0)
  # condition$Normal<-ifelse(condition$Type=='Normal',1,0)
  # condition<-condition[,-c(1,2)]
  # 
  # # library(tidyverse)
  # datTraits<-condition
  # 
  # datExpr[1:4,1:4]
  # 
  # traitColors = numbers2colors(datTraits, signed = FALSE)
  # tree2=hclust(dist(datExpr),method ='complete')   ##可视化表型数据与基因表达量数据的联系，重构样本聚类树  c("single", "complete", "average", "mcquitty", "ward.D", "centroid", "median", "ward.D2")
  # pdf(file='01.sampleClustering_group.pdf',w=10,h=7)
  # par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
  # plotDendroAndColors(tree2,
  #                     traitColors,
  #                     #dendroLabels = T,
  #                     hang = 0.03,
  #                     cex.dendroLabels = 1,
  #                     #addGuide = TRUE,   ##添加网格线
  #                     font=2,
  #                     guideHang = 0.05,
  #                     cex.axis=1.4,  ##坐标轴刻度文字的缩放倍数。类似cex。
  #                     cex.lab=1.6,   ##坐标轴刻度文字的缩放倍数。类似cex。
  #                     cex.main=1.6,   ##标题的缩放倍数。
  #                     font.axis = 2,
  #                     font.lab = 2,
  #                     font.main = 2,
  #                     font.sub =2,
  #                     cex.colorLabels=1,
  #                     groupLabels = colnames(datTraits),
  #                     main = "Sample Clustering and trait heatmap")
  # dev.off()
  # 
  # png(file='01.sampleClustering_group.png',w=10,h=7,units='in',res=600,bg='white')
  # par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
  # plotDendroAndColors(tree2,
  #                     traitColors,
  #                     #dendroLabels = T,
  #                     hang = 0.03,
  #                     cex.dendroLabels = 1,
  #                     #addGuide = TRUE,   ##添加网格线
  #                     font=2,
  #                     guideHang = 0.05,
  #                     cex.axis=1.4,  ##坐标轴刻度文字的缩放倍数。类似cex。
  #                     cex.lab=1.6,   ##坐标轴刻度文字的缩放倍数。类似cex。
  #                     cex.main=1.6,   ##标题的缩放倍数。
  #                     font.axis = 2,
  #                     font.lab = 2,
  #                     font.main = 2,
  #                     font.sub =2,
  #                     cex.colorLabels=1,
  #                     groupLabels = colnames(datTraits),
  #                     main = "Sample Clustering and trait heatmap")
  # dev.off()
  
}

# 一步构建共表达网络和划分模块
##设置网络构建参数选择范围，计算beta值（为了计算power）
powers <- c(seq(1, 10, by=1), seq(12, 20, by=2)) #设置beta值的取值范围
#powers <- c(seq(1, 30, by=2)) #设置beta值的取值范围
#rsquare.cut=0.85  #设置r2阈值为0.85，一般选择0.85或0.9

#开启多线程，加快运行速度
enableWGCNAThreads()
#disbleWGCNAThreads()     #如果运行pickSoftThreshold报call()相关的错误，需要运行禁用多线程
#trace(pickSoftThreshold,edit = T) #148行，2改成4，下边画图改成sft$fitIndices[,4]
#untrace(pickSoftThreshold)
num <- 0.85
# sft = pickSoftThreshold(datExpr1, RsquaredCut = num, powerVector = powers, verbose = 5)
sft = pickSoftThreshold(datExpr1, RsquaredCut = 0.85, powerVector = powers, verbose = 5,nBreaks = 20)
sft$powerEstimate ###### 4
save.image('powerEstimate_0.85.RData')


load('./powerEstimate_0.85.RData')
#软阈值选择绘图
#左图为每个候选beta对应的R2值,
#绘制阈值线0.9（红线）和0.85（蓝线），如没有超过蓝线，则没找到软阈值，sft$powerEstimate为空
#右图为每个候选beta对应的平均连接度
pdf('02.softThreshold.pdf',w=12,h=8)
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1),family='serif')
cex1=1.5
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2,
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=num,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"),
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

png('02.softThreshold.png',w=12,h=8,units='in',res=600,bg='white')
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1),family='serif')
cex1=1.5
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2,
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=num,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"),
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，
# 数值越高，网络越符合无标度特征 (non-scale)。
##自动构建WGCNA模型
#表达矩阵转换成邻接矩阵，然后再将邻接矩阵转换成拓扑矩阵,识别模块

POWER = sft$powerEstimate


datExpr <- datExpr1
cor <- WGCNA::cor
net = blockwiseModules(
  datExpr,
  # power = sft$powerEstimate,
  power = POWER,
  minModuleSize =300,                    #每个模块最少的基因数
  deepSplit = 2,                         #剪切树参数，一般设置为2，取值0-4（最灵敏）
  mergeCutHeight = 0.25,                 #模块合并参数，越大模块越少，默认0.25
  numericLabels = TRUE,                  #T返回数字，F返回颜色
  networkType  = "signed",               #网络类型，一般不需要修改
  maxBlockSize = ncol(datExpr),
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "FPKM-TOM",
  loadTOMs = TRUE,
  verbose = 3
)

cor<-stats::cor
table(net$colors)
#   0    1    2    3    4    5    6    7    8    9   10 
# 6012 2610 1560 1534 1210 1187  869  790  752  500  467 


##模块重新排序
#获取模块颜色
moduleColors = labels2colors(net$colors)
#获取每个模块的特征值
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); #重新排序，相似颜色的挨在一起
useMEs=subset(MEs, select = -c(MEgrey))  #去掉未分组的基因(灰色模块）

# useMEs=subset(MEs)  #去掉未分组的基因(灰色模块）
# 输出每个基因所在的模块，以及与该模块的KME值
if(file.exists('01.All_Gene_KME.txt')) file.remove('01.All_Gene_KME.txt')
for(module in substring(colnames(useMEs),3)){
  ME=as.data.frame(useMEs[,paste("ME",module,sep="")])
  colnames(ME)=module
  datModExpr=datExpr1[,moduleColors==module]
  datKME = signedKME(datModExpr, ME)
  datKME= cbind(datKME,rep(module,length(datKME)))
  write.csv(datKME,quote = F,row.names = T,file = "01.All_Gene_KME.csv",col.names = F)
}

## 模块绘图
mergedColors = labels2colors(net$colors)
mergedColors[net$blockGenes[[1]]]
pdf("03.wgcna.dendroColors.pdf",height = 7,width = 9)
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1),family='serif')
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                    cex.axis=1.6,  ##坐 标轴刻度文字的缩放倍数。类似cex。
                    cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
                    cex.main=1.8,   ##标题的缩放倍数。
                    font.axis = 2,
                    font.lab = 2,
                    font.main = 2,
                    cex.colorLabels=1,
                    font.sub =2)

dev.off()

png("03.wgcna.dendroColors.png",height = 7,width = 9,units='in',res=600)
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1),family='serif')
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                    cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
                    cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
                    cex.main=1.8,   ##标题的缩放倍数。
                    font.axis = 2,
                    font.lab = 2,
                    font.main = 2,
                    cex.colorLabels=1,
                    font.sub =2)

dev.off()

save.image('Modules.RData')
# load('./Modules.RData')
if(F){
  # 模块间相关性图
  # module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
  MEs = net$MEs
  
  ### 不需要重新计算，改下列名字就好
  ### 官方教程是重新计算的，起始可以不用这么麻烦
  MEs_col = MEs
  colnames(MEs_col) = paste0("ME", labels2colors(
    as.numeric(str_replace_all(colnames(MEs),"ME",""))))
  MEs_col = orderMEs(MEs_col)
  
  # 根据基因间表达量进行聚类所得到的各模块间的相关性图
  # marDendro/marHeatmap 设置下、左、上、右的边距
  plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                        marDendro = c(3,3,2,4),
                        marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                        xLabelsAngle = 90)
  
}




#--------------wgcna-------
rm(list = ls())
setwd('E:/LZ/24080')
if (! dir.exists('./03_WGCNA_GEO')){
  dir.create('./03_WGCNA_GEO')
}
setwd('./03_WGCNA_GEO')
# 将结果保存成cytoscape的输入文件格式
# 读入构建时保存的TOM矩阵
load('FPKM-TOM-block.1.RData')
TOM=as.matrix(TOM)
load('./Modules.RData')
library(WGCNA)
library(tidyverse)
if(F){
  ## 模块与表型的关联分析
  # 表型数据：
  datTraits <- read.csv("/data/nas1/huyuchen/project/01_KMZK-40205-2/02_TGRS/02.CIBERSORT_tcell_result.csv",row.names = 1)
  # datTraits <- datTraits %>% 
  # rownames_to_column("sample")
  # save(file="WGCNA.RData",datExpr1,SampleName,useMEs)
  str(datTraits)
  # load('WGCNA.RData')
  end_col = ncol(datTraits)-1
  # moduleTraitCor = cor(useMEs, datTraits[,2:end_col], use = "p")
  # moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples0)
  
  #  使用PEARSON
  if (corType=="pearson") {
    modTraitCor = cor(useMEs, datTraits[,2:end_col], use = "p")
    modTraitP = corPvalueStudent(moduleTraitCor, nrow(datExpr1)) %>% signif(3)
    modTraitPadj = p.adjust(modTraitP, method = "BH")
  } else {
    modTraitCorP = bicorAndPvalue(useMEs, datTraits, robustY=robustY)
    modTraitCor = modTraitCorP$bicor
    modTraitP   = modTraitCorP$p
    modTraitPadj = p.adjust(modTraitP, method = "BH")
  }
}

###-----------------------模块与表型的关联分析：case/normal---------------------
if(F){
  ## 模块与表型的关联分析
  # 表型数据：case = 1
  datTraits <- read.csv(file = "../00_rawdata/group.txt",
                        sep="\t") %>% mutate(Normal = ifelse(group == "normal",1,0),
                                             Case = ifelse(group == "tumor",1,0)) %>% 
    column_to_rownames("sample") %>% select("Normal","Case")
  
  #datTraits <- eset %>% as.data.frame() %>% 
  #  rownames_to_column("sample") %>% 
  #  left_join(.,group,by="sample") %>% 
  #  column_to_rownames("sample")
  
  # datTraits <- datTraits %>% 
  # rownames_to_column("sample")
  # save(file="WGCNA.RData",datExpr1,SampleName,useMEs)
  str(datTraits)
  # load('WGCNA.RData')
  # end_col = ncol(datTraits)-1
  # moduleTraitCor = cor(useMEs, datTraits[,2:end_col], use = "p")
  # moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples0)
  
  #  使用PEARSON
  if(F){
    if (corType=="pearson") {
      modTraitCor = cor(useMEs, datTraits[,2:end_col], use = "p")
      modTraitP = corPvalueStudent(moduleTraitCor, nrow(datExpr1)) %>% signif(3)
      modTraitPadj = p.adjust(modTraitP, method = "BH")
    } else {
      modTraitCorP = bicorAndPvalue(useMEs, datTraits, robustY=robustY)
      modTraitCor = modTraitCorP$bicor
      modTraitP   = modTraitCorP$p
      modTraitPadj = p.adjust(modTraitP, method = "BH")
    }
  }
  
}


#----------------------------基于免疫浸润与表型的关联分析-------------------------------------------------------
if(F){
  datTraits <- read.csv("../02_ssGSEA/Score_result.csv",
                        row.names = 1)
  str(datTraits)
  
  # 读取显著性数据表格
  datSig <- read.csv("../02_ssGSEA/immucell_wilcox_test.csv")
  # 删除sig列的英文数字
  datSig$sig <- gsub('[A-Za-z0-9]',"",datSig$sig)
  # 删除sig列前后的多余空格
  datSig$sig <- gsub('^\\s+|\\s+$',"",datSig$sig)
  
  # 筛选有显著性的
  library(tidyverse)
  mask <- datSig %>% 
    filter(sig != "") %>% 
    pull(immune_cell)
  datTraits <- datTraits[rownames(datTraits) %in% mask,]
  
  # 计算每列的均值
  #mean_values <- colSums(datTraits)
  # 将均值转换为数据框的一行
  #mean_row <- as.data.frame(t(mean_values))
  # 将新的一行添加到原数据框中
  #datTraits <- rbind(datTraits, mean_row)
  # 添加行名"score"以标识新行
  #rownames(datTraits)[nrow(datTraits)] <- "score"
  # 只保留score
  #datTraits <- datTraits["score",]
  
  datTraits <- t(datTraits) 
}

# --------------------------ssGSEA评分作为表形性状----------------------------------
if(T){
  score <- read.csv('../02_ssGSEA/01.ssgsea_score_GEO.csv',check.names = F,row.names = 1)
  #score<-data.frame(score)%>%t(.)
  colnames(score)[1] <- c('NRGs')
  SampleName<-rownames(datExpr1)
  score <- score[SampleName,,drop=F]
  #score1 <- read.csv('../02_ssGSEA/01.ssgsea_score_MP-RGs.csv',check.names = F,row.names = 1)
  #score1 <- score1[,c('MP-RGs'),drop=F]
  #colnames(score1)[1] <- c('MP-RGs')
  #score1 <- score1[SampleName,,drop=F]
  
  ##判断两个样本是否完全一致
  ##identical(rownames(score),rownames(score1))
  datTraits<-data.frame(row.names=rownames(score),NRGs=score$NRGs)
  
  identical(rownames(datExpr),rownames(datTraits))
  datExpr[1:4,1:4]
  datTraits[1:4,]
}






# ----------------------分组------------------------------------------------------
moduleTraitCor = cor(useMEs, datTraits, use = "p") # use=p代表去掉缺失计算
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples0)#
#整理要显示在图中的数字,第一行为相关系数，第二行为p值
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

#绘制关联热图
pdf("04.wgcna.Module-trait.heatmap.pdf", width = 10, height =15)
par(mar = c(14, 12, 3, 1.5),family='serif')  ##下左上右   ##par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
labeledHeatmap(Matrix = moduleTraitCor,   #相关系数
               xLabels = colnames(datTraits), #x轴为表型
               yLabels = names(useMEs),#y轴为模块
               ySymbols = names(useMEs),
               colorLabels = FALSE,
               colors =  blueWhiteRed(50),
               textMatrix = textMatrix, #每个单元格的内容
               setStdMargins = FALSE,
               xLabelsAngle = 45,
               xLabelsAdj = 0.98,
               cex.text =1.1,cex.lab.y = 1.5,cex.legendLabel =1.5,cex.lab.x=1.5,
               zlim = c(-1,1),
               main = "Module-trait relationships",
               cex.lab = 1.2,
               font.lab.x = 2, font.lab.y = 2 )
dev.off()

png("04.wgcna.Module-trait.heatmap.png", width = 10, height =15,unit='in',res=600,bg='white')
par(mar = c(14, 12, 3, 1.5),family='serif')  ##下左上右   ##par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
labeledHeatmap(Matrix = moduleTraitCor,   #相关系数
               xLabels = colnames(datTraits), #x轴为表型
               yLabels = names(useMEs),#y轴为模块
               ySymbols = names(useMEs),
               colorLabels = FALSE,
               colors =  blueWhiteRed(50),
               textMatrix = textMatrix, #每个单元格的内容
               setStdMargins = FALSE,
               xLabelsAngle = 45,
               xLabelsAdj = 0.98,
               cex.text =1.1,cex.lab.y = 1.5,cex.legendLabel =1.5,cex.lab.x=1.5,
               zlim = c(-1,1),
               main = "Module-trait relationships",
               cex.lab = 1.2,
               font.lab.x = 2, font.lab.y = 2 )
dev.off()




if(T){
  ##基因共表达网络热图
  #dissTOM = 1-TOM
  kME=signedKME(datExpr, useMEs, outputColumnName = "kME", corFnc = "cor", corOptions = "use = 'p'")
  # write.table(kME,"kME.txt",quote = F,sep = "\t",col.names = T)
  # if (dim(datExpr1)[2]>=1500) nSelect=1500 else nSelect=dim(datExpr1)[2]
  # set.seed(1)   #For reproducibility, we set the random seed
  # select = sample(nGenes1, size = nSelect)
  # selectTOM = dissTOM[select, select]
  
  # There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
  # 重新画聚类图
  # selectTree = hclust(as.dist(selectTOM),method ='complete')
  # selectColors = moduleColors[select]
}


# save.image(file='temp.RData')

# rm(list=ls())
# setwd('E:/02-OA/03_WGCNA/')
# library("WGCNA")
# library(tidyverse)
# load('temp.RData')
# ----------------------GS-MM:筛选出关键模块后再继续------------------#
# 不同模块的基因显著性图
geneTraitSignificance  = as.data.frame(cor(datExpr1, datTraits, use = "p"))
write.table(geneTraitSignificance,"GS.txt",quote = F,sep = "\t",row.names = T)
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 70))
names(geneTraitSignificance) = paste("GS.", colnames(datTraits), sep="")
names(GSPvalue) = paste("GS.", colnames(sample), sep="")
modNames = substring(names(MEs), 3)

#
# 
# GS-MM分析
# modNames1<-c("turquoise","brown")
modNames1<-c("yellow")
# modNames1 <- c("salmon","magenta","turquoise")

column = match("salmon", modNames1);
moduleGenes = moduleColors=="salmon";
salmon_kme<-kME[moduleGenes, ]
gene_salmon <- rownames(salmon_kme) %>% as.data.frame()
dim(gene_salmon)
# 381        1


column = match("magenta", modNames1);
moduleGenes = moduleColors=="magenta";
magenta_kme<-kME[moduleGenes, ]
gene_magenta <- rownames(magenta_kme) %>% as.data.frame()
dim(gene_magenta)
# 498        1


column = match("turquoise", modNames1);
moduleGenes = moduleColors=="turquoise";
turquoise_kme<-kME[moduleGenes, ]
gene_turquoise <- rownames(turquoise_kme) %>% as.data.frame()
dim(gene_turquoise)
# 2610        1


column = match("brown", modNames1);
moduleGenes = moduleColors=="brown";
brown_kme<-kME[moduleGenes, ]
gene_brown <- rownames(brown_kme) %>% as.data.frame()
dim(gene_brown)
# 2610        1


column = match("greenyellow", modNames1);
moduleGenes = moduleColors=="greenyellow";
greenyellow_kme<-kME[moduleGenes, ]
gene_greenyellow <- rownames(greenyellow_kme) %>% as.data.frame()
dim(gene_greenyellow)
# 442        1



column = match("green", modNames1);
moduleGenes = moduleColors=="green";
green_kme<-kME[moduleGenes, ]
gene_green <- rownames(green_kme) %>% as.data.frame()
dim(gene_green)
# 702        1


column = match("yellow", modNames1);
moduleGenes = moduleColors=="yellow";
yellow_kme<-kME[moduleGenes, ]
gene_yellow <- rownames(yellow_kme) %>% as.data.frame()
dim(gene_yellow)
# 732    1

column = match("blue", modNames1);
moduleGenes = moduleColors=="blue";
blue_kme<-kME[moduleGenes, ]
gene_blue <- rownames(blue_kme) %>% as.data.frame()
dim(gene_blue)
# 1560    1

column = match("black", modNames1);
moduleGenes = moduleColors=="black";
black_kme<-kME[moduleGenes, ]
gene_black <- rownames(black_kme) %>% as.data.frame()
dim(gene_black)
# 790    1

# 绘制关键模块的GS-MM图
for (module in modNames1){
  if(module== "grey"){ next }
  column = match(module, modNames); # col number of interesting modules
  pdf(paste("06.GS_MM.", module, ".pdf", sep=""),height = 6,width = 7)
  moduleGenes = moduleColors==module;
  par(mfrow = c(1,1))
  verboseScatterplot(abs(kME[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance",
                     main = paste("Module membership vs gene significance"),
                     cex.main = 1.2, cex.lab = 1.2, pch=19,cex.axis = 1.2, 
                     col = module,font.lab=2)
  abline(v=0.8,lwd=2,col="red")
  abline(h=0.2,lwd=2,col="red")
  dev.off()
  
  png(paste("06.GS_MM.", module, ".png", sep=""),height = 6,width = 7,units='in',res=600)
  moduleGenes = moduleColors==module;
  par(mfrow = c(1,1))
  verboseScatterplot(abs(kME[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance",
                     main = paste("Module membership vs gene significance"),
                     cex.main = 1.2, cex.lab = 1.2, pch=19,cex.axis = 1.2, 
                     col = module,font.lab=2)
  abline(v=0.8,lwd=2,col="red")
  abline(h=0.2,lwd=2,col="red")
  dev.off()
}


if(T){
  # 获得black模块下的所有基因和关键基因
  column = match("black", modNames);
  moduleGenes = moduleColors=="black";
  black_keep = abs(kME[moduleGenes, column])>0.8&abs(geneTraitSignificance[moduleGenes, 1]) >0.2
  black_kme<-kME[moduleGenes, ]
  black_gene<-rownames(black_kme[black_keep,])
  write.table(rownames(black_kme),file = "07.black_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
  write.table(black_gene,file = "08.black_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)
}



if(F){
  # 获得magenta模块下的所有基因和关键基因
  column = match("salmon", modNames);
  moduleGenes = moduleColors=="salmon";
  salmon_keep = abs(kME[moduleGenes, column])>0.8&abs(geneTraitSignificance[moduleGenes, 1]) >0.2
  salmon_kme<-kME[moduleGenes, ]
  salmon_gene<-rownames(salmon_kme[salmon_keep,])
  write.table(rownames(salmon_kme),file = "07.salmon_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
  write.table(salmon_gene,file = "08.salmon_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)
}


if(F){
  # 获得magenta模块下的所有基因和关键基因
  column = match("magenta", modNames);
  moduleGenes = moduleColors=="magenta";
  magenta_keep = abs(kME[moduleGenes, column])>0.8&abs(geneTraitSignificance[moduleGenes, 1]) >0.2
  magenta_kme<-kME[moduleGenes, ]
  magenta_gene<-rownames(magenta_kme[magenta_keep,])
  write.table(rownames(magenta_kme),file = "07.magenta_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
  write.table(magenta_gene,file = "08.magenta_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)
}

if(T){
  # 获得turquoise模块下的所有基因和关键基因
  column = match("turquoise", modNames);
  moduleGenes = moduleColors=="turquoise";
  turquoise_keep = abs(kME[moduleGenes, column])>0.8&abs(geneTraitSignificance[moduleGenes, 1]) >0.2
  turquoise_kme<-kME[moduleGenes, ]
  turquoise_gene<-rownames(turquoise_kme[turquoise_keep,])
  write.table(rownames(turquoise_kme),file = "07.turquoise_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
  write.table(turquoise_gene,file = "08.turquoise_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)
}



if(F){
  # 获得magenta模块下的所有基因和关键基因
  column = match("greenyellow", modNames);
  moduleGenes = moduleColors=="greenyellow";
  greenyellow_keep = abs(kME[moduleGenes, column])>0.8&abs(geneTraitSignificance[moduleGenes, 1]) >0.2
  greenyellow_kme<-kME[moduleGenes, ]
  greenyellow_gene<-rownames(greenyellow_kme[greenyellow_keep,])
  write.table(rownames(greenyellow_kme),file = "07.greenyellow_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
  write.table(greenyellow_gene,file = "08.greenyellow_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)
}

if(F){
  # 获得turquoise模块下的所有基因和关键基因
  column = match("green", modNames);
  moduleGenes = moduleColors=="green";
  green_keep = abs(kME[moduleGenes, column])>0.8&abs(geneTraitSignificance[moduleGenes, 1]) >0.2
  green_kme<-kME[moduleGenes, ]
  green_gene<-rownames(green_kme[green_keep,])
  write.table(rownames(green_kme),file = "07.green_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
  write.table(green_gene,file = "08.green_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)
}


if(T){
  # 获得turquoise模块下的所有基因和关键基因
  column = match("yellow", modNames);
  moduleGenes = moduleColors=="yellow";
  yellow_keep = abs(kME[moduleGenes, column])>0.8&abs(geneTraitSignificance[moduleGenes, 1]) >0.2
  yellow_kme<-kME[moduleGenes, ]
  yellow_gene<-rownames(yellow_kme[yellow_keep,])
  write.table(rownames(yellow_kme),file = "07.yellow_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
  write.table(yellow_gene,file = "08.yellow_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)
}

if(T){
  # 获得turquoise模块下的所有基因和关键基因
  column = match("blue", modNames);
  moduleGenes = moduleColors=="blue";
  blue_keep = abs(kME[moduleGenes, column])>0.8&abs(geneTraitSignificance[moduleGenes, 1]) >0.2
  blue_kme<-kME[moduleGenes, ]
  blue_gene<-rownames(blue_kme[blue_keep,])
  write.table(rownames(blue_kme),file = "07.blue_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
  write.table(blue_gene,file = "08.blue_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)
}

############################ 合并每个模块的hub_gene ######################################

#hub_gene <- list(c(blue_gene, yellow_gene))
hub_gene <- yellow_gene
# hub_gene <- list(c(magenta_gene, salmon_gene, turquoise_gene))
#write.table(hub_gene[[1]],file = "09.hub_gene.txt",sep = "\t",quote = F,row.names = F,col.names=F)
write.table(hub_gene,file = "09.hub_gene.txt",sep = "\t",quote = F,row.names = F,col.names=F)
# 142



if(F){
  
  modNames1<-c("black",'yellow')
  
  
  #
  # for (module in modNames1){
  #   if(module== "grey"){ next }
  #   column = match(module, modNames); # col number of interesting modules
  #   png(paste("06.GS_MM.", module, ".png", sep=""),height = 6,width = 7,family='Times',units='in',res=600)
  #   #png(paste("06.GS_MM.", module, ".png", sep=""),height = 6,width = 7,family='Times',units='in',res=600)
  #   moduleGenes = moduleColors==module;
  #   par(mfrow = c(1,1))
  #   verboseScatterplot(abs(kME[moduleGenes, column]),
  #                      abs(geneTraitSignificance[moduleGenes, 1]),
  #                      xlab = paste("Module Membership in", module, "module"),
  #                      ylab = "Gene significance",
  #                      main = paste("Module membership vs gene significance
  #
  # "),cex.main = 1.2, cex.lab = 1.2, pch=19,cex.axis = 1.2, col = module,font.lab=2)
  #   abline(v=0.8,lwd=2,col="red")
  #   abline(h=0.2,lwd=2,col="red")
  #   dev.off()
  # }
  #
  
  # column = match("magenta", modNames);
  # moduleGenes = moduleColors=="magenta";
  # magenta_kme<-kME[moduleGenes, ]
  # gene_magenta <- rownames(magenta_kme) %>% as.data.frame()
  
  
  save.image('result.RData')
  
  load('/data/nas3/zenghuanbing/project/01_YQTJ-30103-3-par/04_WGCNA/result.RData')
  
  
  column = match("blue", modNames);
  moduleGenes = moduleColors=="blue";
  blue_kme<-kME[moduleGenes, ]
  gene_blue <- rownames(blue_kme) %>% as.data.frame()
  dim(gene_blue)
  #1778        1
  
  
  column = match("turquoise", modNames);
  moduleGenes = moduleColors=="turquoise";
  turquoise_kme<-kME[moduleGenes, ]
  gene_turquoise <- rownames(turquoise_kme) %>% as.data.frame()
  dim(gene_turquoise)
  # 4225
  
  # column = match("white", modNames);
  # moduleGenes = moduleColors=="white";
  # purple_kme<-kME[moduleGenes, ]
  # gene_purple <- rownames(purple_kme) %>% as.data.frame()
  # 
  # column = match("purple", modNames);
  # moduleGenes = moduleColors=="purple";
  # skyblue_kme<-kME[moduleGenes, ]
  # gene_skyblue <- rownames(skyblue_kme) %>% as.data.frame()
  # 
  # column = match("purple", modNames);
  # moduleGenes = moduleColors=="purple";
  # purple_kme<-kME[moduleGenes, ]
  # gene_skyblu <- rownames(purple_kme) %>% as.data.frame()
  
  
  modlegene <- rbind(gene_turquoise,gene_blue)
  
  # DEG <- read.csv('../01_DEG/01.DEG_sig_padjlogfc1.csv')
  # DEG <- DEG[,1] %>% as.data.frame()
  # 
  # DEG2 <- read.csv('../02_ssGSEA/01.DEG_sig_score_adjplogfc1.5.csv')
  # DEG2 <- DEG2[,1] %>% as.data.frame()
  # intersect <- data.frame(symbol = intersect(DEG$., gene_turquoise$.))
  # intersect <- data.frame(symbol = intersect(DEG$., gene_purple$.))
  #gene_set <- read.csv('../02_ssGSEA/01.RETINOIC_ACID_METABOLIC_PROCESS.csv',row.names = 1)
  #intersect <- data.frame(symbol = intersect(DEG$., gene_set$x))
  #intersect <- data.frame(symbol = intersect(DEG$.,df.mm$gene_symbol))
  #intersect
  #write.csv(gene_turquoise,file = "01.gene_yellow.csv",quote = F,row.names = F)
  modlegene <- data.frame(symbol=modlegene$.)
  write.csv(modlegene,file = "04.modlegene.csv",quote = F,row.names = F)
  
  
  #brown
  # column = match("black", modNames);
  # moduleGenes = moduleColors=="black";
  # brown_keep = abs(kME[moduleGenes, column])>0.8&abs(geneTraitSignificance[moduleGenes, 1]) >0.2
  # brown_kme<-kME[moduleGenes, ]
  # brown_gene<-rownames(brown_kme[brown_keep,])
  # write.table(rownames(brown_kme),file = "07.brown_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
  # write.table(brown_gene,file = "08.brown_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)
  #
  # ##blue
  # column = match("yellow", modNames);
  # moduleGenes = moduleColors=="yellow";
  # yellow_keep = abs(kME[moduleGenes, column])>0.8&abs(geneTraitSignificance[moduleGenes, 1]) >0.2
  # yellow_kme<-kME[moduleGenes, ]
  # yellow_gene<-rownames(yellow_kme[yellow_keep,])
  # write.table(rownames(yellow_kme),file = "09.yellow_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
  # write.table(yellow_gene,file = "10.yellow_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)
  #
}





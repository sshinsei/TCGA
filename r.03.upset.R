rm(list=ls())
setwd("E:/LZ/24080")  #modify work path
if(! dir.exists("03_venn")){
  dir.create("03_venn")
}
setwd("./03_venn")
library(UpSetR)
da1 <- read.table("../02_DEG/diffmRNAExp.txt", header = T, sep = "\t") # 142
# da2 <- read.table("../03_WGCNA/07.black_all.txt", header = F, sep = "\t") # 88
# da2 <- read.table("../03_WGCNA/07.turquoise_all.txt", header = F, sep = "\t") # 790
da2 <- read.table("../03_venn/LRGs.txt", header = F, sep = "\t") # 236
# da2 <- read.csv("../02_DEG/DEG_sig.csv",header = T)
#miRTarBase <- read.table("miRTarBase.txt", header = F, sep = "\t")
#miRDB <- read.table("miRDB.txt", header = F, sep = "\t")
a <- list(LIHC_DEG = da1$id, LRGs = da2$V1) 
inter_ID <- Reduce(intersect,a) # 263
paste0(inter_ID,collapse = ', ')


pdf("upset.pdf", onefile=F)
# upset(fromList(a), sets.bar.color=c("red", "blue", "black", "yellow", "pink"))
upset(fromList(a), sets.bar.color=c("red", "blue"),
      text.scale = c(2, 1.5, 1.5, 1.5, 2, 2))
dev.off()

tiff(file="upset.tiff", width =8, height =7, units ="in", bg="white",
     res= 600)
upset(fromList(a), sets.bar.color=c("red", "blue"),
      text.scale = c(2, 1.5, 1.5, 1.5, 2, 2))
dev.off()


write.table(as.data.frame(inter_ID),'feature_genes.txt',row.names = FALSE,col.names=FALSE, sep = '\t', quote = FALSE)


###############################################################################
###-------------------------------cor plot----------------------------------
rm(list=ls())
setwd("E:/LZ/24080")  #modify work path
if(! dir.exists("03_venn")){
  dir.create("03_venn")
}
setwd("./03_venn")
norExp = read.csv(file = "../02_DEG/normalizeExp.txt", sep = "\t", header = T,
                  check.names = F, row.names = 1)
genes = read.csv("../06_multicox/multiCox.csv", header = T)

rt <- norExp[unlist(sapply(paste("^", paste(genes$id, "$", sep = ""), sep = ""),
                           grep, rownames(norExp))), ]
# rt <- rt[ , -(2:14)]  #delete normal sample
rt = as.matrix(rt)
# rownames(rt) = rt[ , 1]
# exp = rt[ , 2:ncol(rt)]
exp = rt
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
#data=log2(data + 0.001)
data=t(data)
#colnames(data)[17] = "DUPD1"
corr=cor(data)
library(corrplot)
tiff(file="cor.tiff",width = 8, height = 8, units = "in", bg = "white",
     res = 500)
corrplot(corr, type = "upper", method = "ellipse", order = "AOE",
         tl.cex = 1.0, tl.col = "black",
         insig = 'blank', # 不显著留空
         col=colorRampPalette(c("blue", "white", "red"))(300))
dev.off()
#write
write.table(corr, file = paste0("corr",".xls"), sep = "\t",
            row.names = T, quote = F)

pdf(file = "cor.pdf", height = 7, width = 7)
corrplot(corr, type = "upper", method = "ellipse", order = "AOE",
         tl.cex = 1.0, tl.col = "black",
         insig = 'blank', # 不显著留空
         col=colorRampPalette(c("blue", "white", "red"))(300))
dev.off()










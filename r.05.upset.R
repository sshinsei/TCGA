rm(list=ls())
setwd("E:/LZ/24080")  #modify work path
if(! dir.exists("05_venn")){
  dir.create("05_venn")
}
setwd("./05_venn")
library(UpSetR)
# da1 <- read.table("../02_DEG/diffmRNAExp.txt", header = T, sep = "\t") # 4932
da2 <- read.csv("../04_clustDEG/DEG_sig.csv", header = T) # 6868
da1 <- read.table("../03_WGCNA/09.hub_gene.txt", header = F, sep = "\t") # 142
#miRTarBase <- read.table("miRTarBase.txt", header = F, sep = "\t")
#miRDB <- read.table("miRDB.txt", header = F, sep = "\t")
# , WGCNA = da3$V1
a <- list(WGCNA = da1$V1, CLUST_DEGs = da2$X)
inter_ID <- Reduce(intersect,a) # 140
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








rm(list=ls())
setwd("E:/LZ/24105") 
if (!dir.exists("./20_TIDE")) {dir.create("./20_TIDE")}
setwd("20_TIDE")

Expr <- read.csv("../00_rawdata/normalizeExp.txt",sep = "\t",row.names = 1)
# 删除正常样本
colnames(Expr) <- substring(colnames(Expr),1 , 16)
colnames(Expr) <- gsub("\\.", "-", colnames(Expr))
Expr <- Expr[,!grepl("11A",colnames(Expr))]
# 二者选其一
Expr <- t(apply(Expr, 1, function(x)x-(mean(x)))) 
# Expr <- scale(Expr)
write.table(Expr, file="./Expr_tide.txt",sep="\t",row.names=T,quote=F)

rm(list=ls())
setwd("E:/LZ/24105") 
if (!dir.exists("./20_TIDE")) {dir.create("./20_TIDE")}
setwd("20_TIDE")
library(ggpubr)
rt = read.table("../00_rawdata/risk.txt", header = T, sep = "\t", check.names = F)
TIDE = read.table("tide.csv", header = T, sep = ",", check.names = F)
# tide的样本名处理------------
TIDE$Patient <- gsub("\\.","-",TIDE$Patient)

TIDE$Patient <- gsub("11A.*$", "11A", TIDE$Patient)
TIDE$Patient <- gsub("01A.*$", "01A", TIDE$Patient)
TIDE$Patient <- gsub("01B.*$", "01B", TIDE$Patient)

# TIDE <- TIDE[!grepl("\\.1",TIDE$Patient), ]
# 删除正常样本：
TIDE <- TIDE[!grepl("11A",TIDE$Patient), ]

#rt = rt[order(rt$risk, decreasing = T), ]
#inte = intersect(rt$id, TIDE$Patient)
#rows = na.omit(match(inte, TIDE$Patient))
#rownames(rt) = rt$id
#rt1 = cbind(rt[inte, c(1,5)], TIDE[rows, ])
#unique(rownames(rt1) == rt1$Patient)  #check
# rt$id = substring(rt$id, 1, 15)
rt1 = merge.data.frame(rt, TIDE, by.x = "id", by.y = "Patient")
rt1 = rt1[order(rt1$risk, decreasing = T), ]
rt1$risk = factor(rt1$risk, levels = c("high", "low"))
p <- ggviolin(rt1, x = "risk", y = "TIDE", fill = "risk",
              palette = c("#DC0000FF", "#4DBBD5FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("TIDE")+
  stat_compare_means(label = "p.format",label.x = 1.4, size = 5)+
  theme(text = element_text(size = 20))
p
ggsave(filename = "TIDE.tiff", plot = p, dpi = 400, width = 9, height = 8)
ggsave(filename = "TIDE.pdf", plot = p, dpi = 500, width = 9, height = 8)



################################################################################
###################################     TCIA   #################################
################################################################################
rm(list=ls())
library(ggpubr)
setwd('E:/LZ/24105')
if (! dir.exists('./20_TIDE')){
  dir.create('./20_TIDE')
}
setwd('./20_TIDE')
df <- read.csv("../00_rawdata/TCIA-ClinicalData.tsv",sep="\t")
library(tidyverse)
df1 <- df[ ,grepl("ips",colnames(df))]
df1 <- cbind(df$barcode,df1)
colnames(df1) <- c("barcode","ips_ctla4_neg_pd1_neg",	
                   "ips_ctla4_neg_pd1_pos",	"ips_ctla4_pos_pd1_neg",	"ips_ctla4")
write.table(df1,"./TCIA.txt",sep="\t",col.names  = T,
            row.names = F)

rm(list=ls())

rt = read.table("../00_rawdata/risk.txt", header = T, sep = "\t", check.names = F)
tcia = read.table("TCIA.txt", header = T, sep = "\t", check.names = F)


#rt = rt[order(rt$risk), ]
#inte = intersect(rt$id, tcia$barcode)
#rows = na.omit(match(inte, tcia$barcode))
#rownames(rt) = rt$id
#rt1 = cbind(rt[inte, c(1,5)], tcia[rows, ])
#unique(rt1$id == rt1$barcode) #check
# 修改样本名 ：与tcia中保持一致
rt$id = substring(rt$id, 1, 12)

#--------------------ips_ctla4_neg_pd1_neG--------------------------
rt1 = merge.data.frame(rt, tcia, by.x = "id", by.y = "barcode")
rt1 = rt1[order(rt1$risk, decreasing = T), ]
rt1$risk = factor(rt1$risk, levels = c("high", "low"))
p <- ggviolin(rt1, x = "risk", y = "ips_ctla4_neg_pd1_neg", fill = "risk",
              palette = c("#DC0000FF", "#4DBBD5FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("ips_ctla4_neg_pd1_neg")+
  stat_compare_means(label = "p.format",label.x = 1.25, size = 5, label.y = 11)+
  theme(text = element_text(size = 20))
p
ggsave(filename = "ips_ctla4_neg_pd1_neg.tiff", plot = p, dpi = 400, width = 8, height = 7)
ggsave(filename = "ips_ctla4_neg_pd1_neg.pdf", plot = p, dpi = 500, width = 8, height = 8)

#-------------------ips_ctla4_neg_pd1_pos-------------------------
p <- ggviolin(rt1, x = "risk", y = "ips_ctla4_neg_pd1_pos", fill = "risk",
              palette = c("#DC0000FF", "#4DBBD5FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("ips_ctla4_neg_pd1_pos")+
  stat_compare_means(label = "p.format",label.x = 1.25, size = 5, label.y = 11)+
  theme(text = element_text(size = 20))
p
ggsave(filename = "ips_ctla4_neg_pd1_pos.tiff", plot = p, dpi = 400, width = 8, height = 7)
ggsave(filename = "ips_ctla4_neg_pd1_pos.pdf", plot = p, dpi = 500, width = 8, height = 8)

#------------------ips_ctla4_pos_pd1_neg-----------------------------------
p <- ggviolin(rt1, x = "risk", y = "ips_ctla4_pos_pd1_neg", fill = "risk",
              palette = c("#DC0000FF", "#4DBBD5FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("ips_ctla4_pos_pd1_neg")+
  stat_compare_means(label = "p.format",label.x = 1.25, size = 5, label.y = 11)+
  theme(text = element_text(size = 20))
p
ggsave(filename = "ips_ctla4_pos_pd1_neg.tiff", plot = p, dpi = 400, width = 8, height = 7)
ggsave(filename = "ips_ctla4_pos_pd1_neg.pdf", plot = p, dpi = 500, width = 8, height = 8)


#----------------------------ips_ctla4----------------------------------
p <- ggviolin(rt1, x = "risk", y = "ips_ctla4", fill = "risk",
              palette = c("#DC0000FF", "#4DBBD5FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("ips_ctla4")+
  stat_compare_means(label = "p.format",label.x = 1.25, size = 5, label.y = 11)+
  theme(text = element_text(size = 20))
p
ggsave(filename = "ips_ctla4.tiff", plot = p, dpi = 400, width = 8, height = 7)
ggsave(filename = "iips_ctla4.pdf", plot = p, dpi = 500, width = 8, height = 8)













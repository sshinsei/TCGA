rm(list=ls())
setwd("G:/LZ/24080") 
if (!dir.exists("./20_TIDE")) {dir.create("./20_TIDE")}
setwd("20_TIDE")

Expr <- read.csv("../02_DEG/normalizeExp.txt",sep = "\t",row.names = 1)
# 删除正常样本
colnames(Expr) <- substring(colnames(Expr),1 , 16)
colnames(Expr) <- gsub("\\.", "-", colnames(Expr))
Expr <- Expr[,!grepl("11A",colnames(Expr))]
# 二者选其一
Expr <- t(apply(Expr, 1, function(x)x-(mean(x)))) 
# Expr <- scale(Expr)
write.table(Expr, file="./Expr_tide.txt",sep="\t",row.names=T,quote=F)

rm(list=ls())
setwd("G:/LZ/24080") 
if (!dir.exists("./20_TIDE")) {dir.create("./20_TIDE")}
setwd("20_TIDE")
library(ggpubr)
rt = read.table("../00_rawdata/risk.txt", header = T, sep = "\t", check.names = F)
TIDE = read.table("TIDE.csv", header = T, sep = ",", check.names = F)
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

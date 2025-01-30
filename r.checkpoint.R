rm(list=ls())
setwd("G:/LZ/24080") 
if (!dir.exists("./20_TIDE")) {dir.create("./20_TIDE")}
setwd("20_TIDE")
imm_checkp = read.table("check_point_gene.txt",sep="\t",header=F,check.names=F)
rt = read.csv(file = "../02_DEG/normalizeExp.txt", sep="\t", 
              header=T, check.names=F, row.names = 1)


#rt = rt[ , -c(1:51)]  #delete normal sample

# 去掉正常样本：
colnames(rt) <- substring(colnames(rt),1,16)
rt <- rt[,!grepl("11A",colnames(rt)[2:ncol(rt)])]

rows = na.omit(match(imm_checkp$V1, rownames(rt)))
rt_imm_checkp = rt[rows, ]

write.table(rt_imm_checkp, file="Immune_checkpoint_exp.txt", sep="\t",
            quote=F, col.names=T, row.names=T)


library(ggpubr)
library(tidyverse)

rt=rt_imm_checkp
#colnames(rt) = substring(colnames(rt), 1, 12)
Type=read.table("../18_estimate/cluster.txt", 
                sep = "\t", check.names = F, row.names = 1,
                header = F)
#Type=Type[order(Type[ , 2]), ]
# 处理rt的样本名
# getImmune中已经处理过样本，所以不用处理
if(F){
  # ----------------rt的列名-------------------------- #
  rt <- as.data.frame(t(rt))
  rt <- rt %>% 
    rownames_to_column("id")
  
  rt$id <- gsub("11A.*$", "11A", rt$id)
  rt$id <- gsub("01A.*$", "01A", rt$id)
  rt$id <- gsub("01B.*$", "01B", rt$id)
  
  # 去重
  rt <- rt %>% 
    distinct(id, .keep_all = TRUE) %>% 
    column_to_rownames("id")
  
  # 判断列名是否都唯一
  !any(duplicated(colnames(rt)))
  
  # 转置
  rt <- as.data.frame(t(rt))
  # --------------------------------------------- #
}


# 先取type的子集
Type <- Type[rownames(Type) %in% colnames(rt),]

cols = na.omit(match(rownames(Type), colnames(rt)))
rt=t(rt[ , cols])
unique(rownames(rt) == rownames(Type)) #check

#׼??????ͼ???????ļ?
data = data.frame()
rows = na.omit(match(rownames(rt), rownames(Type)))
for(i in colnames(rt)){
  a = as.data.frame(cbind(expression=rt[ , i], gene=i, Subtype=as.vector(Type[rows, 2])))
  a[ , 1] = as.numeric(a[ , 1])
  test <- wilcox.test(expression ~ Subtype, data = a)
  pValue = test$p.value
  if(!is.na(pValue) & pValue < 0.01){
    data = rbind(data, a)
    print(pValue)
  }
  #data=rbind(data,cbind(expression=rt[ , i], gene=i, Subtype=as.vector(Type[,2])))
}
write.table(data, file="data.txt",sep="\t",row.names=F,quote=F)

#????????ͼ
data=read.table("data.txt",sep="\t",header=T,check.names=F)  
data$expression = log2(data$expression + 0.01)
data$Subtype = factor(data$Subtype, levels = c("risk_low", "risk_high"))
p = ggboxplot(data, x="gene", y="expression", fill = "Subtype",
            ylab = "Gene expression",
            xlab = "", palette = c("#4DBBD5FF", "#DC0000FF"))+
  theme(text = element_text(size = 21))
p = p+rotate_x_text(60)
p = p+stat_compare_means(aes(group=Subtype),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")
p
ggsave(filename = "boxplot.tiff", plot = p, dpi = 400, width = 17, height = 8)
ggsave(filename = "boxplot.pdf", plot = p, dpi = 500, width = 17, height = 8)


rm(list=ls())
library(limma)
library(GEOquery)
library(tidyverse)


setwd("E:/LZ/24080/00_rawdata/01_GEO")

GEO_data <- 'GSE10186'
gene_annotation <- 'GPL5474'
gset <- getGEO(GEO_data,
               destdir = '.',
               GSEMatrix = T,
               getGPL = T) 

expr <- as.data.frame(exprs(gset[[1]]))
gpl <- getGEO(gene_annotation, destdir = '.') 
gpl <- Table(gpl)

# 提取基因symbol
colnames(gpl)

# 将GENE SYMBOL中以///连接的symbol进行处理
probe2symbol <- gpl %>%
  dplyr::select('ID','Symbol')%>%
  filter('Symbol'!='')%>%
  separate('Symbol',c('symbol','drop'),sep = '//')%>%
  dplyr::select(-drop)
probe2symbol=probe2symbol[probe2symbol$symbol!='' & probe2symbol$symbol!='---',]


names(probe2symbol) <- c('ID', 'symbol')

dat <- expr
dat$ID <- rownames(dat)
dat$ID <- as.character(dat$ID)


probe2symbol$ID <- as.character(probe2symbol$ID)

dat <- dat %>%
  inner_join(probe2symbol, by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol, everything())%>%     ## 重新排列
  mutate(rowMean = rowMeans(.[grep('GSM', names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol, .keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%  ## 反向选择去除rowMean这一列
  # filter(symbol != "---" & grepl("NM_",symbol)) %>% ##去除比对后基因名为---的行
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除


# 删除表达量全为0的基因
dat=dat[rowMeans(dat)>0,]


#--------------------------------group----------------------------------------------

# 样本信息
a <- gset[[1]]
pd <- pData(a)

# 样本分类，样本号对应分类
group <- data.frame(sample = pd$geo_accession, group = pd$title)
str(group$group)
group <- group %>% 
  mutate(group = case_when(# grepl("Non-Tumor",group) ~ "normal",
                           grepl("hcc",group) ~ "tumor",
                           TRUE ~ NA))
group <- group[,1:2]
group <- na.omit(group)
table(group$group)
# normal   tumor 
#  0      118
write.table(group,file=paste0(GEO_data,".group.txt"),sep="\t",quote=F,col.names=T,row.names = F)
dat <- dat[, group$sample]
write.table(dat,file=paste0(GEO_data,".txt"),sep="\t",quote=F,col.names=T)



#--------------------------生存信息---------------------------------------------#
###------------提取临床数据---------------
# 样本信息
a <- gset[[1]]
pd <- pData(a)
pdata1 <- pd
# 剔除正常样本
# 增加一列分组
pdata1$group <- ifelse(grepl("hcc", pdata1$title), "tumor", "normal")
pdata1 <- pdata1[pdata1$group == "tumor",]
pdata <- pdata1[,c("characteristics_ch1.1","characteristics_ch1.2")]
pdata$characteristics_ch1.1 <- gsub("\\D", "", 
                                    pdata$characteristics_ch1.1)
pdata$characteristics_ch1.2 <- gsub("\\D", "", 
                                    pdata$characteristics_ch1.2)
pdata$characteristics_ch1.2 <- substring(pdata$characteristics_ch1.2,3,3)

colnames(pdata) <- c("fustate","futimes")
pdata <- pdata %>% filter(fustate != "NA" & futimes != "NA")
# 确保所有列都为数值型
pdata[] <- lapply(pdata, function(x) as.numeric(as.character(x)))
pdata <- pdata %>% 
  rownames_to_column("sample")
write.table(pdata,file=paste0(GEO_data,".clinic.txt"),sep="\t",quote=F,col.names=T)



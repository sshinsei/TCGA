rm(list = ls())
setwd('E:/LZ/24105')
if (! dir.exists("./03_ssGSEA")){dir.create("./03_ssGSEA")}
setwd("03_ssGSEA")
#2.GSVA#####
library(GSVA)
library(ggpubr)
library(ggplot2)
library(devtools)
library(reshape2)


expr_region <- read.csv(file ="../00_rawdata/normalizeExp.txt",
                        sep="\t",row.names = 1)
group <- read.csv(file = "../00_rawdata/group.txt",
                  sep="\t")
colnames(expr_region)<-gsub('.','-',colnames(expr_region),fixed = T)

group$group <- factor(group$group,levels = c("normal", "tumor"))

case <- group[group$group == "tumor",]$sample
expr_region <- expr_region[,colnames(expr_region) %in% case]

condition_region <- group[group$group == "tumor",]
# 将数据框重新命名行名
#rownames(expr_region) <- expr_region[,1]
#expr_region <- expr_region[ , -1] # 删除第1列
# 确保所有列都为数值型
# expr_region[] <- lapply(expr_region, function(x) as.numeric(as.character(x)))
##转置
#expr_region <- t(expr_region)

#inner_gene1
gene <- read.table("../00_rawdata/genes.txt", header = F, sep ="\t")

###Cuproptosis score ####
gene.set.list<-data.frame(gene$V1)#


ssgsea_par <- ssgseaParam(as.matrix(expr_region), gene.set.list)  # all other values are default values
ssgsea_res <- gsva(ssgsea_par)
# ERROR:Calling gsva(expr=., gset.idx.list=., method=., ...) is defunct; use a method-specific parameter object (see '?gsva').
if(F){
  ssgsea_res <-
    gsva(
      as.matrix(expr_region),
      gene.set.list,
      method = "gsva",#“gsva”, “ssgsea”, “zscore”, “plage”
      kcdf = "Gaussian",
      abs.ranking = T
    )
}

ssgsea_score<-data.frame(ssgsea_res)%>%t(.)
# write.csv(ssgsea_score,'01.ssgsea_score.csv',quote=F)
med <- median(ssgsea_score)
ssgsea_score1 <- as.data.frame(ssgsea_score)
ssgsea_score1 <- ssgsea_score1 %>% rownames_to_column("sample")
colnames(ssgsea_score1) <- c("sample","score")
ssgsea_score1$sample <- gsub("\\.","-",ssgsea_score1$sample)

ssgsea_score1$clust <- ifelse(ssgsea_score1$score > med, "high","low")
clust <- ssgsea_score1[,c(1,3)]
colnames(clust)[2] <- "group"
write.table(clust,"cluster.txt",quote = F,row.names = F,col.names = T,sep="\t")


##箱线图
group<-data.frame(sample=condition_region$sample,group=clust$group)  
dat <- ssgsea_score  %>% as.data.frame %>% tibble::rownames_to_column(var = "sample")
dat$sample <- gsub("\\.","-",dat$sample)
identical(dat$sample,group$sample)
dat <- merge(dat, group, by = "sample")
dat2 <- tidyr::gather(dat, Genes, Score, -c(sample, group))

##p值检验
library(rstatix)
stat_res <- dat2 %>%
  group_by(Genes) %>%
  wilcox_test(Score ~ group) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p.adj")
stat_res
write.table(stat_res, file = "02.ssgsea_wilcoxon_res.xls", sep = "\t", row.names = F, quote = F)

##箱线图
clust$group <- as.character(clust$group)
boxplot_dat <-data.frame(cbind(ssgsea_score,clust$group))

colnames(boxplot_dat)[2]<-'group'
head(boxplot_dat)
boxplot_dat$id <- rownames(boxplot_dat)
boxplot_dat <- melt(boxplot_dat,id = c("id","group"))
boxplot_dat$value <-as.numeric(boxplot_dat$value)
boxplot_dat$group<- factor(boxplot_dat$group)
colnames(boxplot_dat)<-c('id','group','Gene','value')

library(ggsci)
library(ggplot2)

library(ggpubr)

p <- ggboxplot(
  boxplot_dat,
  x = "group",
  y = "value",
  fill = "group", 
  palette = c('#e65100','#1b5e20'),
  add = "jitter", # 添加数据点
  add.params = list(alpha = 0.7) # 设置数据点的透明度
) +
  stat_pvalue_manual(
    stat_res,
    y.position = max(boxplot_dat$value)+0.5, # 设置p值标注的位置
    size = 5,
    color = 'black',
    family = "Times",
    label = 'p.adj.signif',
    face = "bold"
  ) +
  theme_bw() + 
  xlab("") + 
  ylab("ssGSEA score") +
  theme(
    axis.title.x = element_text(size = 20, family = "Times", face = "bold"),
    axis.text.x = element_text(size = 16, family = "Times", face = "bold"),
    axis.title.y = element_text(size = 20, family = "Times", face = "bold"),
    axis.text.y = element_text(size = 16, family = "Times", face = "bold"),
    strip.text = element_text(size = 14, family = "Times", face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 15), 
    legend.text = element_text(size = 14)
  ) +
  # 增加箱型图边框
  theme(
    panel.border = element_rect(color = "black", size = 1.5)
  ) +
  # 修改图例的标题和文本样式
  guides(fill = guide_legend(title = "Group", title.position = "top", title.hjust = 0.5)) +
  theme(
    legend.title = element_text(size = 16, family = "Times", face = "bold"),
    legend.text = element_text(size = 14, family = "Times"),
    legend.background = element_rect(fill = "white", color = "black")
  )
p
ggsave(filename = "03.ssgsea_score_boxplot.pdf", width = 5.5, height = 5,p)
ggsave(filename = "03.ssgsea_score_boxplot.png", width = 5.5, height = 5,p)









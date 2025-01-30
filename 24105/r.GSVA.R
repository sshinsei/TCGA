rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)  #install.packages("msigdbr")
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)

setwd("E:/LZ/24105")
setwd("14_high_low_enrich")
# load("gsva_kegg.rdata")

## msigdbr包提取下载 一般尝试KEGG和GO做GSVA分析
##KEGG
KEGG_df_all <-  msigdbr(species = "Homo sapiens", # Homo sapiens or Mus musculus
                        category = "C2",
                        subcategory = "CP:KEGG") 
KEGG_df <- dplyr::select(KEGG_df_all,gs_name,gs_exact_source,gene_symbol)
kegg_list <- split(KEGG_df$gene_symbol, KEGG_df$gs_name) ##按照gs_name给gene_symbol分组


####  GSVA  ####
#GSVA算法需要处理logCPM, logRPKM,logTPM数据或counts数据的矩阵####
dat <- read.csv("E:/LZ/24080/01_TCGA/files/mRNA_count.txt", sep="\t", header=T, 
                check.names=F, row.names=1) 
# dat <-  read.csv("../00_rawdata/normalizeExp.txt", sep = "\t", row.names=1)
colnames(dat) <- substring(colnames(dat), 1, 16)
colnames(dat) <- gsub("\\.", "-", colnames(dat))
dat <- dat[,!grepl("11A",colnames(dat))]

#dat <- as.matrix(log2(edgeR::cpm(counts))+1)
#dat <- as.matrix(log2(tpm+1))
group <- read.csv(file = "../00_rawdata/risk.txt",
                  sep="\t")
group <- group[,c(1, 5)]
colnames(group) <- c("sample","group")

# subset
dat <- dat[,colnames(dat) %in% group$sample]
dat <- as.matrix(dat)
# geneset <- read.gmt("E:/02-OA/00_rawdata/00_temp/c2.cp.kegg.v7.5.symbols.gmt") 
geneset <- kegg_list

GSVAPARAM_C2 <- gsvaParam(exprData = dat, 
                                geneSets = geneset,
                                kcdf="Poisson" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                                minSize = 2,
                                maxSize = 500,
)
gsva_C2<- gsva(GSVAPARAM_C2)

write.csv(gsva_C2,"gsva_kegg_matrix.csv")
save.image("gsva_kegg.rdata")

#### 进行limma差异处理 ####
##设定 实验组exp / 对照组ctr
gl
exp="high"
ctr="low"
group_list <- group$group
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(gsva_C2)
contrast.matrix <- makeContrasts(contrasts=paste0(exp,'-',ctr),  #"high相对于low上下调
                                 levels = design)
contrast.matrix

fit1 <- lmFit(gsva_C2,design)                 #拟合模型
fit2 <- contrasts.fit(fit1, contrast.matrix) #统计检验
efit <- eBayes(fit2)                         #修正

summary(decideTests(efit,lfc=1, p.value=0.05)) #统计查看差异结果
tempOutput <- topTable(efit, coef=paste0(exp,'-',ctr), n=Inf)
degs <- na.omit(tempOutput) 
write.csv(degs,"gsva_kegg_degs.results.csv")


#### 发散条形图绘制 ####
library(tidyverse)  # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(ggthemes)
library(ggprism)
p_cutoff=0.05
# degs <- gsva_kegg_degs  #载入gsva的差异分析结果
Diff <- rbind(subset(degs,logFC>0)[1:20,], subset(degs,logFC<0)[1:20,]) #选择上下调前20通路     
dat_plot <- data.frame(id  = row.names(Diff),
                       p   = Diff$P.Value,
                       lgfc= Diff$logFC)
dat_plot$group <- ifelse(dat_plot$lgfc>0 ,1,-1)    # 将上调设为组1，下调设为组-1
dat_plot$lg_p <- -log10(dat_plot$p)*dat_plot$group # 将上调-log10p设置为正，下调-log10p设置为负
dat_plot$id <- str_replace(dat_plot$id, "KEGG_","");dat_plot$id[1:10]
dat_plot$threshold <- factor(ifelse(abs(dat_plot$p) <= p_cutoff,
                                    ifelse(dat_plot$lgfc >0 ,'Up','Down'),'Not'),
                             levels=c('Up','Down','Not'))

dat_plot <- dat_plot %>% arrange(lg_p)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)

## 设置不同标签数量
low1 <- dat_plot %>% filter(lg_p < log10(p_cutoff)) %>% nrow()
low0 <- dat_plot %>% filter(lg_p < 0) %>% nrow()
high0 <- dat_plot %>% filter(lg_p < -log10(p_cutoff)) %>% nrow()
high1 <- nrow(dat_plot)

p <- ggplot(data = dat_plot,aes(x = id, y = lg_p, 
                                fill = threshold)) +
  geom_col()+
  coord_flip() + 
  scale_fill_manual(values = c('Up'= '#36638a','Not'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-log10(p_cutoff),log10(p_cutoff)),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('-log10(P.Value) of GSVA score') + 
  guides(fill="none")+
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'black') + #黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 黑色标签
p
ggsave("GSVA_barplot_pvalue.pdf",p,width = 15,height  = 15)
ggsave("GSVA_barplot_pvalue.png",p,width = 15,height  = 15)




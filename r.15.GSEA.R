# GSEA ---------------------------------------------------------------------
rm(list = ls())
setwd('E:/LZ/24080')
if (!dir.exists("15_GSEA")) {dir.create("15_GSEA")}
setwd("15_GSEA")
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(limma)
library(tidyverse)
library(ggplot2)


# kegggmt<-read.gmt("E:/02-OA/00_rawdata/00_temp/c2.cp.v2023.2.Hs.symbols.gmt") #读gmt文件
kegggmt<-read.gmt("E:/02-OA/00_rawdata/00_temp/c2.cp.kegg.v7.5.symbols.gmt") #读gmt文件
# 分组信息
group <- read.csv(file = "../00_rawdata/group.txt",
                  sep="\t")
# 表达矩阵
expr <- read.csv(file ="../02_DEG/normalizeExp.txt",
                sep="\t",row.names=1)
# 关键基因
hub_gene <- read.csv('../06_multicox/multiCox.csv', header = T)
hub_gene <- hub_gene$id

hub.exp <- expr[rownames(expr) %in% hub_gene,] 
# hub.exp<-hub.exp[order(rownames(hub.exp),decreasing = F),]  ##按照字母进行排序

library(patchwork)


gsea.plot = function(res.kegg, top.hall, gene){
  # 获取选定路径的p.adjust值
  selected_pathways <- res.kegg@result[res.kegg@result$ID %in% top.hall, ]
  p_values <- paste0("p.adjust = ", formatC(selected_pathways$p.adjust, format = "e", digits = 2))
  
  # 构建数据
  gsdata <- do.call(rbind, lapply(top.hall, enrichplot:::gsInfo, object = res.kegg))
  gsdata$Description <- factor(gsdata$Description, levels = top.hall)
  
  p1 = ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + 
    # theme_bw() +
    theme_classic(14) +
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_brewer(palette = "Paired") +
    ggtitle(paste0("Gene Set Enrichment Analysis:",gene)) +
    # geom_hline(yintercept = 0, color = "black", size = 0.8) +
    geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) +
    theme(legend.position = "top", 
          legend.justification = c(0,1),
          legend.title = element_blank(), 
          legend.background = element_rect(fill = "transparent")) +
    guides(
      color = guide_legend(
        ncol = 1, 
        byrow = TRUE,
        reverse = T)
    )+
    ylab("Running Enrichment Score") + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          text = element_text(face = "bold", family = "Times"))+
    theme(
      axis.title.y =element_text(size=15,family = "Times", face = "bold"),
      axis.text.y=element_text(size=15,family = "Times", face = "bold"),
      panel.grid.major=element_blank(),panel.grid.minor=element_blank()
    )
  
  # 添加p值到图中
  for (i in seq_along(top.hall)) {
    p1 <- p1 + annotate(
      "text", 
      x = max(gsdata$x) * 0.75, # 调整x坐标位置
      y = max(gsdata$runningScore) + 0.15 * i, # 调整y坐标位置
      label = p_values[i], 
      size = 5, 
      hjust = 0,
      color = "black"
    )
  }
  
  i = 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1 }
  p2 = ggplot(gsdata, aes_(x = ~x)) + 
    geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) + 
    xlab(NULL) + ylab(NULL) + 
    # theme_bw() +
    theme_classic(14) +
    theme(legend.position = "none", 
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank(),
          text = element_text(face = "bold", family = "Times")) + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    scale_color_brewer(palette = "Paired")+
    theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      # legend.title=element_text(size=15) , 
      legend.text=element_text(size=14))
  p3 <- gseaplot2(res.kegg,top.hall,subplots = 3)+
    
    theme(axis.title.x =element_text(size=15,family = "Times", face = "bold"),
          axis.text.x =element_text(angle=0,size=15,hjust = 0,family = "Times", face = "bold"),
          axis.title.y =element_text(size=15,family = "Times", face = "bold"),
          axis.text.y=element_text(size=15,family = "Times", face = "bold"))+
    theme(
      panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
      legend.title=element_text(size=15) , legend.text=element_text(size=14))
  # p = aplot::insert_bottom(p1, p2, height = 0.15)
  p <- p1/p2/p3+plot_layout(ncol = 1, height = c(4,1,2))
  return(p)
}



dat <- expr
#
num <- nrow(hub.exp)
#每个基因富集一次
for (i in c(1:num)){
  # i <- 11
  train<-t(hub.exp)%>%as.data.frame()
  ## 以hub基因将样本进行高低表达量分析
  group=as.vector(ifelse(train[,i]>median(train[,i]),'high','low'))  ##根据TNFRSF1A表达量的大小将样本分组
  group <- factor(group,levels = c("high","low")) 
  condition<-data.frame(group)
  rownames(condition)<-rownames(train)
  design<-model.matrix(~0+group)   
  contrast.matrix<-makeContrasts(contrasts = "grouphigh-grouplow", levels = design) 
  fit <- lmFit(dat,design)     
  fit1 <- contrasts.fit(fit, contrast.matrix)    
  fit2 <- eBayes(fit1)  
  tempOutput <- topTable(fit2, coef=1, n=nrow(fit2),adjust="fdr")  ###所有基因检验结果
  genelist<-data.frame(SYMBOL=rownames(tempOutput ),logFC=tempOutput $logFC)
  #开始ID转换
  gene=bitr(genelist$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  ## 去重
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  gene_df <- data.frame(ENTREZID=gene$ENTREZID ,#可以是foldchange
                        SYMBOL = gene$SYMBOL) #记住你的基因表头名字
  gene_df <- merge(gene_df,genelist,by="SYMBOL")
  geneList<-gene_df $logFC #第二列可以是folodchange，也可以是logFC
  names(geneList)=gene_df$SYMBOL #使用转换好的ID
  geneList=sort(geneList,decreasing = T) #从高到低排序
  
  KEGG<-GSEA(geneList,TERM2GENE = kegggmt,pvalueCutoff = 0.05) #GSEA分析
  sortKEGG <-data.frame(KEGG)
  sortKEGG <- subset(sortKEGG,sortKEGG$qvalue<0.25 & abs(sortKEGG$NES) >1)
  sortKEGG<-sortKEGG[order(sortKEGG$p.adjust, decreasing = F),]#按照p.adjust从高到低排序
  write.table(sortKEGG,file=paste('0',i,'.',rownames(hub.exp)[i],'.xls',sep=''),sep="\t",quote=F,row.names=F)
  paths <- rownames(sortKEGG[c(1:5),])#选取你需要展示的通路ID
  
  #特定通路作图
  #trace(gseaplot2,edit=T)
  #p = gseaplot2(KEGG,c(1:5),base_size =10,color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),rel_heights = c(1.5, 0.3, 0.5),title=rownames(hub.exp)[i])
  p <- gsea.plot(KEGG, paths,rownames(hub.exp)[i])
  fn1 = paste0("0", i, ".", rownames(hub.exp)[i], ".png")
  fn2 = paste0("0", i, ".", rownames(hub.exp)[i], ".pdf")
  png(fn1,w=9,h=8,units = "in",res = 600)
  print(p)
  dev.off()
  pdf(fn2,w=9,h=8)
  print(p)
  dev.off()
  i<-i+1
}




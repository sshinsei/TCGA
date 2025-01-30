rm(list = ls())
setwd("E:/LZ/24080") 
if (!dir.exists("./13_immune/cibersort")) {dir.create("./13_immune/cibersort")}
setwd("13_immune")

library(GSVA)
library(magrittr)
library(IOBR)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(ggcorrplot)
library(ggthemes)
library(tidyverse)

# mark1 char set ----------------------------------------------------------
ssgsea = 0#选择ssGSEA算法 or CIBERSORT
xcell = 0
ciber= 1
epic = 0
datType = F#CIBERSORT算法中是测序数据还是芯片数据,QN芯片为T(no_log)，测序为F(FPKM) 
bar.cell = 1#是否绘制柱状堆积图
box.pic = 1#是否绘制箱线图展示差异结果
cor.map = 1#是否绘制基因与显著细胞的相关性热图。
singlecor = 1#是否绘制风险评分与显著细胞的相关性散点图
cell.cor = 1#是否绘制细胞间的相关性热图
group.case <- "high"
group.ctrl <- "low"
dat.type = 0

data<- read.csv("../02_DEG/normalizeExp.txt", 
                sep = "\t", row.names=1)
group <- read.csv(file = "../00_rawdata/risk.txt",
                  sep="\t")
colnames(data)<-gsub('.','-',colnames(data),fixed = T)
colnames(data)<- substring(colnames(data),1,16)
data <- data[,colnames(data) %in% group$id]
group$group <- factor(group$risk,levels = c("high", "low")) 
group <- group[,c(1,ncol(group))]
colnames(group) <- c("sample","group")
hubgene <- read.csv('../06_multicox/multiCox.csv')$id



# mark2 datprpare ---------------------------------------------------------------
if(T){
 
  expr <-  data
  group <- group
  rownames(group) <- group$sample
  Treat_sample<- group$sample[which(group$group == paste0(group.case))]
  Control_sample <- group$sample[which(group$group == paste0(group.ctrl))]
  if(ciber == 1){
    tryCatch({
      tiics_result <- readRDS('./ciber_score.rds')
    },error = function(e){
      cibersort <- deconvo_tme(eset = expr, method = "cibersort", arrays = datType, perm = 200 )#perm = 1000
      cibersort <- column_to_rownames(cibersort,var = "ID")
      colnames(cibersort) %<>% gsub("_CIBERSORT","",.)
      colnames(cibersort) %<>% gsub("_"," ",.)
      cibersort_result <- cibersort[1:22] %>% t
      saveRDS(cibersort_result,"./ciber_score.rds")
      tiics_result <- readRDS('./ciber_score.rds')
      # keep<-rowSums(tiics_result>0)>=floor(0.75*ncol(tiics_result))
      # tiics_result<-tiics_result[keep,]
    })
  }
  
  if(ssgsea == 1){
    gene_set <- read.table("E:/04/00_rawdata/00_temp/mmc3.txt", header = T, sep ="\t")
    gene_list <- split(as.matrix(gene_set)[,1], gene_set[,2])
    ssgsea_score <-  gsva(as.matrix(expr), gene_list, method = "ssgsea", kcdf='Gaussian',abs.ranking=TRUE)
    saveRDS(ssgsea_score,'./ssgsea_score.rds')
    tiics_result <- ssgsea_score  %>% as.matrix()
  }
  
  if(xcell == 1){
    xcell.dat<-deconvo_tme(eset = expr, method = "xcell", arrays = datType, perm = 200 )
    xcell.dat <- column_to_rownames(xcell.dat,var = "ID")
    colnames(xcell.dat) %<>% gsub("_xCell","",.)
    colnames(xcell.dat) %<>% gsub("_"," ",.)
    xcell.dat <- xcell.dat[1:22] %>% t
    saveRDS(xcell.dat,"./xcell_score.rds")
    tiics_result <- readRDS('./xcell_score.rds')}
  
  if(epic == 1){
	epic.dat <- deconvo_tme(eset = expr, method = "epic", arrays = datType, perm = 200)
	epic.dat <- column_to_rownames(epic.dat,var = "ID")
	colnames(epic.dat) %<>% gsub("_EPIC","",.)
	colnames(epic.dat) %<>% gsub("_"," ",.)
	epic.dat <- epic.dat[1:8] %>% t
	saveRDS(epic.dat,"./epic_score.rds")
	tiics_result <- readRDS('./epic_score.rds')
  }
  
  
  
  pvalue = padj = log2FoldChange <- matrix(0, nrow(tiics_result), 1)
  for (i in 1:nrow(tiics_result)){
    pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, Treat_sample],tiics_result[i, Control_sample])$p.value
    log2FoldChange[i, 1] = mean(tiics_result[i, Treat_sample]) - mean(tiics_result[i, Control_sample])}
  padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
  rTable <- data.frame(log2FoldChange, 
                       pvalue, 
                       padj,
                       row.names = rownames(tiics_result))
  Treat <- signif(apply(tiics_result[rownames(rTable), Treat_sample], 
                        1,
                        mean), 4)
  Control <- signif(apply(tiics_result[rownames(rTable), Control_sample], 
                          1, 
                          mean), 4)
  rTable <- data.frame(Treat, 
                       Control,
                       rTable[, c("pvalue", "pvalue", "log2FoldChange")])
  rTable$immune_cell <- rownames(rTable)
  rTable$sig <- ifelse(rTable$pvalue < 0.05,
                       ifelse(rTable$pvalue < 0.01, 
                              ifelse(rTable$pvalue < 0.001,
                                     ifelse(rTable$pvalue < 0.0001,
                                            paste(rTable$immune_cell, "****",  sep = ""),
                                            paste(rTable$immune_cell, "***", sep = "")),
                                     paste(rTable$immune_cell, "**", sep = "")),
                              paste(rTable$immune_cell, "*",  sep = "")), 
                       rTable$immune_cell)
  write.table(rTable,
              file = "immucell_wilcox_test.csv",
              quote = F,
              row.names = F,
              sep = ',')
  write.table(tiics_result,file = 'Score_result.csv',sep = ',',row.names = T,quote = F)
  diff_cibersort_Table<-rTable[which(rTable$pvalue<0.05),]
  sigcell <- rownames(diff_cibersort_Table)
  saveRDS(sigcell,'./sigCell.rds')
  write.table(diff_cibersort_Table,
              file = 'DiffTable.csv',
              quote = F,
              sep = ',')
}

#---------------------------------heatmap---------------------------------------
if(heat.map == 1){
  # tiics_result是得分结果
  rt <- as.data.frame(t(tiics_result))
  cluster <- group
  rt = rt[ , cluster[ , 1]]
  unique(cluster$sample == colnames(rt))  #check
  # rownames(cluster) = cluster$sample
  # cluster = cluster[ , -c(1, 2)]
  
  # cluster <- cluster %>% column_to_rownames("sample")
  cluster <- as.data.frame(cluster[, 2])
  rownames(cluster) = group$sample
  colnames(cluster) = "Subtype"
  #??????ͼ
  rt <- as.matrix(rt)
  library(ggsci)
  ann_color <- list(Subtype = c(normal = "#4DBBD5FF", tumor = "#DC0000FF"))
  tiff(file="HM.tiff",width = 8, height = 5, units = "in", bg = "white",
       res = 400)
  pheatmap(rt, annotation_col=cluster, 
           color = colorRampPalette(c("blue", "white", "red"))(50),
           cluster_cols =F,
           fontsize=8,
           fontsize_row=8,
           scale="row",
           show_colnames=F,
           fontsize_col=3,
           annotation_colors = ann_color)
  dev.off()
  
  pdf("HM.pdf", height = 5, width = 8)
  pheatmap(rt, annotation_col=cluster, 
           color = colorRampPalette(c("blue", "white", "red"))(50),
           cluster_cols =F,
           fontsize=8,
           fontsize_row=8,
           scale="row",
           show_colnames=F,
           fontsize_col=3,
           annotation_colors = ann_color)
  dev.off()
  
}  






# 柱状堆叠图 -------------------------------------------------------------------
# 要进入代码块一步一步运行
if(bar.cell == 1){
  bar.dat<-as.data.frame(tiics_result)
  bar.dat$cell<-rownames(bar.dat)
  box.dat <- gather(bar.dat, key=sample, value='score', -c("cell"))
  box.dat$group<-ifelse(box.dat$sample %in% Control_sample,paste0(group.ctrl),paste0(group.case))
  library(RColorBrewer)
  mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
  box_plot <- ggplot(box.dat, aes(x=sample,y=100*score,fill=cell)) +
    geom_bar(position = 'stack',stat = 'identity')+
    scale_y_continuous(expand = c(0,0))+
    theme_bw()+
    labs(x='',
         y='Relative Percent',
         fill='')+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'top') +
    # scale_color_fish(option = "Hypsypops_rubicundus", direction = -1) +
    scale_fill_manual(values = mypalette(38))+
    facet_grid(~box.dat$group,scales= "free",space= "free")
  # scale_fill_npg(22)
  box_plot
  pdf(file = paste0("cell.proportion.pdf"),width = 9,height = 6)
  a <- dev.cur()   
  png(file = paste0("cell.proportion.png"),width= 9, height= 6, units="in", res=300)
  dev.control("enable")
  par(mar = c(2,2,2,2),cex=1,family="Times")
  box_plot
  dev.copy(which = a) 
  dev.off()
  dev.off()
}

# 细胞间相关性图 -----------------------------------------------------------------
if(cell.cor == 1){
  
  if(ciber == 1){
    # 有显著性的细胞
    sigcell <- readRDS('./sigCell.rds')
    #keep <- rowSums(tiics_result > 0) >= floor(0.30 * ncol(tiics_result))
    #cellcor.dat <- tiics_result[keep,]
    #keep <- rowSums(cellcor.dat > 0) >= floor(0.30 * nrow(cellcor.dat))
    #cellcor.dat <- cellcor.dat[keep,]
    # 筛选具有显著性的细胞进行绘图
    cellcor.dat <- tiics_result
    cellcor.dat <- cellcor.dat[rownames(cellcor.dat) %in% sigcell,]
    cellcor.dat <- cellcor.dat %>% t
  } 
  if(ssgsea == 1){
    cellcor.dat <- tiics_result
    cellcor.dat <- cellcor.dat %>% t
  }
  
  
  corr<-round(cor(cellcor.dat),3)
  # p.mat<-as.data.frame(cor_pmat(cellcor.dat))
  p.mat<-cor_pmat(cellcor.dat)
  
  
  library(corrplot)
  tiff(file="cor.tiff",width = 12, height = 12, units = "in", bg = "white",
       res = 500)
  corrplot(corr, type = "upper", method = "ellipse", order = "AOE",
           tl.cex = 1.0, tl.col = "black",p.mat = p.mat,
           insig = 'blank', # 不显著留空
           col=colorRampPalette(c("blue", "white", "red"))(300))
  dev.off()
  
  
  pdf(file = "cor.pdf", height = 12, width = 12)
  corrplot(corr, type = "upper", method = "ellipse", order = "AOE",
           tl.cex = 1.0, tl.col = "black",p.mat = p.mat,
           insig = 'blank', # 不显著留空
           col=colorRampPalette(c("blue", "white", "red"))(300))
  dev.off()
  
  corr1 <- as.data.frame(corr)
  cor_p1 <- as.data.frame(p.mat)
  
  openxlsx::write.xlsx(corr1, file = "cellcor.index.xlsx",rowNames = T)
  openxlsx::write.xlsx(cor_p1, file = "cellcor.pvalue.xlsx",rowNames = T)
}

# 箱线图 ---------------------------------------------------------------------
if(box.pic == 1){
  library(ggpubr)
  dat <- tiics_result %>% t %>% as.data.frame() %>% tibble::rownames_to_column(var = "sample")
  group <- group[,"group",drop = F] %>% rownames_to_column(.,var = "sample")
  colnames(group)[2] <- "group"
  dat <- merge(group, dat, by = "sample")
  dat2 <- tidyr::gather(dat, ImmuneCell, Score, -c("sample", "group"))
  dat2$Score <- as.numeric(dat2$Score)
  dat2$group <- factor(dat2$group,levels = c(paste0(group.ctrl),paste0(group.case)))
  x1 = 1:2
  if(sample(x1,1) == 1){
    p <- ggboxplot(dat2, x = "ImmuneCell", y = "Score", 
                   color = "group",
                   palette = c("seagreen4","darkorange3"),
                   outlier.shape = NA,
                   bxp.errorbar = T
    ) +
      stat_compare_means(aes(group = group),label = "p.signif") +
      labs(x = "", y = "Score", color = "") +
      theme(axis.title = element_text(size = 20, face = "bold", family = "Times"),
            axis.text.x = element_text(size = 15, face = "bold", family = "Times", angle = 45, vjust = 1, hjust = 1),
            axis.text.y = element_text(size = 15, face = "bold", family = "Times"),
            legend.text = element_text(size = 15, family = "Times"),
            text = element_text(family = "Times"))
  }else{
    p <-  ggplot(dat2, 
                 aes(x=ImmuneCell,y=Score,
                     fill=group,
                     outlier.shape = NA,
                     bxp.errorbar = T)) + 
      geom_boxplot(width=0.8,
                   alpha=1,
                   position = position_dodge(0.9))+ 
      # annotate(geom = "text", x = rownames(rTable), y = 1, size = 5, family = "Times",
      #          label =as.character(rTable$pvalue.signif)) +
      stat_compare_means(aes(group = group),label = "p.signif") +
      theme_bw()+
      theme(legend.position = "top")+
      theme(axis.title.x =element_text(size=16,family = "Times", face = "bold"),
            axis.text.x =element_text(angle=45,size=16,hjust = 1,family = "Times", face = "bold"),
            axis.title.y =element_text(size=20,family = "Times", face = "bold"),
            axis.text.y=element_text(size=16,family = "Times", face = "bold"))+
      theme(
        legend.title=element_text(size=15) , legend.text=element_text(size=14))+
      labs(x="",y="Score")+
      ggsci::scale_fill_d3()
  }
  ggsave(filename = "ImmuneCell_boxplot.pdf", width = 14, height = 8, p)
  ggsave(filename = "ImmuneCell_boxplot.png", width = 14, height = 8, p)
}

# Correlatin with SigCell ------------------------------------------------


if(cor.map == 1){
  ssgsea_score <- readRDS("./ciber_score.rds") %>% as.data.frame()
  library(tibble)
  library(magrittr)
  library(corrplot)
  library(reshape2)
  filte_dat<- expr
  filte_dat <- filte_dat[,match(group$sample,colnames(filte_dat))]
  filte_dat <- filte_dat %>% as.data.frame() %>% rownames_to_column(.,var = "Gene.name")
  tem_12 <- filte_dat[match(hubgene,filte_dat$Gene.name),] 
  keyRNA_exp <- tem_12[,-1]
  rownames(keyRNA_exp) <- tem_12$Gene.name
  keyRNA_exp <- keyRNA_exp %>% t %>% as.data.frame()
  # sigcell --------------------------------------------------------------------=
  tiics_result <- readRDS("./ciber_score.rds") %>% t %>% as.data.frame()
  sigcell_exp <- tiics_result  %>% as.data.frame()
  sigcell_exp <- sigcell_exp[group$sample,]
  cor<-stats::cor
  #基因在列，样品在行
  corr <- round(cor(sigcell_exp,keyRNA_exp, method = "spearman"), 3) %>% na.omit() 
  cor_p <- WGCNA::corPvalueStudent(corr, nrow(sigcell_exp)) %>% round(.,4) 
  corr2 <- as.data.frame(corr)
  cor_p2 <- as.data.frame(cor_p)
  openxlsx::write.xlsx(corr2, file = "gene_with_cellcor_index.xlsx",rowNames = T)
  openxlsx::write.xlsx(cor_p2, file = "gene_with_cellcor_pvalue.xlsx",rowNames = T)
  
  library(reshape2)
  temcor <- melt(corr,value.name="cor")
  temcor <-  temcor[abs(temcor$cor) > 0.3,]
  temcorp <- melt(cor_p,value.name="pvalue")
  temcorp <- temcorp[temcorp$pvalue < 0.05,]
  colnames(temcorp)[1:2] <- c("Var1","Var2")
  tem <- merge(temcor,temcorp,all = T) %>% na.omit()
  tem <- tem[order(tem$Var2,decreasing = F),]
  tem$pvalue[tem$pvalue == 0] <- "<0.001"
  
  library(RColorBrewer)
  library(pheatmap)
  col1 <- colorRampPalette(c("#53a7db" ,"#ddeff6","#EE0000B2"),alpha = TRUE)(100)
  cor <- corr
  pvalue <- cor_p
  myfun <- function(pval) {
    stars = ""
    if(pval <= 0.001)
      stars = "***"
    if(pval > 0.001 & pval <= 0.01)
      stars = "**"
    if(pval > 0.01 & pval <= 0.05)
      stars = "*"
    if(pval > 0.05 & pval <= 0.1)
      stars = ""
    stars
  }
  heatmap <- melt(cor) %>%
    mutate(pvalue=melt(pvalue)[,3]) %>%
    mutate(signif = sapply(pvalue, function(x) myfun(x)))
  cor_char <- as.character(round(heatmap$value,2))
  heatmap$signif <- paste0(cor_char,"\n",heatmap$signif)
  colnames(heatmap) <- c("Cell","Gene","cor","pvalue","signif")  
  
  p1 <- ggplot(heatmap,aes(Gene,Cell,fill = cor))+
    geom_tile()+theme_minimal()+
    geom_text(aes(label=signif),size=5,color="black",
              hjust=0.5,vjust=0.5)+
    labs(x = NULL,y = NULL,color=NULL) + 
    scale_fill_viridis_c()+
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme(text=element_text(family="Times"),
          axis.ticks.x = element_blank(),
          axis.ticks.y=element_blank(),
          panel.border = element_rect(fill=NA,color="grey70",
                                      size=1, linetype="solid")) +
    # scale_size(range=c(1,15),guide=NULL)+
    guides(color = guide_colorbar(direction = "vertical",
                                  reverse = F,barwidth = unit(.5, "cm"),
                                  barheight = unit(9, "cm")))+
    scale_fill_gradientn(colours=col1, limits=c(-1, 1),
                         name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"),
                         guide=guide_colourbar(ticks=T, 
                                               nbin=50,
                                               barheight=10, 
                                               label=T,
                                               barwidth=1))+
    theme(legend.position = "right")+
    theme(axis.title.x =element_text(size=16,family = "Times", face = "bold"),
          axis.text.x =element_text(angle=45,size=16,hjust = 1,family = "Times", face = "bold"),
          axis.title.y =element_text(size=16,family = "Times", face = "bold"),
          axis.text.y=element_text(size=16,family = "Times", face = "bold"))+
    theme(
      legend.title=element_text(size=15) , legend.text=element_text(size=14))+
    labs(x="",y="")+coord_flip()
  pdf(file = paste0("gene_correlation.pdf"),width=18, height=8)
  a <- dev.cur()   
  png(file = paste0("gene_correlation.png"),width=18, height=8, units="in", res=300)
  dev.control("enable")
  print(p1)
  dev.copy(which = a) 
  dev.off()
  dev.off()
}

sigcell <- readRDS('./sigCell.rds')
length(sigcell)
# ciber 8

paste0(sigcell,collapse = '，')
# [1] "T cells CD4 naive，T cells CD4 memory resting，T cells CD4 memory activated，T cells follicular helper，T cells regulatory (Tregs)，Monocytes，Macrophages M0，Mast cells resting"
save.image('./immu.Rdata')


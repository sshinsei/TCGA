rm(list=ls())
foldChange = 1
padj = 0.05
setwd('E:/LZ/24080')
if (! dir.exists('./17_high_low_enrich')){
  dir.create('./17_high_low_enrich')
}
setwd('./17_high_low_enrich')
library(edgeR)
library(tidyverse)
rt = read.csv("../01_TCGA/files/mRNA_count.txt", sep="\t", header=T, check.names=F) 
cluster = read.table("../00_rawdata/risk.txt", 
                      sep="\t", header=T, check.names=F)
cluster = cluster[order(cluster$risk, decreasing = T), ]
# 修改样本名 ：与risk中保持一致
colnames(rt)[2:ncol(rt)] = substring(colnames(rt)[2:ncol(rt)], 1, 16)
# 提取rt的cluster
inte = intersect(colnames(rt)[2:ncol(rt)], cluster$id)
cluster1 = cluster[cluster$id %in% inte,]


#cluster = cluster[cluster$cluster == "1" | cluster$cluster == "2", ]
# rt = rt[ ,!grepl("11A",colnames(rt)[2:ncol(rt)])]  

#delete normal sample
# rt = rt[ ,!grepl("11A",colnames(rt)[2:ncol(rt)])] 
# rt <- rt[,colnames(rt)[2:ncol(rt)] %in% cluster1$id]
# 转置
rt <- as.data.frame(t(rt))
colnames(rt) <- rt[1,]
rt <- rt[-1,]
# 提取rt的cluster
inte = intersect(rownames(rt), cluster$id)
cluster1 = cluster[cluster$id %in% inte,]


cols = na.omit(match(cluster1$id, rownames(rt)))
rt1 = rt[cols, ]
# rm(list=c("rt"))


# cols = na.omit(match(cluster1$id, colnames(rt)[2:ncol(rt)]))
# rt1 = rt[ , c(1, cols)]
unique(rownames(rt1) == cluster1$id)  #check
table(cluster1$risk)
# high  low 
# 172   172

# 转置
rt1 <- as.data.frame(t(rt1))

  
rt1 = as.matrix(rt1)
# rownames(rt1) = rt1[ , 1]
exp = rt1
rm(list=c("rt", "rt1"))
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rm(list=c("exp"))
data=avereps(data)
data=data[rowMeans(data)>1, ]

# 标准化
data <- log2(data + 1)

group=c(rep("low", 172), rep("high", 172))  
design <- model.matrix(~group)
y <- DGEList(counts = data, group = group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y, pair = c("low", "high")) # 高风险组相对于低风险组上下调
topTags(et)
ordered_tags <- topTags(et, n = 100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR) == FALSE, ]
diff=allDiff
newData=y$pseudo.counts

write.table(diff,file="edgerOut.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp_highlow.txt",sep="\t",quote=F,col.names=F)   #???????л???У?????ı???ֵ??normalizeExp.txt??
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)         #????????????У?????ı???ֵ??diffmRNAExp.txt??

heatmapData <- newData[rownames(diffSig),]

#volcano

tiff(file="vol.tiff",width = 10,height =10,units ="cm",compression="lzw",
     bg="white",res=400)
xMax=max(-log10(allDiff$FDR))+1
yMax=12
plot(-log10(allDiff$FDR), allDiff$logFC, xlab="-log10(FDR)",ylab="logFC",
     xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC>foldChange,]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC<(-foldChange),]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()

pdf(file="vol.pdf")
xMax=max(-log10(allDiff$FDR))+1
yMax=12
plot(-log10(allDiff$FDR), allDiff$logFC, xlab="-log10(FDR)",ylab="logFC",
     xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC>foldChange,]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC<(-foldChange),]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()

#heatmap
hmExp=log10(heatmapData+0.001)
library('gplots')
hmMat=as.matrix(hmExp)
#pdf(file="heatmap.pdf",width=60,height=90)
tiff(file="heatmap.tiff",width =60,height =90,units ="cm",compression="lzw",bg="white",res=400)
par(oma=c(10,3,3,7))
heatmap.2(hmMat,col='greenred',trace="none")
dev.off()


pdf(file="heatmap.pdf",width=60,height=90)
par(oma=c(10,3,3,7))
heatmap.2(hmMat,col='greenred',trace="none")
dev.off()






##################################################################################
#############################  enrich  ###########################################
##################################################################################
rm(list = ls())
####富集数据文件准备####
#载入DESeq2差异分析结果表：
setwd('E:/LZ/24080/17_high_low_enrich')
if (! dir.exists('./all')){
  dir.create('./all')
}
setwd('./all')
library(clusterProfiler)
library(stringr)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(ggplot2)
options(stringsAsFactors = F)
load("go.Rdata")
rt1 <- read.table('../down.xls', 
                  sep="\t", header=T, check.names=F)
rt2 <- read.table('../up.xls', 
                  sep="\t", header=T, check.names=F)
rt1 <- rt1 %>% 
  rownames_to_column("id")
rt2 <- rt2 %>% 
  rownames_to_column("id")
DESeq2 <- rbind(rt1,rt2)

head(DESeq2)
length(rownames(DESeq2)) # 299
#添加entrez ID列：
##symbol转entrez ID：
enrich_ID <- DESeq2$id

gene_symbol <- clusterProfiler::bitr(
  geneID = enrich_ID,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db")

######  03.1_GO######
enrich_go <- enrichGO(gene = gene_symbol[,2],
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "ALL",
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)
enrich_go <- mutate(enrich_go, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

#保存GO注释文件
GO_ALL <- enrich_go@result
GO_ALL <- GO_ALL[GO_ALL$pvalue <0.05,]
BP <- GO_ALL[GO_ALL$ONTOLOGY=='BP', ]
BP <- BP[order(BP$Count,decreasing = T),]
CC <- GO_ALL[GO_ALL$ONTOLOGY=='CC', ]
CC <- CC[order(CC$Count,decreasing = T),]
MF <- GO_ALL[GO_ALL$ONTOLOGY=='MF', ]
MF <- MF[order(MF$Count,decreasing = T),]
paste0("得到",dim(GO_ALL)[[1]],"个结果，其中",dim(BP)[[1]],"个生物学过程，",dim(CC)[[1]],"个细胞组分，",dim(MF)[[1]],"个分子功能")
# [1] "得到72个结果，其中64个生物学过程，4个细胞组分，4个分子功能"

write.table(GO_ALL,file = "00.go_ALL.txt",sep = "\t",quote = F,row.names = F)
write.table(as.data.frame(BP), '00.GO_BP.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(as.data.frame(CC), '00.GO_CC.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(as.data.frame(MF), '00.GO_MF.txt', sep = '\t', row.names = FALSE, quote = FALSE)


display_number = c(64, 4, 4)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
ego_result_BP <- as.data.frame(BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(MF)[1:display_number[3], ]
paste0(ego_result_BP$Description,collapse = "（）；")
# "anterior/posterior pattern specification（）；regionalization（）；pattern specification process（）；blastocyst development"

paste0(ego_result_CC$Description,collapse = "（）；")
# [1] "synaptic membrane（）；postsynaptic membrane（）；postsynaptic specialization membrane（）；postsynaptic density membrane"

paste0(ego_result_MF$Description,collapse = "（）；")
# [1] "DNA-binding transcription activator activity, RNA polymerase II-specific（）；DNA-binding transcription activator activity（）；neurotransmitter receptor activity（）；glutamate receptor activity"

#GO柱状图-合并
# BP&MF:top25
if(T){
  go_enrich_df <- data.frame(
    ID=c(ego_result_BP[1:25,]$ID, ego_result_CC$ID, ego_result_MF$ID),
    Description=c(ego_result_BP[1:25,]$Description,ego_result_CC$Description,ego_result_MF$Description),
    GeneNumber=c(ego_result_BP[1:25,]$Count, ego_result_CC$Count,ego_result_MF$Count),
    type=factor(c(rep("BP:biological process", 25), 
                  rep("CC:cellular component", display_number[2]),
                  rep("MF:molecular function", display_number[3])), 
                levels=c("BP:biological process", "CC:cellular component","MF:molecular function" )))
  
  
  ##开始绘制GO柱状图
  ###横着的柱状图
  go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
  COLS <- c("#EFC000B2", "#4DBBD5B2", "#00A087B2")#设定颜色
  
  p <- ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
    geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
    scale_fill_manual(values = COLS) + ###颜色
    coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
    xlab("GO term") + 
    ylab("Gene Number") + 
    labs(title = "GO Terms")+
    theme_bw()+
    theme(
      legend.background = element_rect(fill = "white", color = "black", size = 0.2),
      legend.text = element_text(face="bold",color="black",family = "Times",size=14),
      plot.title = element_text(hjust = 0.5, face = "bold",color = "black",family = "Times",size = 16),
      axis.text.x = element_text(face = "bold",color = "black",size = 12),
      axis.text.y = element_text(face = "bold",color = "black",size = 10),
      axis.title.x = element_text(face = "bold",color = "black",family = "Times",size = 16),
      axis.title.y = element_text(face = "bold",color = "black",family = "Times",size = 16),
      plot.subtitle = element_text(hjust = 0.5,family = "Times", size = 14, face = "italic", colour = "black")
    )
  ggsave(filename = "01.go_barplot.pdf", height = 8, width = 14, p)
  ggsave(filename = "01.go_barplot.png", height = 8, width = 14, p)
} 
library(openxlsx)
# circ
if(T){
  select_DERNA <- read.table("../diffSig.xls",row.names = 1)
  select_DERNA <- select_DERNA[rownames(select_DERNA) %in% enrich_ID,]
  
  # 根据p.value排序
  #ego_result_BP <- ego_result_BP[order(ego_result_BP$p.adjust,decreasing = F),]
  #ego_result_CC <- ego_result_CC[order(ego_result_CC$p.adjust,decreasing = F),]
  #ego_result_MF <- ego_result_MF[order(ego_result_MF$p.adjust,decreasing = F),]
  
  
  go_tem <- rbind(ego_result_BP[1:10,],ego_result_CC,ego_result_MF[1:10,])
  go <-  data.frame(Category = go_tem[,'ONTOLOGY'],ID = go_tem[,'ID'],Term = go_tem[,'Description'], 
                    Genes = gsub("/", ", ", go_tem[,'geneID']), adj_pval = go_tem[,'p.adjust'])
  genelist <- data.frame(ID = rownames(select_DERNA), logFC = select_DERNA$logFC)
  
  row.names(genelist)=genelist[,1]
  library(GOplot)
  circ <- circle_dat(go, genelist)
  # p <- GOCircle(circ,rad1=2.5,rad2=3.5,label.size= 5,nsub=10)
  go_tem$p.adjust <- -log10(go_tem$p.adjust)
  ##将以上我们摘取的部分通路重新组合成数据框
  go_enrich_df <- data.frame(
    category=go_tem$ONTOLOGY, 
    gene_num.min = 0,
    gene_num.max = 200,
    gene_num.rich=go_tem$Count,
    "-log10.p" = go_tem$p.adjust,
    up.regulated = 0,
    down.regulated = 0,
    rich.factor = go_tem$richFactor,
    Description=go_tem$ID)
  circ$sig <- ifelse(circ$logFC > 0,"up","down")
  tem1 <- go_enrich_df$Description
  
  for(i in 1:28){
    x = tem1[i]
    circ1 <- circ[grep(x,circ$ID),]
    num1 <- nrow(circ1)
    tem2 <- as.data.frame(table(circ1$sig))
    if(nrow(tem2) == 2){
      go_enrich_df[match(x,go_enrich_df$Description),7] = table(circ1$sig)[[1]]
      go_enrich_df[match(x,go_enrich_df$Description),6] = table(circ1$sig)[[2]]
    }else if(tem2[1,1] == "up"){
      go_enrich_df[match(x,go_enrich_df$Description),7] = table(circ1$sig)[[1]]
    }else{
      go_enrich_df[match(x,go_enrich_df$Description),6] = table(circ1$sig)[[1]]}
  }
  
  #读取示例数据，富集结果
  #首先给 ko id 排个序，默认按照原表格中的排列顺序
  go_enrich_df$Description <- factor(go_enrich_df$Description, levels = go_enrich_df$Description)
  rownames(go_enrich_df) <- go_enrich_df$Description
  pdf(file = paste0("01_GOcirclize_plot.pdf"),width = 12,height = 6)
  a <- dev.cur()   
  png(file = paste0("01_GOcirclize_plot.png"),width= 12, height= 6, units="in", res=300)
  dev.control("enable")
  par(mar = c(2,2,2,2),cex=1.5)
  circle_size = unit(1, 'snpc')
  #加载 circlize 包
  library(circlize)
  ##整体布局
  circos.par(gap.degree = 2, start.degree = 90)
  ##第一圈，绘制 ko
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max')]  #选择作图数据集，定义了 ko 区块的基因总数量范围
  ko_color <- c(rep('#EFC000B2', 7), rep('#4DBBD5B2', 7), rep('#00A087B2', 7))  #定义分组颜色
  circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)  #一个总布局
  circos.track(
    ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = ko_color,  #圈图的高度、颜色等设置
    panel.fun = function(x, y) {
      ylim = get.cell.meta.data('ycenter')  #ylim、xlim 用于指定 ko id 文字标签添加的合适坐标
      xlim = get.cell.meta.data('xcenter')
      sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称
      circos.axis(h = 'top', labels.cex = 0.4, major.tick.length = 0.4, labels.niceFacing = FALSE)  #绘制外周的刻度线
      circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  #将 ko id 文字标签添加在图中指定位置处
    } )
  colnames(go_enrich_df)[5] <- c("-log10(p.adj)")
  ##第二圈，绘制富集的基因和富集 p 值
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.rich', '-log10(p.adj)')]  #选择作图数据集，包括富集基因数量以及 p 值等信息
  label_data <- go_enrich_df['gene_num.rich']  #标签数据集，仅便于作图时添加相应的文字标识用
  p_max <- round(max(go_enrich_df$`-log10(p.adj)`)) + 1  #定义一个 p 值的极值，以方便后续作图
  colorsChoice <- colorRampPalette(c('white','#4393C3'))  #这两句用于定义 p 值的渐变颜色
  color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #区块的长度反映了富集基因的数量，颜色与 p 值有关
      ylim = get.cell.meta.data('ycenter')  #同上文，ylim、xlim、sector.name 等用于指定文字标签（富集基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),1]
      circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  #将文字标签添（富集基因数量）加在图中指定位置处
    } )
  
  ##第三圈，绘制上下调基因
  #首先基于表格中上下调基因的数量，计算它们的占比
  go_enrich_df$all.regulated <- go_enrich_df$up.regulated + go_enrich_df$down.regulated
  go_enrich_df$up.proportion <- go_enrich_df$up.regulated / go_enrich_df$all.regulated
  go_enrich_df$down.proportion <- go_enrich_df$down.regulated / go_enrich_df$all.regulated
  
  #随后，根据上下调基因的相对比例，分别计算它们在作图时的“区块坐标”和“长度”
  go_enrich_df$up <- go_enrich_df$up.proportion * go_enrich_df$gene_num.max
  plot_data_up <- go_enrich_df[c('Description', 'gene_num.min', 'up')]
  names(plot_data_up) <- c('Description', 'start', 'end')
  plot_data_up$type <- 1  #分配 1 指代上调基因
  
  go_enrich_df$down <- go_enrich_df$down.proportion * go_enrich_df$gene_num.max + go_enrich_df$up
  plot_data_down <- go_enrich_df[c('Description', 'up', 'down')]
  names(plot_data_down) <- c('Description', 'start', 'end')
  plot_data_down$type <- 2  #分配 2 指代下调基因
  
  #选择作图数据集（作图用）、标签数据集（添加相应的文字标识用），并分别为上下调基因赋值不同颜色
  plot_data <- rbind(plot_data_up, plot_data_down)
  label_data <- go_enrich_df[c('up', 'down', 'up.regulated', 'down.regulated')]
  color_assign <- colorRamp2(breaks = c(1, 2), col = c('#FDDBC7', '#9ACD32'))
  #继续绘制圈图
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #这里紫色代表上调基因，蓝色代表下调基因，区块的长度反映了上下调基因的相对占比
      ylim = get.cell.meta.data('cell.bottom.radius') - 0.5  #同上文，ylim、xlim、sector.name 等用于指定文字标签（上调基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),3]
      circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  #将文字标签（上调基因数量）添加在图中指定位置处
      xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),4]
      circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  #类似的操作，将下调基因数量的标签也添加在图中
    } )
  # go_enrich_df$`-log10(pvalue)` <- 10*go_enrich_df$`-log10(pvalue)`
  ##第四圈，绘制富集得分
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max', 'rich.factor')]  #选择作图数据集，标准化后的富集得分
  label_data <- go_enrich_df['category']  #将通路的分类信息提取出，和下一句一起，便于作图时按分组分配颜色
  color_assign <- c('BP' = '#EFC000B2', 'CC' = '#4DBBD5B2', 'MF' = '#00A087B2')
  #E64B35B2" "#4DBBD5B2" "#00A087B2
  circos.genomicTrack(
    plot_data, ylim = c(0, 0.5), track.height = 0.3, bg.col = 'gray95', bg.border = NA,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称，并添加在下一句中匹配 ko 对应的高级分类，以分配颜色
      circos.genomicRect(region, value, col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...)  #绘制矩形区块，高度代表富集得分，颜色代表 ko 的分类
      circos.lines(c(0, max(region)), c(0.05, 0.05), col = 'gray', lwd = 0.3)  #可选在富集得分等于 0.5 的位置处添加一个灰线
    } )
  ##绘图完毕后，不要忘了清除痕迹，以免影响下一次作图
  circos.clear()
  
  category_legend <- Legend(#E64B35B2" "#4DBBD5B2" "#00A087B2
    labels = c("BP:biological process", "CC:cellular component","MF:molecular function" ),
    type = 'points', pch = NA, background = c('#EFC000B2', '#4DBBD5B2', '#00A087B2'), 
    labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  
  updown_legend <- Legend(
    labels = c('Up-regulated', 'Down-regulated'), 
    type = 'points', pch = NA, background = c('#FDDBC7', '#9ACD32'), 
    
    labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  
  pvalue_legend <- Legend(
    col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                         colorRampPalette(c('white','#4393C3'))(6)),
    legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
    title_gp = gpar(fontsize = 9), title_position = 'topleft', title = '-Log10(p.adj)')
  
  lgd_list_vertical <- packLegend(updown_legend, pvalue_legend)
  pushViewport(viewport(x = 0.85, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  
  lgd_list_vertical <- packLegend(category_legend)
  pushViewport(viewport(x = 0.5, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  dev.copy(which = a) 
  dev.off()
  dev.off()
}

#--------------------------chord plot:25097--------------------------------------
library(GOplot)
library(RColorBrewer)
library(openxlsx)
library(tidyverse)
df <- GO_ALL
df$geneID <- str_replace_all(df$geneID,'/',',')
df <- df[,c(2,3,10,12,1)] # 
names(df) <- c('ID','Term','adj_pval',"Genes",'Category')
# 选取富集结果p.adj最显著的前7可视化
df <- df[order(df$adj_pval,decreasing = F),]
df <- df[1:10,]
select_DERNA <- read.table("../diffSig.xls",row.names = 1)


genelist <- data.frame(ID = rownames(select_DERNA), logFC = select_DERNA$logFC)
use_gene <- genelist %>% filter(ID %in% enrich_ID)

# 格式化基因信息
df$Genes <- str_replace_all(df$Genes, "\\s+", "")  # 去除空格
use_gene$ID <- toupper(trimws(use_gene$ID))       # 转为大写，去除空格

circ <- circle_dat(df,use_gene)
chord <- chord_dat(circ,use_gene) #logFC NA
p <- GOChord(data = chord,title="GO enrichment",space=0.01,limit=c(1,1),
             gene.order="logFC",gene.space = 0.25,gene.size = 5,
             # lfc.col = c('red',"white",'blue'),
             # ribbon.col = c('#f5e56b','#c0da6e','#6fc198','#727aa3','#765280'),
             ribbon.col = brewer.pal(10, "Set3"),
             process.label = 10) +
  theme(plot.title = element_text(hjust = 0.5,size=10))

print(p)

pdf(file = paste0("01.GO_chord.pdf"),width = 15,height = 15)
p
dev.off()

png(file = paste0("01.GO_chord.png"),width = 15,height = 15, units="in", res=300)
p
dev.off()
save.image('go.Rdata')


###### 03.2_KEGG ---------------------------------------------------------------
rm(list = ls())
library(clusterProfiler)
library(stringr)
library(ComplexHeatmap)
library(org.Hs.eg.db)
options(stringsAsFactors = F)

rt1 <- read.table('../down.xls', 
                  sep="\t", header=T, check.names=F)
rt2 <- read.table('../up.xls', 
                  sep="\t", header=T, check.names=F)
rt1 <- rt1 %>% 
  rownames_to_column("id")
rt2 <- rt2 %>% 
  rownames_to_column("id")
DESeq2 <- rbind(rt1,rt2)

head(DESeq2)
length(rownames(DESeq2)) # 299
#添加entrez ID列：
##symbol转entrez ID：
enrich_ID <- DESeq2$id

gene_symbol <- clusterProfiler::bitr(
  geneID = enrich_ID,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db")
KEGG <- enrichKEGG( gene = gene_symbol$ENTREZID,   #基因列表 
                    organism = "hsa",  #物种
                    keyType = "kegg",  #指定的基因ID类型，
                    minGSSize = 1, 
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,  
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05
)
KEGG <- mutate(KEGG, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

kk <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg_result <- kk@result
hh <- as.data.frame(kegg_result)
rownames(hh) <- 1:nrow(hh)
hh <- hh[hh$pvalue <= 0.05,]
dim(hh)# 23 15
hh <- hh[order(hh$Count,decreasing = T),]
write.table(hh,file = "00.KEGG.txt",sep = "\t",quote = F,row.names = F)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
kk <- hh
paste0(kk$Description,collapse = "（）；")
# "Neuroactive ligand-receptor interaction（）；Calcium signaling pathway（）；Retrograde endocannabinoid signaling（）；Ras signaling pathway（）；

str(hh)
select_DERNA <- read.table("../diffSig.xls")
go <-  data.frame(Category = "KEGG",ID = kk[,'ID'],Term = kk[,'Description'], 
                  Genes = gsub("/", ", ", kk[,'geneID']), adj_pval = kk[,'p.adjust'])
genelist <- data.frame(ID = rownames(select_DERNA), logFC = select_DERNA$logFC)

row.names(genelist)=genelist[,1]
library(GOplot)
circ <- circle_dat(go, genelist)
p <- GOCircle(circ,rad1=2.5,rad2=3.5,label.size= 5,nsub=10)
kk$p.adjust <- -log10(kk$p.adjust)
##将以上我们摘取的部分通路重新组合成数据框
go_enrich_df <- data.frame(
  category="KEGG", 
  gene_num.min = 0,
  gene_num.max = 100,
  gene_num.rich=kk$Count,
  "-log10.p" = kk$p.adjust,
  up.regulated = 0,
  down.regulated = 0,
  rich.factor = kk$richFactor,
  Description=kk$ID)
circ$sig <- ifelse(circ$logFC > 0,"up","down")

tem1 <- go_enrich_df$Description

for(i in 1:10){
  x = tem1[i]
  circ1 <- circ[grep(x,circ$ID),]
  num1 <- nrow(circ1)
  tem2 <- as.data.frame(table(circ1$sig))
  if(nrow(tem2) == 2){
    go_enrich_df[match(x,go_enrich_df$Description),7] = table(circ1$sig)[[1]]
    go_enrich_df[match(x,go_enrich_df$Description),6] = table(circ1$sig)[[2]]
  }else if(tem2[1,1] == "up"){
    go_enrich_df[match(x,go_enrich_df$Description),7] = table(circ1$sig)[[1]]
  }else{
    go_enrich_df[match(x,go_enrich_df$Description),6] = table(circ1$sig)[[1]]}
}

#读取示例数据，富集结果
#首先给 ko id 排个序，默认按照原表格中的排列顺序
go_enrich_df$Description <- factor(go_enrich_df$Description, levels = go_enrich_df$Description)
rownames(go_enrich_df) <- go_enrich_df$Description
library(ggplot2)

# circ
if(T){
  ####创建一个 pdf 画板
  
  pdf(file = paste0("02.KEGG_cic.pdf"),width = 13,height = 6)
  a <- dev.cur()   
  png(file = paste0("02.KEGG_cic.png"),width= 13, height= 6, units="in", res=300)
  dev.control("enable")
  par(family="Times")
  circle_size = unit(1, 'snpc')
  
  #加载 circlize 包
  library(circlize)
  
  ##整体布局
  circos.par(gap.degree = 2, start.degree = 90)
  ##第一圈，绘制 ko
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max')]  #选择作图数据集，定义了 ko 区块的基因总数量范围
  ko_color <- c(rep('#F7CC13', 10))  #定义分组颜色
  circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)  #一个总布局
  circos.track(
    ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = ko_color,  #圈图的高度、颜色等设置
    panel.fun = function(x, y) {
      ylim = get.cell.meta.data('ycenter')  #ylim、xlim 用于指定 ko id 文字标签添加的合适坐标
      xlim = get.cell.meta.data('xcenter')
      sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称
      circos.axis(h = 'top', labels.cex = 0.4, major.tick.length = 0.4, labels.niceFacing = FALSE)  #绘制外周的刻度线
      circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  #将 ko id 文字标签添加在图中指定位置处
    } )
  colnames(go_enrich_df)[5] <- c("-log10(padj)")
  ##第二圈，绘制富集的基因和富集 p 值
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.rich', '-log10(padj)')]  #选择作图数据集，包括富集基因数量以及 p 值等信息
  label_data <- go_enrich_df['gene_num.rich']  #标签数据集，仅便于作图时添加相应的文字标识用
  p_max <- round(max(go_enrich_df$`-log10(padj)`)) + 1  #定义一个 p 值的极值，以方便后续作图
  colorsChoice <- colorRampPalette(c('white','#4393C3'))  #这两句用于定义 p 值的渐变颜色
  color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))
  
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #区块的长度反映了富集基因的数量，颜色与 p 值有关
      ylim = get.cell.meta.data('ycenter')  #同上文，ylim、xlim、sector.name 等用于指定文字标签（富集基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),1]
      circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  #将文字标签添（富集基因数量）加在图中指定位置处
    } )
  
  ##第三圈，绘制上下调基因
  #首先基于表格中上下调基因的数量，计算它们的占比
  go_enrich_df$all.regulated <- go_enrich_df$up.regulated + go_enrich_df$down.regulated
  go_enrich_df$up.proportion <- go_enrich_df$up.regulated / go_enrich_df$all.regulated
  go_enrich_df$down.proportion <- go_enrich_df$down.regulated / go_enrich_df$all.regulated
  
  #随后，根据上下调基因的相对比例，分别计算它们在作图时的“区块坐标”和“长度”
  go_enrich_df$up <- go_enrich_df$up.proportion * go_enrich_df$gene_num.max
  plot_data_up <- go_enrich_df[c('Description', 'gene_num.min', 'up')]
  names(plot_data_up) <- c('Description', 'start', 'end')
  plot_data_up$type <- 1  #分配 1 指代上调基因
  
  go_enrich_df$down <- go_enrich_df$down.proportion * go_enrich_df$gene_num.max + go_enrich_df$up
  plot_data_down <- go_enrich_df[c('Description', 'up', 'down')]
  names(plot_data_down) <- c('Description', 'start', 'end')
  plot_data_down$type <- 2  #分配 2 指代下调基因
  
  #选择作图数据集（作图用）、标签数据集（添加相应的文字标识用），并分别为上下调基因赋值不同颜色
  plot_data <- rbind(plot_data_up, plot_data_down)
  label_data <- go_enrich_df[c('up', 'down', 'up.regulated', 'down.regulated')]
  color_assign <- colorRamp2(breaks = c(1, 2), col = c('#FDDBC7', '#9ACD32'))
  
  #继续绘制圈图
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #这里紫色代表上调基因，蓝色代表下调基因，区块的长度反映了上下调基因的相对占比
      ylim = get.cell.meta.data('cell.bottom.radius') - 0.5  #同上文，ylim、xlim、sector.name 等用于指定文字标签（上调基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),3]
      circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  #将文字标签（上调基因数量）添加在图中指定位置处
      xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),4]
      circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  #类似的操作，将下调基因数量的标签也添加在图中
    } )
  # go_enrich_df$`-log10(padj)` <- 10*go_enrich_df$`-log10(pvalue)`
  ##第四圈，绘制富集得分
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max', 'rich.factor')]  #选择作图数据集，标准化后的富集得分
  label_data <- go_enrich_df['category']  #将通路的分类信息提取出，和下一句一起，便于作图时按分组分配颜色
  color_assign <- c('KEGG' = '#FFA500')
  
  circos.genomicTrack(
    plot_data, ylim = c(0, 0.5), track.height = 0.3, bg.col = 'gray95', bg.border = NA,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称，并添加在下一句中匹配 ko 对应的高级分类，以分配颜色
      circos.genomicRect(region, value, col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...)  #绘制矩形区块，高度代表富集得分，颜色代表 ko 的分类
      circos.lines(c(0, max(region)), c(0.05, 0.05), col = 'gray', lwd = 0.3)  #可选在富集得分等于 0.5 的位置处添加一个灰线
    } )
  
  ##绘图完毕后，不要忘了清除痕迹，以免影响下一次作图
  circos.clear()
  
  category_legend <- Legend(
    labels = c("KEGG Pathway"),
    type = 'points', pch = NA, background = c('#FFA500'), 
    labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  
  updown_legend <- Legend(
    labels = c('Up-regulated', 'Down-regulated'), 
    type = 'points', pch = NA, background =c('#FDDBC7', '#9ACD32'), 
    labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  
  pvalue_legend <- Legend(
    col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                         colorRampPalette(c('white','#4393C3'))(6)),
    legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
    title_gp = gpar(fontsize = 9), title_position = 'topleft', title = '-Log10(padj)')
  
  lgd_list_vertical <- packLegend(updown_legend, pvalue_legend)
  pushViewport(viewport(x = 0.85, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  
  lgd_list_vertical <- packLegend(category_legend)
  pushViewport(viewport(x = 0.5, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  dev.copy(which = a) 
  dev.off()
  dev.off()
}

if(T){
  # 气泡图
  paste0(kk$Description,collapse = "（）；")#显示出10个功能的名称
  kk <- hh
  kk <- na.omit(kk)
  # 绘图，结果如右图所示，保存结果，BP的分析就完成了，然后就是CC、MF的绘制，过程和BP的绘制过程一致
  p2 <-  ggplot(kk,aes(y=order,x=Count))+
    geom_point(aes(size=Count,color= p.adjust))+# 修改点的大小
    scale_color_gradient(high="#EE0000B2",low = "#008B45B2")+
    labs(color=expression(p.Value,size="Count"), 
         x="Gene Number",y="Pathways",title="KEGG Enrichment")+
    theme_bw()+
    theme(
      legend.background = element_rect(fill = "white", color = "black", size = 0.2),
      legend.text = element_text(face="bold",color="black",family = "Times",size=14),
      plot.title = element_text(hjust = 0.5, face = "bold",color = "black",family = "Times",size = 18),
      axis.text.x = element_text(face = "bold",color = "black",size = 14),
      axis.text.y = element_text(face = "bold",color = "black",size = 12),
      axis.title.x = element_text(face = "bold",color = "black",family = "Times",size = 18),
      axis.title.y = element_text(face = "bold",color = "black",family = "Times",size = 18),
      plot.subtitle = element_text(hjust = 0.5,family = "Times", size = 14, face = "bold", colour = "black")
    )+scale_y_discrete(labels=function(x) str_wrap(x, width=40))
  ggsave(filename = "02.KEGG_dot.pdf", height = 8, width = 9, p2)
  ggsave(filename = "02.KEGG_dot.png", height = 8, width = 9, p2)
}
save.image('kegg.Rdata')




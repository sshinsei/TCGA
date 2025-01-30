###??TCGA???ݿ????? ѡ???֣?program-TCGA????simple nucleotide variation????Masked Somatic Mutation
###????cart
rm(list=ls())
setwd('E:/LZ/24105')
if (! dir.exists('./14_TMB')){
  dir.create('./14_TMB')
}
setwd('./14_TMB')
library(R.utils)
library(tidyverse)
library(readxl)
library(writexl)
library(maftools)
# BiocManager::install("GenVisR")
# txt_files <- list.files(pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
if(F){
  files <- list.files(pattern = '*.gz',recursive = TRUE)
  all_mut <- data.frame()
  for (file in files) {
    mut <- read.delim(file,skip = 7, header = T, fill = TRUE,sep = "\t")
    all_mut <- rbind(all_mut,mut)
  }
  write.table(all_mut, file="input.maf", sep = "\t", quote = F, row.names = F)
}



#####
###??
# library(GenVisR)
library(maftools)
library(openxlsx)


################################################################################
######################## 交集候选基因TMB #######################################
################################################################################

rt <- read.table("input.maf",sep = "\t",header = T,check.names = F)
#outTab <- substr(rt$Tumor_Sample_Barcode, 1, 12)
#rt <- data.frame(id = outTab,rt)
gene <- read.table("../03_upset/feature_genes.txt",sep = "\t",header = F,check.names = F)
#rownames(gene) = gene$gene
data = merge.data.frame(rt, gene, by.x = "Hugo_Symbol", by.y = "V1")
write.table(data,"gene.maf", quote = F, sep = "\t", row.names = F,
            col.names = T)

laml <- read.maf("gene.maf")
#pdf("high_summary.pdf",width = 10)
#plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#dev.off()

pdf("waterfull.pdf", height = 6, width = 7)
oncoplot(maf = laml, titleText = "", fontSize = 0.75)
dev.off()

tiff(file = "waterfull.tiff", width = 7, height = 6, units = "in",
     compression = "lzw", bg = "white", res = 400)
oncoplot(maf = laml, titleText = "", fontSize = 0.75)
dev.off()


library(ggplot2)
library(ggbeeswarm)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(ggsci)
#BiocManager::install("ggbeeswarm")
laml <- read.maf("E:/LZ/24080/14_TMB/input.maf")
tmb_table_wt_log = tmb(maf = laml)
# laml <- read.maf("input.maf")
# getClinicalData(laml)
# mafSurvival(maf = laml, genes = 'DNMT3A', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)
colnames(tmb_table_wt_log)[1] = "id"
tmb_table_wt_log$id = substring(tmb_table_wt_log$id, 1, 12)
st <- read.table("../00_rawdata/risk.txt",sep = "\t",header = T)
st$id = substring(st$id, 1, 12)
#st <- cbind(id = st$id, risk = st$risk)
tmb <- merge(tmb_table_wt_log, st, by = "id")
tmb$risk <- factor(tmb$risk,levels = c("low","high"))
p <- ggviolin(tmb, x = "risk", y = "total_perMB_log", fill = "risk",
              palette = c("#4DBBD5FF", "#DC0000FF"), add = "boxplot",
              add.params = list(fill="white"))+
  ylab("Tumor Burden Mutation")+
  stat_compare_means(label.x = 1.25)  #根据图片修改
p
ggsave(filename = "TMB.tiff", plot = p, dpi = 400, width = 6, height = 5)
ggsave(filename = "TMB.pdf", plot = p, dpi = 500, width = 6, height = 5)




################################################################################
######################## 高低风险TMB #######################################
################################################################################

#TCGA-LIHC.somaticmutation_wxs.tsv
# rt <- all_mut
rt <- read.table("E:/LZ/24080/14_TMB/input.maf",sep = "\t",header = T,check.names = F)
outTab <- substr(rt$Tumor_Sample_Barcode, 1, 12)

rt <- data.frame(id = outTab,rt)
st <- read.table("../00_rawdata/risk.txt",sep = "\t",header = T)
st$id = substring(st$id, 1, 12)
st <- cbind(id = st$id, risk = st$risk)

risk_input <- merge(st, rt, by = "id")
risk_input$Tumor_Sample_Barcode <- risk_input$Matched_Norm_Sample_Barcode

num <- grep("low", risk_input$risk)

low_risk <- risk_input[num, ]

num <- grep("high",risk_input$risk)
high_risk <- risk_input[num,]

low_risk <- low_risk[,-c(1,2)]
high_risk <- high_risk[,-c(1,2)]

write.table(low_risk,"low_risk_input.maf", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(high_risk,"high_risk_input.maf", quote = F, sep = "\t", row.names = F, col.names = T)
###################high#################
setwd('E:/LZ/24105')
if (! dir.exists('./14_TMB')){
  dir.create('./14_TMB')
}
setwd('./14_TMB')
laml <- read.maf("high_risk_input.maf")
pdf("high_summary.pdf",width = 8)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
pdf("high_waterfull.pdf",height = 5)
oncoplot(maf = laml, top = 20)
dev.off()

tiff("high_summary.tiff", width = 8,  height = 7, units = "in", bg = "white",
     res = 500)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median',
               dashboard = TRUE, titvRaw = FALSE)
dev.off()

tiff("high_waterfull.tiff", width = 5,  height = 5, units = "in", bg = "white",
     res = 500)
oncoplot(maf = laml, top = 20)
dev.off()

# pdf("high_modelgene_waterfall.pdf")
# oncostrip(maf = laml, genes = gene$V1[-1])
# dev.off()
#################low#####################
laml2 <- read.maf("low_risk_input.maf")
pdf("low_summary.pdf",width = 8)
plotmafSummary(maf = laml2, rmOutlier = TRUE, addStat = 'median',
               dashboard = TRUE, titvRaw = FALSE)
dev.off()
pdf("low_waterfull.pdf",height = 5)
oncoplot(maf = laml2, top = 20)
dev.off()

tiff("low_summary.tiff", width = 8,  height = 7, units = "in", bg = "white",
     res = 500)
plotmafSummary(maf = laml2, rmOutlier = TRUE, addStat = 'median',
               dashboard = TRUE, titvRaw = FALSE)
dev.off()

tiff("low_waterfull.tiff", width = 5,  height = 5, units = "in", bg = "white",
     res = 500)
oncoplot(maf = laml2, top = 20)
dev.off()



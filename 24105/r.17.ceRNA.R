rm(list = ls())
setwd("G:/LZ/24080")
if (!dir.exists("18_ceRNA")) {dir.create("18_ceRNA")}
setwd("18_ceRNA")

#利用r包筛选mirna
library(multiMiR)
hubgene <- read.csv('../06_multicox/multiCox.csv')

gene2mir <- get_multimir(org     = 'hsa',
                         target  = hubgene$id,
                         table   = 'predicted',
                         summary = TRUE,
                         predicted.cutoff.type = 'n',
                         predicted.cutoff= 500000)
table(gene2mir@data$database)
# diana_microt        elmmo    microcosm      miranda        mirdb       pictar         pita 
#   142          286           205           255           180           207           262 
# targetscan 
#   133
targetscan = gene2mir@data[gene2mir@data$database=="targetscan",] # 133
targetscan$miRNA <- paste0(targetscan$mature_mirna_id,"_",targetscan$target_symbol)
targetscan2<- rstatix::filter(targetscan,!duplicated(targetscan$miRNA))
write.csv(targetscan2,"01.targetscan_mirna_filter.csv",row.names = F)  # 118

mirdb = gene2mir@data[gene2mir@data$database=="mirdb",] 
mirdb$miRNA <- paste0(mirdb$mature_mirna_id,"_",mirdb$target_symbol)
mirdb2<- rstatix::filter(mirdb,!duplicated(mirdb$miRNA))
write.csv(mirdb2,"01.mirdb_mirna_filter.csv",row.names = F)   # 180

hubgene <- data.frame(symbol=base::intersect(targetscan2$miRNA,mirdb2$miRNA))
mirna<- rbind(targetscan2,mirdb2)
mirna <- mirna[which(mirna$miRNA%in%hubgene$symbol),]
mirna2<- rstatix::filter(mirna,!duplicated(mirna$miRNA))
mirna2<-mirna2[,c("mature_mirna_id","target_symbol","miRNA")]
colnames(mirna2)<-c("mirna","mrna","mirna_mrna")
write.csv(mirna2,"02.mirna.csv",row.names = F)
# 42
mirna2$mirna
#[1] "hsa-miR-506-3p"   "hsa-miR-130a-5p"  "hsa-miR-129-2-3p" "hsa-miR-129-1-3p"

#[5] "hsa-miR-519d-3p" "hsa-miR-130b-3p" "hsa-miR-4295"    "hsa-miR-93-5p"  
#[9] "hsa-miR-454-3p"  "hsa-miR-17-5p"   "hsa-miR-20b-5p"  "hsa-miR-3666"   
#[13] "hsa-miR-106a-5p" "hsa-miR-106b-5p" "hsa-miR-20a-5p"  "hsa-miR-19b-3p" 
#[17] "hsa-miR-19a-3p"  "hsa-miR-148a-3p" "hsa-miR-152-3p"  "hsa-miR-148b-3p"


#lncRNA-------------------------------------------------------------------
# 读取mi-lnc互作的数据
if(F){
  rm(list = ls())
  setwd("E:/02-OA/11_ceRNA/lncRNA")
  # lncRNA_miRNA_interaction下载自ENCORI
  file_list <- list.files(pattern = "\\.txt$")
  combined_df <- do.call(rbind, lapply(file_list, function(file) {
    read.table(file, header=T,sep="\t",comment.char ="#",stringsAsFactors=F)
  }))
  
  # combined_df <- combined_df %>% filter(geneType == "lncRNA")
  write.csv(combined_df,"mi_lnc_interact.csv")
}


####----------------------------------------------------------------------------
rm(list = ls())
setwd("G:/LZ/24080")
if (!dir.exists("18_ceRNA")) {dir.create("18_ceRNA")}
setwd("18_ceRNA")



starbase <- read.csv('E:/04/06_ceRNA/mi_lnc_interact.csv', row.names = 1)
mirna<-read.csv('./02.mirna.csv')
lncRNA1 <- merge(starbase,mirna,by.x="miRNAname",by.y="mirna")
lncRNA1 <-lncRNA1[which(lncRNA1$clipExpNum>10),]
clinical1 <- lncRNA1 # 224
clinical1 <- clinical1[,c(21,1,4)]
colnames(clinical1) <- c("mRNA","miRNA","lncRNA")

clinical1$miRNA_lncRNA <- paste0(clinical1$miRNA,"_",clinical1$lncRNA)
clinical2<- rstatix::filter(clinical1,!duplicated(clinical1$miRNA_lncRNA))
write.csv(clinical2,"02.lncrna.csv",row.names = F)


#数据合并方便画图
mi<-read.csv('02.mirna.csv')
lnc<-read.csv('02.lncrna.csv')
mi<-mi[,2:1]
colnames(mi)<-c('val1','val2')
lnc<-lnc[,2:3]
colnames(lnc)<-c('val1','val2')
ce<-rbind(mi,lnc)
#检查有没有重复的lnc
{
ce2<-ce
ce2$val3 <- paste0(ce2$val1,"_",ce2$val2)
ce2<- rstatix::filter(ce2,!duplicated(ce2$val3))    #检查有没有重复的lnc
}
write.csv(ce,'03.cerna_cytoscape.csv',row.names = F)
# 190
ce2<- rstatix::filter(lnc,!duplicated(lnc$val2))
ce2$val2 # 20
# [1] "H19"        "AC021078.1" "SNHG16"     "MALAT1"     "LINC01618"  "XIST"      
# [7] "NEAT1"      "OIP5-AS1"   "MIR497HG"   "SNHG22"     "MIR17HG"    "THBS3-AS1" 
# [13] "SNHG7"      "DANT2"      "KCNQ1OT1"   "Z93241.1"   "AC079781.5" "AL356488.2"
# [19] "SNHG1"      "LINC00632" 




#########-------------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)

mi<-read.csv('02.mirna.csv')
lnc<-read.csv('02.lncrna.csv')
colnames(mi) <- c("miRNA","mRNA","mRNA_miRNA")
df <- left_join(mi,lnc,by="miRNA")
# 去除重复的miRNA_lncRNA
# df<- rstatix::filter(df,!duplicated(df$miRNA_lncRNA)) 
df <- df[,c(1,2,5,6)]
colnames(df) <- c("miRNA","mRNA","lncRNA","clipExpNum")

# df1 = make_long(df, miRNA, mRNA, lncRNA)

df1 = make_long(df, lncRNA, miRNA, mRNA)

mycol<- c4a('rainbow_wh_rd',68)
mycol2<- sample(mycol,length(mycol)) #随机打乱配色顺序
mycol2 <- scales::hue_pal()(length(unique(df$clipExpNum)))
custom_palette <- RColorBrewer::brewer.pal(n = length(unique(df$clipExpNum)), 
                                           name = "Set3")
#图表美化：
p1<- ggplot(df1, aes(x = x,
                    next_x= next_x,
                    node= node,
                    next_node= next_node,
                    fill= node,
                    label= node)) +
  scale_fill_manual(values = mycol) + #更改配色
  geom_sankey(flow.alpha = 0.5, #条带不透明度
              smooth= 8, #条带弯曲度
              width= 0.2 #节点宽度
              ) +  
  geom_sankey_text(size = 3.2,
                   color= "black") +
  #theme_void()+
  ggsci::scale_fill_aaas() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))
p1


#将条带变为灰色：
p2<- ggplot(df1, aes(x = x,
                    next_x= next_x,
                    node= node,
                    next_node= next_node,
                    fill= node,
                    label= node)) +
  geom_sankey(flow.alpha = 0.5,
              flow.fill = 'grey', #条带填充色
              flow.color = 'grey80', #条带描边色
              node.fill = mycol, #节点填充色
              smooth= 8,
              width= 0.2) +
  geom_sankey_text(size = 7,
                   color= "black")+
  #theme_void()+
  ggsci::scale_fill_aaas() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))

p2
ggsave(filename = "ceRNA.pdf", width = 20, height = 18, p2,dpi=400)
ggsave(filename = "ceRNA.png", width = 20, height = 18, p2,dpi=400)




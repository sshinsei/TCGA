rm(list=ls())
library(survminer)
library(survival)
setwd("E:/LZ/24105")  #modify work path
if(! dir.exists("07_COMBATmulticox")){
  dir.create("07_COMBATmulticox")
}
setwd("./07_COMBATmulticox")



library(tidyverse)
##----------要删除正常样本----------------##
norExp = read.table(file = "../00_rawdata/normalizeExp.txt", 
                    sep = "\t", header = T, check.names = F) # 430,有一列是id
cli_data = read.table(file = "../00_rawdata/survival.csv", sep = ",", header = T)
#fea_genes = read.table(file = "../04_clustDEG/diffSig.csv", sep = "\t", header = T, check.names = F)
#fea_genes <- fea_genes %>% rownames_to_column("V1")
fea_genes <- read.table(file = "../03_upset/feature_genes.txt",sep = ",", header = F, check.names = F) #28
#fea_genes <- read.table(file = "../00_rawdata/ERGs_2.txt", 
#                        sep = "\t", header = F, check.names = F)
colnames(fea_genes)[1] <- "V1"
rows = na.omit(match(fea_genes$V1, norExp$id)) # 332
norExp_1 = norExp[rows, ]
unique(norExp_1$id == fea_genes$V1)  #check
rownames(norExp_1) = norExp_1$id
norExp_1 = norExp_1[ , -1]
norExp_1 = t(norExp_1)
norExp_1 = as.data.frame(norExp_1)
norExp_1$id = rownames(norExp_1)
#rows = na.omit(match(rownames(norExp_1), cli_data$Id))
#norExp_2 = cbind(norExp_1, cli_data[rows, ])
#norExp_2 = norExp_2[ , c((ncol(norExp_2)-2):ncol(norExp_2), 1:(ncol(norExp_2)-3))]
#unique(rownames(norExp_2) == norExp_2$Id)  #check

norExp_1$id <- substring(norExp_1$id, 1, 16) # 424


norExp_2 = merge(cli_data, norExp_1, by.x = "sample", by.y = "id") # 418
# 删除正常样本
norExp_2 = norExp_2[!grepl("11A",norExp_2$sample),] # 367
# 去重
norExp_2 <- norExp_2 %>% 
  distinct(sample, .keep_all = TRUE) # 367
write.table(norExp_2, file = "../00_rawdata/tumor_clinical.txt", sep = "\t", row.names = F, quote = F)


rm(list=ls())
rt = read.table("../00_rawdata/tumor_clinical.txt", 
                header = T, sep = "\t", row.names = 1,
                check.names = FALSE) # 367
rt = na.omit(rt)
rt = rt[rt$futimes > 30, ] # 344
geneCoef = read.table("../06_lasso/geneCoef.txt",header = T, sep = "\t", check.names = FALSE) # 2
#geneCoef = read.csv('../06_rf/rf_genes.csv')
#geneCoef = read.table("./gene.txt",header = F, sep = "\t", check.names = FALSE) 
colnames(geneCoef)[1] <- "V1"
cols = match(geneCoef[ , 1], colnames(rt))
rt1 = rt[ , c(1, 2, cols)]
colnames(rt1)[3:ncol(rt1)] == geneCoef$V1  #check colnames

#rt = read.table("multiInput.txt", header = T, sep = "\t", check.names = F, row.names = 1)
#rt1[ , 3:ncol(rt1)] = log2(rt1[ , 3:ncol(rt1)] + 0.001)
rt1[ , "futimes"] = rt1[ , "futimes"]/365
cox <- coxph(Surv(futimes, fustate) ~ ., data = rt1)
cox = step(cox, direction = "both", steps = 5000)

# 根据Riscscore对患者进行高低风险分组
riskScore = predict(cox, type = "risk", newdata = rt1)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
write.table(cbind(id = rownames(rt), cbind(rt1[ , 1:2], riskScore, risk)),
            file = "../00_rawdata/risk.txt", sep = "\t", quote = F, row.names = F)
cox

#Output model parameters
multiCoxSum = summary(cox)
outTab = data.frame()
outTab = cbind(
  coef = multiCoxSum$coefficients[ , "coef"],
  HR = multiCoxSum$conf.int[ , "exp(coef)"],
  HR.95L = multiCoxSum$conf.int[ , "lower .95"],
  HR.95H = multiCoxSum$conf.int[ , "upper .95"],
  pvalue = multiCoxSum$coefficients[ , "Pr(>|z|)"])
outTab = cbind(id = row.names(outTab), outTab)
outTab = gsub("`", "", outTab)
write.table(outTab, file = "multiCox.csv", sep = ",", row.names = F, quote = F)

# 图要修改
if(F){
  pdf(file="forest.pdf", width = 10, height = 7, pointsize = 12)
  ggforest(cox,
           main = "Hazard ratio",
           cpositions = c(-4.5, -2.5, -0.8), 
           fontsize = 0.8, 
           refLabel = "reference", 
           noDigits = 2)
  dev.off()
  
}





if(F){
  # --------------------------foretploter-----------------------------------------------
  # devtools::install_github("adayim/forestploter")
  library(forestploter)
  rt <- read.csv(file = "multiCox.xls", sep = "\t",row.names = 1)
  str(rt)
  rt[, c(2, 3, 4)] <- round(rt[, c(2, 3, 4)], 5)
  rt$`HR (95%CI)` = paste0(rt$HR, " (", rt$HR.95L, "-", rt$HR.95H, ")", sep = "")
  rt$`HR (95%CI)`
  rt$id <- rownames(rt)
  # modify arrange of columns
  rt <- rt[,c(7,1:6)]
  # 添加样本数量
  rt$Number <- rep(344,nrow(rt))
  
  c_index <- multiCoxSum$concordance[1]
  # ------------------------------------plot-----------------------------------------------
  # 定义一个简单的主题，大家可以随意发挥自己的审美！
  tm <- forest_theme(base_size = 10,           # 设置基础字体大小
                     refline_col = "red4",     # 设置参考线颜色为红色
                     arrow_type = "closed",    # 设置箭头类型为闭合箭头
                     footnote_col = "black",   # 设置脚注文字颜色为蓝色
                     footnote_cex = 0.8          # 设置注脚文字大小
  )
  
  # 绘制森林图
  p <- forest(rt[,c(1,8,7,6)],   # 选择要在森林图中使用的数据列，这里包括变量名列、患者数量列、绘图要用的空白列和HR（95%CI）列
              est = rt$HR,          # 效应值，也就是HR列
              lower = rt$`HR.95L`,  # 置信区间下限
              upper = rt$`HR.95H`,  # 置信区间上限
              sizes = 0.5,        # 黑框框的大小
              ci_column = 3,             # 在第3列（可信区间列）绘制森林图
              ref_line = 1,              # 添加参考线
              arrow_lab = c("Low risk", "High Risk"),  # 箭头标签，用来表示效应方向，如何设置取决于你的样本情况
              xlim = c(0.995, 1.003),          # 设置x轴范围
              ticks_at = c(0.995, 0.997, 0.999, 1.001, 1.003),  # 在指定位置添加刻度
              theme = tm,                # 添加自定义主题
              footnote = paste("Concordance Index (C-index):", round(c_index, 3)))  # 添加脚注信息
  p
  
}

# ------------------------forestplot----------------------------------------------
library(forestplot)

rt <- read.csv(file = "multiCox.csv", row.names = 1)
str(rt)
rt[, c(2, 3, 4,5)] <- round(rt[, c(2, 3, 4,5)], 5)
rt$`HR (95%CI)` = paste0(rt$HR, " (", rt$HR.95L, "-", rt$HR.95H, ")", sep = "")
rt$`HR (95%CI)`
rt$id <- rownames(rt)
# modify arrange of columns
rt <- rt[,c(7,1:6)]
# 添加样本数量
rt$Number <- rep(344,nrow(rt))
# 新增一行
message1 <- c("Gene", "Coef", NA, NA, NA, "P.value","HR (95%CI)","Number")
rt1 <- rbind(message1,rt)
c_index <- multiCoxSum$concordance[1]
round(c_index, 3)
line <- paste0(nrow(rt)+2)

fig <- forestplot(
  rt1[, c(1,8,7,6)],                # 需要显示在森林图中的列
  mean = rt1$HR,             # 均值列（HR），它将显示为森林图的小方块或其他形状哈哈哈哈哈
  lower = rt1$`HR.95L`,            # 95%置信区间的下限数据列
  upper = rt1$`HR.95H`,            # 95%置信区间的上限数据列
  zero = 1,                        # 均值为1时的参考线，也就是零线
  boxsize = 0.15,                   # 方框的大小
  graph.pos = 3,             # 森林图在右侧显示
  hrzl_lines = list(               # 水平线样式的设置
    "1" = gpar(lty = 1, lwd = 2),  # 均值线
    "2" = gpar(lty = 2),           # 下限和上限之间的虚线
    "4" = gpar(lwd = 2, lty = 1) # 下限和上限线
  ),
  graphwidth = unit(.25, "npc"),   # 森林图的宽度
  xlab = paste0("Concordance Index (C-index): ",round(c_index, 3)), # x轴标签
  xticks = c(0.7,0.8, 0.9,  1.0, 1.1, 1.2,1.3,1.4), # x轴刻度
  # 判断是否为汇总行，汇总行就是连续变量或者分类变量名所在的行，可以加粗让它显眼一点，好看的！
  # is.summary = c(T, F, T, F, F, T, F, F, T, F, F, F),  
  txt_gp = fpTxtGp(                # 文本样式的设置
    label = gpar(cex = 1),       # 标签的大小
    ticks = gpar(cex = 0.8),         # 刻度标记的大小
    xlab = gpar(cex = 1),        # x轴标签的大小
    title = gpar(cex = 1.2)        # 标题的大小
  ),
  lwd.zero = 1,                    # 均值线的宽度
  lwd.ci = 1.5,                    # 置信区间线的宽度
  lwd.xaxis = 2,                   # x轴线的宽度
  lty.ci = 1.5,                    # 置信区间线的样式
  ci.vertices = T,                 # 是否显示置信区间的顶点
  ci.vertices.height = 0.2,        # 置信区间顶点的高度
  clip = c(0.1, 8),                # x轴的截断范围
  ineheight = unit(8, 'mm'),       # 森林图内线的高度
  line.margin = unit(8, 'mm'),     # 森林图线的边距
  colgap = unit(6, 'mm'),          # 森林图中方框的间距
  fn.ci_norm = "fpDrawDiamondCI",  # 置信区间的形状（这里是钻石形）
  title = "MultivariateCox Forest Plot",   # 森林图的标题
  col = fpColors(                  # 颜色设置
    box = "#e31a1c",                 # 方框颜色 #e31a1c 绿色：#76bf08
    lines = "#e31a1c",               # 线条颜色
    zero = "black"                 # 均值为0时的颜色
  )
)
fig

tiff(file="forest2.tiff",width=8,height=6,units ="in",compression="lzw",
     bg="white",res=400)
fig
dev.off()


pdf(file="forest2.pdf", width = 8, height = 6)
fig
dev.off()



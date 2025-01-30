##机器学习筛选特征基因
rm(list = ls())
setwd('E:/LZ/24105') 
if(!dir.exists("20_model")){dir.create("20_model")}
setwd("20_model")

library(readr)
library(magrittr)
library(caret)
# 加载表达矩阵数据
exp <- read.csv(file ="../00_rawdata/GSE104580/02normalize/GSE104580.txt",
                sep="\t",row.names=1)
# 加载分组变量
group_list <- read.csv(file = "../00_rawdata/GSE104580/02normalize/GSE104580.group.txt",
                       sep="\t")
Candidate_genes <- read.csv('../07_multicox/multiCox.csv', header = T) # 54
# Candidate_genes <- read.csv('../04_DEG_clust/DEG_sig.csv', header = T) # 1461

# 提取样本子集和特征基因
expr_region <- exp %>% subset(., select = group_list$sample)
expr_top <- data.frame(t(expr_region[Candidate_genes$id, ]))
expr_top$Group <- group_list$group
expr_top$Group <- factor(expr_top$Group, levels = c('Normal', 'Case'))

# 数据划分为训练集和验证集
set.seed(123)
train_index <- createDataPartition(expr_top$Group, p = 0.7, list = FALSE)
train_data <- expr_top[train_index, ]
test_data <- expr_top[-train_index, ]

# 评估指标函数
if(F){
	evaluate_model <- function(predictions, true_labels) {
		conf_matrix <- confusionMatrix(as.factor(predictions), as.factor(true_labels))
		precision <- conf_matrix$byClass["Precision"]
		recall <- conf_matrix$byClass["Recall"]
		f1_score <- 2 * (precision * recall) / (precision + recall)
		return(list(Precision = precision, Recall = recall, F1_Score = f1_score))
	}
}



# new metrics function
# 加载必要的包
library(pROC)

evaluate_model <- function(predictions, true_labels, probabilities = NULL) {
  conf_matrix <- confusionMatrix(as.factor(predictions), as.factor(true_labels))
  precision <- conf_matrix$byClass["Precision"]
  recall <- conf_matrix$byClass["Recall"]
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  # 计算 AUC 值
  auc_value <- NA
  if (!is.null(probabilities)) {
    roc_obj <- roc(as.numeric(true_labels) - 1, as.numeric(probabilities))
    auc_value <- auc(roc_obj)
  }else{
    roc_obj <- roc(as.numeric(true_labels) - 1, as.numeric(predictions))
    auc_value <- auc(roc_obj)
  }
  
  return(list(Precision = precision, Recall = recall, F1_Score = f1_score, AUC = auc_value))
}




###----------------------------LASSO------------------------------------------------####

library(glmnet)
set.seed(11123123)         ##set.seed(121434342)

#-------------10折交叉验证LASSO---------------------------------#
# LASSO 模型训练
res.en <- cv.glmnet(as.matrix(train_data[-ncol(train_data)]), train_data$Group, 
                    family = "binomial", nfolds = 10, alpha = 1)



# 在验证集上进行预测并评估
lasso_probabilities <- predict(res.en, 
                               newx = as.matrix(test_data[-ncol(test_data)]), 
                               s = res.en$lambda.min, type = "response")
lasso_predictions <- ifelse(lasso_probabilities > 0.5, "Case", "Normal")
lasso_metrics <- evaluate_model(lasso_predictions, test_data$Group, lasso_probabilities)
lasso_metrics


coefficient <- coef(res.en, s=res.en$lambda.min)
coefficient
Active.Index <- which(as.numeric(coefficient)!=0)
astive.coefficients <- as.numeric(coefficient)[Active.Index]


##
res.en$lambda.1se # 0.0113915
res.en$lambda.min # 0.004493043
# 找出那些回归系数没有被惩罚为0的
png("LASSO.CV.png", width = 7, height = 7,units = "in", bg = "white",res=600,family='Times')
par(pin = c(4,4), mar = c(6,6,6,1))  
plot(res.en, 
     #cex.axis=1.8,  ##。cex。
     cex.lab=1.6,   ##。cex。
     cex.main=2.0,   ##。
     #font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)
dev.off()

pdf("LASSO.CV.pdf", width = 7, height = 7,family='Times')
par(pin = c(4,4), mar = c(6,6,6,1))  
plot(res.en, 
     #cex.axis=1.8,  ##。cex。
     cex.lab=1.6,   ##。cex。
     cex.main=2.0,   ##。
     #font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)
dev.off()


png("LASSO.Coef.png", width = 7, height = 7,units = "in", bg = "white",res=600,family='Times')
par(pin = c(4,4), mar = c(6,6,6,1))  
plot(res.en$glmnet.fit, xvar = 'lambda',
     #cex.axis=1.8,  ##。cex。
     cex.lab=1.6,   ##。cex。
     cex.main=2.0,   ##。
     #font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)
dev.off()

pdf("LASSO.Coef.pdf", width = 7, height = 7,family='Times')
par(pin = c(4,4), mar = c(6,6,6,1))  
plot(res.en$glmnet.fit, xvar = 'lambda',
     cex.axis=1.8,  ##。cex。
     cex.lab=2.0,   ##。cex。
     cex.main=2.0,   ##。
     #font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)
dev.off()

lasso_geneids <- rownames(coefficient)[Active.Index] 
lasso_geneids<-lasso_geneids[-1]   
lasso_geneids # 5
# "SLC7A11" "FAP"     "COL1A1"  "ALOX5"   "LYN"    
write.table(lasso_geneids,'lasso_geneids.txt',sep='\t',quote=F,row.names=F,col.names=F)





#---------------------------SVM-REF --------------------------------------------#
set.seed(123)
library(mlbench)
library(caret)
#data(PimaIndiansDiabetes)
#class(PimaIndiansDiabetes[,9])
#class(dat1[,1:159])
train_data_rfe <- train_data
train_data_rfe$Group <- as.factor(train_data_rfe$Group)
control <- rfeControl(functions = caretFuncs, method = "cv", number = 10)

svm_results <- rfe(train_data_rfe[, -ncol(train_data_rfe)], 
                   train_data_rfe$Group, 
                   sizes = c(1:(ncol(train_data_rfe) - 1)),
                   rfeControl = control,
                   method = "svmRadial")
#saveRDS(svm_results,"SVM_RFE.rds")
# 运行一次后直接读取
#library(mlbench)
#library(caret)
#results <- readRDS("SVM_RFE.rds")

# 在验证集上预测
# 提取 RFE 最优模型的特征子集
optimal_features <- predictors(svm_results)
test_data_rfe <- test_data[, optimal_features]

# 使用训练好的 SVM 模型预测
svm_model <- svm_results$fit
svm_predictions <- predict(svm_model, test_data_rfe)
svm_probabilities <- NULL
svm_metrics <- evaluate_model(svm_predictions, test_data$Group, svm_probabilities)
svm_metrics


svm_geneids <- data.frame(id = predictors(svm_results))
svm_geneids <- svm_geneids[[1]]
svm_geneids
# FAP
write.table(svm_geneids,"svm_geneids.txt",sep = "\t",row.names = F,col.names = T,quote = F)


#------------------------只选择n个特征时预测效果最好-------------------------------------#
pdf("SVM_RFE.pdf",w = 5, h = 5,family='Times')
plot(svm_results, type="o",cex.lab = 3,cex.axis = 4,main='SVM-RFE')
dev.off()
png("SVM_RFE.png",w = 5, h = 5,units = "in",res = 600,bg='white')
plot(svm_results, type="o",cex.lab = 3,cex.axis = 4)
dev.off()






#-------------------- RF -----------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------#
library(tidyverse)
library(caret)
# library(DALEX)
library(randomForest)

# 随机森林模型训练
set.seed(1002342)
train_control <- trainControl(method = "cv", number = 10)
rf_model <- train(Group ~ ., data = train_data,
                  method = "rf",
                  ntree = 400,
                  trControl = train_control)

# 在验证集上预测并评估
rf_probabilities <- predict(rf_model, test_data, type = "prob")[, "Case"]
rf_predictions <- ifelse(rf_probabilities > 0.5, "Case", "Normal")
rf_metrics <- evaluate_model(rf_predictions, test_data$Group, rf_probabilities)


# 最终模型
rf = rf_model$finalModel
# 获取模型的混淆矩阵
confusionMatrix(rf_model)

# 查看模型摘要信息
print(rf_model)

# 查看模型的交叉验证结果
print(rf_model$results)




# 累积残差分布图
png("cumulative.png", width = 6.5, height = 6, bg = "white", units = "in", res = 600,family='Times')
par(pin = c(4,4), mar = c(6,6,6,1))
plot(rf)
dev.off()

pdf("cumulative.pdf", width = 6.5, height = 6,family='Times')
par(pin = c(4,4), mar = c(6,6,6,1))
plot(rf)
dev.off()


pdf('Feature_Importance.pdf',w=6.5,h=7,family='Times')
# 绘制特征重要性图
varImpPlot(rf)
dev.off()

png('Feature_Importance.png',w=6.5,h=7,bg = "white", units = "in", res = 600,family='Times')
# 绘制特征重要性图
varImpPlot(rf)
dev.off()


#-----------------------取特征重要性top n----------------------------------------#
# 查看特征重要性

importance <- as.data.frame(importance(rf)) 
importance <- importance %>% rownames_to_column("gene")

# 按特征重要性降序排列
importance_rf1 <- importance %>% arrange(desc(importance$MeanDecreaseGini))
importance_rf1

top10 <- head(importance_rf1, 5)
# top10 <- importance_rf1[importance_rf1$MeanDecreaseGini >1,]
# importance_rf1 <- data.frame(importance(rf))
rf_geneids <- top10$gene
rf_geneids
# [1] "RCC2"     "FAP"      "TNFRSF21" "SLC7A11"  "COL1A1"   "CSK"      "CD300A"  
# [8] "ADAM12"   "NINJ2"    "ALOX5"  

write.table(rf_geneids,'RF_geneids.txt',sep='\t',quote=F,row.names=F,col.names=F)




#----------------------------------整合模型评价----------------------------------#
# 汇总评估结果
eval_metrics <- data.frame(
  Model = c("LASSO", "SVM-RFE", "RF"),
  Precision = c(lasso_metrics$Precision, svm_metrics$Precision, rf_metrics$Precision),
  Recall = c(lasso_metrics$Recall, svm_metrics$Recall, rf_metrics$Recall),
  F1_Score = c(lasso_metrics$F1_Score, svm_metrics$F1_Score, rf_metrics$F1_Score),
  AUC = c(lasso_metrics$AUC, svm_metrics$AUC, rf_metrics$AUC)
)

# 打印评估结果
print(eval_metrics)
write.table(eval_metrics,"metrics.csv")

library(reshape2)
# 可视化评估结果
eval_metrics_long <- melt(eval_metrics, id.vars = "Model")
p <- ggplot(eval_metrics_long, aes(x = Model, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "", x = "model", y = "value", fill = "") +
  scale_fill_manual(values = c('#eba58d', '#f2c0c1', '#a0b9e0','#e0e2ed'))+###修改颜色
  theme_classic()+
  theme(legend.position='top',
        panel.grid = element_blank(),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size = 20,color ='black'),
        axis.text.y = element_text(size=20,color ='black'),#坐标轴数字黑色加粗
        legend.text = element_text(size = 16),
        axis.line = element_line(colour = "black",linewidth =1.0),#坐标轴1.0磅，黑色
        axis.ticks.length.x.bottom = unit (0.2, "cm"),#修改坐标轴刻度线长度
        axis.ticks.length.y.left = unit(0.2,'cm'),
        axis.ticks = element_line(colour = "black",linewidth =1.0)
  )+
  scale_y_continuous(expand = c(0, 0))
p
ggsave("score.pdf", p, width = 8, height = 6)
ggsave("score.png", p, width = 8, height = 6)

save.image("model.RData")





rm(list=ls())
setwd('~/project/24064') 
if(!dir.exists("05_model")){dir.create("05_model")}
setwd("05_model")
load("model.RData")

if(F){
  ###-----------------LASSO---SVM-----rf取交集-----------------------------------------###
  
  lasso_geneids
  # "SLC7A11" "FAP"     "COL1A1"  "ALOX5"   "LYN" 
  rf_geneids
  # "RCC2"     "FAP"      "TNFRSF21" "SLC7A11"  "COL1A1" 
  boruta_geneids
  svm_geneids; # "ARL4C"
  xgb_geneids <- xgb_importance$Feature;
  xgb_geneids
  rf_geneids #10
  
  
  vennlist <- list(`LASSO` = lasso_geneids, `SVM` = svm_geneids,
                   `RF` = rf_geneids)
  
  inter_ID <- Reduce(intersect,vennlist) # 2
  paste0(inter_ID,collapse = ', ')
  #  "TREM2, ALOX5, HAVCR2"
  saveRDS(inter_ID,"lasso_svm_xgb_gene.txt")
  write.table(as.data.frame(inter_ID),'lasso_svm_xgb_gene.txt',
              row.names = FALSE,col.names=FALSE, sep = '\t', quote = FALSE)
  
  library(VennDiagram)
  library(RColorBrewer)
  mycolor <-brewer.pal(12, 'Paired')[c(1, 3, 5,7)]
  library(ggvenn)
  opar <- par(family = "Roboto Condensed")
  
  pdf(file = paste0("lasso_svm_xgb.pdf"),width = 8,height = 6)
  a <- dev.cur()
  png(file = paste0("lasso_svm_xgb.png"),width= 8, height= 6, units="in", res=300)
  dev.control("enable")
  par(mar = c(2,2,2,2),cex=1.5,family="Times")
  ggvenn(vennlist,fill_color=mycolor,fill_alpha = .4,
         stroke_linetype = "longdash",set_name_size = 5,
         stroke_color = "transparent",
         show_percentage = T,
         text_size= 4)
  # mycolor <- c("#FF4040", "#00BFFF")
  dev.copy(which = a)
  dev.off()
  dev.off()
  
  
  
  
  
  
  ###-----------------LASSO--RF取交集-----------------------------------------###
  
  
  
  vennlist <- list(`LASSO` = lasso_geneids,
                   `RF` = rf_geneids)
  
  inter_ID <- Reduce(intersect,vennlist)
  paste0(inter_ID,collapse = ', ') # 4
  # SLC7A11, FAP, COL1A1
  saveRDS(inter_ID,"11.lasso_rf.txt")
  write.table(as.data.frame(inter_ID),'11.lasso_rf.txt',
              row.names = FALSE,col.names=FALSE, sep = '\t', quote = FALSE)
  
  library(VennDiagram)
  library(RColorBrewer)
  mycolor <-brewer.pal(12, 'Paired')[c(1, 3, 5,7)]
  library(ggvenn)
  opar <- par(family = "Roboto Condensed")
  
  pdf(file = paste0("11.lasso_rf.pdf"),width = 8,height = 6)
  a <- dev.cur()
  png(file = paste0("11.lasso_rf.png"),width= 8, height= 6, units="in", res=300)
  dev.control("enable")
  par(mar = c(2,2,2,2),cex=1.5,family="Times")
  ggvenn(vennlist,fill_color=mycolor,fill_alpha = .4,
         stroke_linetype = "longdash",set_name_size = 5,
         stroke_color = "transparent",
         show_percentage = T,
         text_size= 4)
  # mycolor <- c("#FF4040", "#00BFFF")
  dev.copy(which = a)
  dev.off()
  dev.off()
  
  
  
  
  
  
  ###-----------------SVMRFE-RF取交集-----------------------------------------###
  
  lasso_geneids;
  svm_geneids; 
  rf_geneids # top5
  
  vennlist <- list(`SVM` = svm_geneids,
                   `RF` = rf_geneids)
  inter_ID <- Reduce(intersect,vennlist)
  paste0(inter_ID,collapse = '，') # 8
  saveRDS(inter_ID,"12.svm_rf.txt")
  write.table(as.data.frame(inter_ID),'12.svm_rf.txt',
              row.names = FALSE,col.names=FALSE, sep = '\t', quote = FALSE)
  
  library(VennDiagram)
  library(RColorBrewer)
  mycolor <-brewer.pal(12, 'Paired')[c(1, 3, 5,7)]
  library(ggvenn)
  opar <- par(family = "Roboto Condensed")
  
  pdf(file = paste0("12.svm_rf.pdf"),width = 8,height = 6)
  a <- dev.cur()
  png(file = paste0("12.svm_rf.png"),width= 8, height= 6, units="in", res=300)
  dev.control("enable")
  par(mar = c(2,2,2,2),cex=1.5,family="Times")
  ggvenn(vennlist,fill_color=mycolor,fill_alpha = .4,
         stroke_linetype = "longdash",set_name_size = 5,
         stroke_color = "transparent",
         show_percentage = T,
         text_size= 4)
  # mycolor <- c("#FF4040", "#00BFFF")
  dev.copy(which = a)
  dev.off()
  dev.off()
  
  
  
  
  ###-----------------SVM-LASSO取交集-----------------------------------------###
  
  lasso_geneids;
  svm_geneids; 
  # rf_geneids
  
  vennlist <- list(`SVM-RFE` = svm_geneids,
                   `LASSO` = lasso_geneids)
  inter_ID <- Reduce(intersect,vennlist)
  paste0(inter_ID,collapse = '，')
  # "KLF7，FAM162A，FOXO3，EGFR，HK1"
  saveRDS(inter_ID,"13.svm_lasso.txt")
  write.table(as.data.frame(inter_ID),'13.svm_lasso.txt',
              row.names = FALSE,col.names=FALSE, sep = '\t', quote = FALSE)
  
  library(VennDiagram)
  library(RColorBrewer)
  mycolor <-brewer.pal(12, 'Paired')[c(1, 3, 5,7)]
  library(ggvenn)
  opar <- par(family = "Roboto Condensed")
  
  pdf(file = paste0("13.svm_lasso.pdf"),width = 8,height = 6)
  a <- dev.cur()
  png(file = paste0("13.svm_lasso.png"),width= 8, height= 6, units="in", res=300)
  dev.control("enable")
  par(mar = c(2,2,2,2),cex=1.5,family="Times")
  ggvenn(vennlist,fill_color=mycolor,fill_alpha = .4,
         stroke_linetype = "longdash",set_name_size = 5,
         stroke_color = "transparent",
         show_percentage = T,
         text_size= 4)
  # mycolor <- c("#FF4040", "#00BFFF")
  dev.copy(which = a)
  dev.off()
  dev.off()
  
  
  
  ###-----------------SVM-XGB取交集-----------------------------------------###
  xgb_geneids <- xgb_importance$Feature;
  xgb_geneids; # 21
  svm_geneids; 
  # rf_geneids
  
  vennlist <- list(`SVM-RFE` = svm_geneids,
                   `XGB` = xgb_geneids)
  inter_ID <- Reduce(intersect,vennlist) # 19
  paste0(inter_ID,collapse = '，')
  saveRDS(inter_ID,"14.svm_xgb.txt")
  write.table(as.data.frame(inter_ID),'14.svm_xgb.txt',
              row.names = FALSE,col.names=FALSE, sep = '\t', quote = FALSE)
  
  library(VennDiagram)
  library(RColorBrewer)
  mycolor <-brewer.pal(12, 'Paired')[c(1, 3, 5,7)]
  library(ggvenn)
  opar <- par(family = "Roboto Condensed")
  
  pdf(file = paste0("14.svm_lasso.pdf"),width = 8,height = 6)
  a <- dev.cur()
  png(file = paste0("14.svm_lasso.png"),width= 8, height= 6, units="in", res=300)
  dev.control("enable")
  par(mar = c(2,2,2,2),cex=1.5,family="Times")
  ggvenn(vennlist,fill_color=mycolor,fill_alpha = .4,
         stroke_linetype = "longdash",set_name_size = 5,
         stroke_color = "transparent",
         show_percentage = T,
         text_size= 4)
  # mycolor <- c("#FF4040", "#00BFFF")
  dev.copy(which = a)
  dev.off()
  dev.off()
  
  
  
  ###-----------------LASSO--XGB取交集-----------------------------------------###
  
  
  lasso_geneids;
  xgb_geneids
  
  
  vennlist <- list(`LASSO` = lasso_geneids,
                   `XGB` = xgb_geneids)
  
  inter_ID <- Reduce(intersect,vennlist)
  paste0(inter_ID,collapse = '，') # 11
  saveRDS(inter_ID,"15.lasso_xgb.txt")
  write.table(as.data.frame(inter_ID),'15.lasso_xgb.txt',
              row.names = FALSE,col.names=FALSE, sep = '\t', quote = FALSE)
  
  library(VennDiagram)
  library(RColorBrewer)
  mycolor <-brewer.pal(12, 'Paired')[c(1, 3, 5,7)]
  library(ggvenn)
  opar <- par(family = "Roboto Condensed")
  
  pdf(file = paste0("15.lasso_xgb.pdf"),width = 8,height = 6)
  a <- dev.cur()
  png(file = paste0("15.lasso_xgb.png"),width= 8, height= 6, units="in", res=300)
  dev.control("enable")
  par(mar = c(2,2,2,2),cex=1.5,family="Times")
  ggvenn(vennlist,fill_color=mycolor,fill_alpha = .4,
         stroke_linetype = "longdash",set_name_size = 5,
         stroke_color = "transparent",
         show_percentage = T,
         text_size= 4)
  # mycolor <- c("#FF4040", "#00BFFF")
  dev.copy(which = a)
  dev.off()
  dev.off()
}



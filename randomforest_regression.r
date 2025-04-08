#清空所有的环境变量
rm(list = ls())
#加载R包
library(randomForest)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggstatsplot)
library(patchwork)
#读取实验设计表
tc_map = read.table("design.txt",header = T, row.names = 1,
                    sep = '\t')
dim(tc_map)
class(tc_map)
tc_map$Group = as.factor(tc_map$Group)
class(tc_map)
class(tc_map$Group)
#读取domain矩阵
otu_table =read.csv("pfam.txt",sep = '\t',header = T, row.names = 1)
#筛选Group1组数据作为训练集数据
#筛选符合条件的行（group1当作训练集）
sub_map = tc_map[tc_map$Group %in% c("group1"),]
dim(sub_map)
sub_map
#####一共筛选了？？个数据作为训练集数据
class(sub_map)
# 筛选OTU
rownames(sub_map)
###逻辑索引，在索引数据
idx = rownames(sub_map) %in% colnames(otu_table)
idx
sub_map = sub_map[idx,]
sub_map
dim(sub_map)
sub_otu = otu_table[, rownames(sub_map)] 
dim(sub_otu)

## 随机森林回归
set.seed(315)
rf = randomForest(t(sub_otu), sub_map$Topt, importance=TRUE, 
                  proximity=TRUE,ntree = 1000)
print(rf)
####这里的参数有时间要好好琢磨一下


## 导出feature重要性
imp= as.data.frame(rf$importance)
class(imp)
head(imp)
###这里的两个参数也要好好研究一下是什么意思，以便更好的理解随机森林模型
###排序贡献度
imp = imp[order(imp[,1],decreasing = T),]
head(imp,n=10)
write.table(imp,file = "importance_class.txt",quote = F,sep = '\t',
            row.names = T, col.names = F)

## 模型预测效率
# 模型评估
train.p = predict(rf, type = "response")
train.p
df = data.frame(observed = sub_map$Topt, predict = train.p)
head(df)
# 保存预测结果与真实结果比较
write.table(df,file = "all_train_predict.txt",quote = F,sep = '\t',
            row.names = T, col.names = NA)
# spearman or pearson(统计语言的相关系数自己选择)
df
cor = cor.test(df[,1], df[,2], method = "spearman") 
cor$p.value
df <- as.data.frame(df)
m = lm(observed ~ predict, df)
m
p <- ggplot(df, aes(x=predict, y=observed)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = paste("rho = " , round(cor$estimate, digits = 3),
                     ", P = " , signif(cor$p.value, digits = 3),
                     ", R2 = ", round(summary(m)$r.squared, digits = 3),
                     sep = "")) +
  theme_bw()
ggsave(paste("train_predict.pdf", sep=""), p, width = 8, height = 5)






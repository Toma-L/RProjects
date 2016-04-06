#資料挖礦與大數據分析==========

library(MASS)
library(RSNNS)
data("Pima.tr")
data("Pima.tr")
set.seed(1111)

Pima <- rbind(Pima.tr, Pima.te)
level_name <- {}

for(i in 1:7) { #貝氏僅支援類別變數，要將連續型屬性離散化
        Pima[, i] = cut(Pima[, i], breaks = 2, ordered_result = TRUE, include.lowest = TRUE)
        level_name <- rbind(level_name, levels(Pima[, i]))
}

level_name <- data.frame(level_name)
row.names(level_name) <- colnames(Pima)[1:7]
colnames(level_name) <- paste("L", 1:2, sep = "")
level_name
Pima.tr <- Pima[1:200, ]
Pima.te <- Pima[201:nrow(Pima), ]


install.packages("bnlearn")
library(bnlearn)

bn <- naive.bayes(Pima.tr, "type") #naive.bayes建立簡單貝氏分類法
plot(bn)
bn
pred <- predict(bn, Pima.te)
tab <- table(pred, Pima.te[, "type"])
tab
acc <- sum(diag(tab)) / sum(tab)
acc


tan <- tree.bayes(Pima.tr, "type") #tree.bayes建構貝氏網路
plot(tan)
tan
fitted <- bn.fit(tan, Pima.tr, method = "bayes")
pred <- predict(fitted, Pima.te)
tab <- table(pred, Pima.te[, "type"])
tab
acc <- sum(diag(tab)) / sum(tab)
acc
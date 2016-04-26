#利用R語言打通大數據的經脈==================================================

#08判別分析==================================================

##費希爾判別
#選出適當的投影軸，使類別內離差盡可能小，不同類別間離差盡可能大

#線性判別分析LDA：將樣本點投影到一維空間，實際常常更複雜
#二次判別分析QDA：投影至二維空間，使用許多二次曲面，常用的非線性判別函數


##貝氏判別：利用貝氏選出最大後驗機率的類別
#不怕雜訊和無關變數，但假設各特徵屬性之間是無關的（現實常常並非如此）


#距離判別：根據待判斷樣本與已知類別樣本之間的距離遠近做判別
#對類別域的交換或重疊較多的待分樣本集來說，此方法較為適合
#常用的有k最近鄰法（k-Nearest Neighbor; kNN）、有權數的k最近鄰（Weighted k-Nearest Neighbor）


##資料準備

library(kknn)
data(miete)
head(miete)
dim(miete)
summary(miete)

library(sampling)
n <- round(2 / 3 * nrow(miete) / 5) #訓練集佔資料總量2/3，計算每一等級中應取出的樣本數
n
sub_train <- strata(miete, stratanames = "nmkat", size = rep(n, 5), method = "srswor") #分層抽樣
head(sub_train)

data_train <- getdata(miete[, c(-1, -3, -12)], sub_train$ID_unit)
data_test <- getdata(miete[, -1, -3, -12], -sub_train$ID_unit)
dim(data_train)
dim(data_test)
head(data_test)


##線性判別分析==================================================

install.packages("MASS")
library(MASS)
fit_lda1 <- lda(nmkat ~., data_train) #執行線性判別
names(fit_lda1) #輸出項目名稱

fit_lda1$prior #先驗機率
fit_lda1$counts #各種類別的樣本數

fit_lda1$means #各變數在每一種類中的平均值

fit_lda1


fit_lda2 <- lda(data_train[, -12], data_train[, 12])
fit_lda2

plot(fit_lda1)

plot(fit_lda1, dimen = 1) #輸出1個判別式的圖形
plot(fit_lda1, dimen = 2)


pre_lda1 <- predict(fit_lda1, data_test) #預測
pre_lda1$class #各樣本的預測結果
pre_lda1$posterior #每一種類別的後驗機率

table(data_test$nmkat, pre_lda1$class) #混淆矩陣

error_lda1 <- sum(as.numeric(as.numeric(pre_lda1$class) != as.numeric(data_test$nmkat)))/nrow(data_test)
error_lda1


##單純貝氏分類==================================================

install.packages("klaR")
library(klaR)
fit_Bayes1 <- NaiveBayes(nmkat ~., data_train)
names(fit_Bayes1)
fit_Bayes1$apriori
fit_Bayes1$tables #用於建立判別規則的所有變數在各種類別下的條件機率
fit_Bayes1$levels
fit_Bayes1$call
fit_Bayes1$usekernel
fit_Bayes$varnames

plot(fit_Bayes1, vars = "wfl", n = 50, col = c(1, "darkgrey", 1, "darkgrey", 1))
plot(fit_Bayes1, vars = "mvdauer", n = 50, col = c(1, "darkgrey", 1, "darkgrey", 1))
plot(fit_Bayes1, vars = "nmqm", n = 50, col = c(1, "darkgrey", 1, "darkgrey", 1))

fit_Bayes1 <- NaiveBayes(data_train[, -12], data_train[, 12])
pre_Bayes1 <- predict(fit_Bayes1, data_test)
pre_Bayes1

table(data_test$nmkat, pre_Bayes1$class)
error_Bayes1 <- sum(as.numeric(as.numeric(pre_Bayes1$class) != as.numeric(data_test$nmkat)))/nrow(data_test)
error_Bayes1


##K最近鄰==================================================

install.packages("class")
library(class)
fit_pre_knn <- knn(data_tain[, -12], data_test[, -12], cl = data_train[, 12])
fit_pre_knn
table(data_test$nmkat, fit_pre_knn)
error_knn <- sum(as.numeric(as.numeric(fit_pre_knn) != as.numeric(data_test$nmkat)))/nrow(data_test)
error_knn

error_knn <- rep(0, 20)
for(i in 1:20){
        fit_pre_knn <- knn(data_train[, -12], data_test[, -12], cl = data_train[, 12], k = i)
        error_knn[i] <- sum(as.numeric(as.numeric(fit_pre_knn) != as.numeric(data_test$nmkat)))/nrow(data_test)
}
error_knn

plot(error_knn, type = "l", xlab = "K")


##有權數的K最近鄰演算法==================================================

install.packages("kknn")
library(kknn)
fit_pre_kknn <- kknn(nmkat ~., data_train, data_test[, -12], k = 5)
summary(fit_pre_kknn)
fit <- fitted(fit_pre_kknn)
fit

table(data_test$nmkat, fit)
error_kknn <- sum(as.numeric(as.numeric(fit) != as.numeric(data_test$nmkat))) / nrow(data_test)
error_kknn
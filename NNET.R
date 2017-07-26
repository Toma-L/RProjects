#資料挖礦與大數據分析==================================================

##倒傳遞類神經網路==================================================

library(MASS)
# install.packages("RSNNS")
library(RSNNS)
data(Pima.tr)
set.seed(1111)
Pima.tr <- Pima.tr[sample(1:nrow(Pima.tr), length(1:nrow(Pima.tr))), ] #將資料順序重新排列
PimaValues <- Pima.tr[, 1:7]
PimaTargets <- decodeClassLabels(Pima.tr[, 8]) #目標屬性重新編碼，因為目標屬性為二分類變數
Pima.tr <- splitForTrainingAndTest(PimaValues, PimaTargets, ratio = .1) #資料切割
Pima.tr <- normTrainingAndTestSet(Pima.tr) #避免7個屬性不同尺度影響分析結果，進行標準化

model <- mlp(Pima.tr$inputsTrain, Pima.tr$targetsTrain, size = 14, learnFuncParams = .01, maxit = 100, inputsTest = Pima.tr$imputsTest, targetsTest = Pima.tr$targetsTest)
plotIterativeError(model)
weightMatrix(model)


p_table <- expand.grid(size = c(12, 13, 14, 15, 16), learning.rate = c(0.001, 0.01, 0.1))
for(i in 1:nrow(p_table)) {
        model <- mlp(Pima.tr$inputsTrain, Pima.tr$targetsTrain, size = p_table[i, 1], 
                     learnFuncParams = p_table[i, 2], 
                     maxit = 100, inputsTest = Pima.tr$inputsTest, targetsTest = Pima.tr$targetsTest)
        p_table$TestError[i] = model$IterativeTestError[100]
}
p_table

Pima.te[, 1:7] <- normalizeData(Pima.te[, 1:7])
predictions <- predict(model, Pima.te[, 1:7])
table <- confusionMatrix(Pima.te[, 8], predictions)
accuracy <- sum(diag(table)) / sum(table)
accuracy


##自我組織映射網路==================================================

library(MASS)
install.packages("kohonen")
library(kohonen)
data("Pima.tr")
Pima_class <- rbind(Pima.tr, Pima.te)[, 8]
Pima <- scale(rbind(Pima.tr, Pima.te)[, -8]) #將屬性標準化以避免不同尺度影響分群結果

set.seed(1111)
Pima_som <- som(data = Pima, grid = somgrid(4, 4, "hexagonal"), rlen = 1000, alpha = c(.05, .01))
plot(Pima_som, type = "changes")
plot(Pima_som, type = "dist.neighbours")
plot(Pima_som, type = "codes")
plot(Pima_som, type = "counts")


##適應性共振理論類神經網路==================================================

library(RSNNS)
data(snnsData)
patterns <- snnsData$art1_letters.pat
inputMaps <- matrixToActMapList(patterns, nrow = 7)
par(mfrow = c(3, 3))
for(i in 1:9) plotActMap(inputMaps[[i]])

model <- art1(patterns, dimX = 7, dimY = 5, learnFuncParams = c(.5, 0, 0), maxit = 100)
table(encodeClassLabels(model$fitted.values))


#利用R語言打通大數據的經脈==================================================

#13神經網路==================================================

install.packages("nnet")
library(nnet)

vector1 <- c("a", "b", "a", "c")
vector2 <- c(1, 2, 1, 3)
class.ind(vector1) #對模型中的y進行前置處理
class.ind(vector2)

url <- "http://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv"
wine <- read.csv(url, header = TRUE, sep = ";")

cha = 0
for(i in 1:4898) {
        if(wine[i, 12] > 6) cha[i] = "good"
        else if(wine[i, 12] > 5) cha[i] = "mid"
        else cha[i] = "bad"
}
wine[, 12] = factor(cha)
summary(wine$quality)

scale01 <- function(x) { #資料歸一化，NNet常見的前置處理，將所有資料轉化為[0, 1]之間的數，取消個維度數據間數量級的差別
        ncol = dim(x)[2] - 1
        nrow = dim(x)[1]
        new = matrix(0, nrow, ncol)
        for(i in 1:ncol) {
                max = max(x[, i])
                min = min(x[, i])
                for(j in 1:nrow) {
                        new[j, i] = (x[j, i] - min) / (max - min)
                }
        }
        new
}

names(wine)
set.seed(71)
samp <- sample(1:4898, 3000)
wine[samp, 1:11] <- scale01(wine[samp, ])
r = 1/max(abs(wine[samp, 1:11])) #參數rang的變化範圍
set.seed(101)
model1 <- nnet(quality ~., data = wine, subset = samp, #rang為初始隨機權重的範圍
               size = 4, rang = r, decay = 5e-4, maxit = 200)

x <- subset(wine, select = -quality)
y <- wine[, 12]
y <- class.ind(y)
set.seed(101)
model2 <- nnet(x, y, decay = 5e-4, maxit = 200, size = 4, rang = r) #decay為權重衰減精度

summary(model1)
summary(model2)

x <- wine[, 1:11]
pred <- predict(model1, x, type = "class") #用於預測的變數個數要跟建立模型時一樣！
set.seed(110)
pred[sample(1:4898, 8)]


xt <- wine[, 1:11]
pred <- predict(model2, xt) #第二種模型
dim(pred)
pred[sample(1:4898, 4), ] #呈現3個輸出結果的值，表示樣本為某種類別的機率

name <- c("bad", "good", "mid")
prednew <- max.col(pred) #max.col()確定每列中最大值的那行
prednewn <- name[prednew]
set.seed(201)
prednewn[sample(1:4898, 8)]


true <- max.col(y)
table(true, prednewn) #檢查預測精度


data(iris)

x <- iris[, -5]
y <- iris[, 5]

x <- subset(iris, select = -Species)
y <- class.ind(y)

model1 <- nnet(x, y, rang = 1/max(abs(x)), size = 4, maxit = 500, decay = 5e-4)
model2 <- nnet(x, y, rang = 1/max(abs(x)), size = 4, maxit = 500, decay = 5e-4)

#每次建模型使用的反覆運算初值都是不同的

model1$convergence #結果為0表示反覆運算會停止「並非」因為達到最大反覆運算次數
model2$convergence #所以最大反覆運算次數並非造成兩模型不同的主因

model1$value #最後值為模型擬合標準同模型權數衰減值的和，越小表示擬合效果越好
model2$value

name <- c("setosa", "versicolor", "virginica")
pred1 <- name[max.col(predict(model1, x))]
pred2 <- name[max.col(predict(model2, x))]
table(iris$Species, pred1) #展示模型精度
table(iris$Species, pred2)

#要得到最佳模型可以多嘗試不同的模型、測試每一節點數目下模型的誤判率


wine <- read.table("wine.txt", sep = ";")
names(wine) <- c("fixed", "volatile", "citric", "residual", "chlorides", 
                 "free", "total", "density", "PH", "sulphates", "alcohol", "quality")
set.seed(71)
wine <- wine[sample(1:4898, 3000), ]
nrow.wine <- dim(wine)[1]
scale01 <- function(x) { #歸一化程式
        ncol = dim(x)[2] - 1
        nrow = dim(x)[1]
        new <- matrix(0, nrow, ncol)
        for(i in 1:ncol) {
                max = max(x[, i])
                min = min(x[, i])
                for(j in 1:nrow) new[j, i] = (x[j, i] - min)/(max - min)
        }
        new
}

cha = 0
for(i in 1:nrow.wine) {
        if(wine[i, 12] > 6) cha[i] = "good"
        else if(wine[i, 12] > 5) cha[i] = "mid"
        else cha[i] = "bad"
}
wine[, 12] <- factor(cha)

set.seed(444)
samp <- sample(1:nrow.wine, nrow.wine * .7)
wine[samp, 1:11] <- scale01(wine[samp, ])
wine[-samp, 1:11] <- scale01(wine[-samp, ])
r = 1/max(abs(wine[samp, 1:11]))
n <- length(samp)

#嘗試不同隱藏層節點個數

err1 = 0
err2 = 0
for(i in 1:17) {
        set.seed(111)
        model <- nnet(quality ~., data = wine, maxit = 400, rang = r, size = i, subset = samp, decay = 5e-4)
        err1[i] <- sum(predict(model, wine[samp, 1:11], type = "class") != wine[samp, 12]) / n
        err2[i] <- sum(predict(model, wine[-samp, 1:11], type = "class") != wine[-samp, 12]) / (nrow.wine - n)
}

par(family = "Songti TC Light")
plot(1:17, err1, "l", col = 1, lty = 1, ylab = "模型誤判率", xlab = "隱藏層節點個數", 
     ylim = c(min(min(err1), min(err2)), max(max(err1), max(err2))))
lines(1:17, err2, col = 1, lty = 3)
points(1:17, err1, col = 1, pch = "+")
points(1:17, err2, col = 1, pch = "o")
legend(1, 0.53, "測試集誤判率", bty = "n", cex = 1.0) #後面有過度擬合的情形
legend(1, 0.35, "訓練集誤判率", bty = "n", cex = 1.0) #隱藏層節點3個就夠了

#測試不同訓練週期

errl1 <- 0
errl2 <- 0
for(i in 1:500) {
        set.seed(111)
        model <- nnet(quality ~., data = wine, maxit = i, rang = r, size = 3, subset = samp)
        errl1[i] <- sum(predict(model, wine[samp, 1:11], type = "class") != wine[samp, 12]) / n
        errl2[i] <- sum(predict(model, wine[-samp, 1:11], type = "class") != wine[-samp, 12]) / (nrow.wine - n)
}
plot(1:length(errl1), errl1, "l", ylab = "模型誤判率", xlab = "訓練週期", col = 1, 
     ylim = c(min(min(errl1), min(errl2)), max(max(errl1), max(errl2))))
lines(1:length(errl1), errl2, col = 1, lty = 3)
legend(250, .47, "測試集誤判率", bty = "n", cex = .8)
legend(250, .425, "訓練集誤判率", bty = "n", cex = .8) #模型的誤判率會趨於平穩

set.seed(111)
model <- nnet(quality ~., data = wine, maxit = 300, rang = r, size = 3, subset = samp)
x <- wine[-samp, 1:11]
pred <- predict(model, x, type = "class")
table(wine[-samp, 12], pred)


#應用R語言於資料分析

#監督式學習 6.3 人工神經網路==================================================

install.packages("neuralnet")
library(neuralnet)

traininginput <- as.data.frame(runif(100, min = 0, max = 100))
trainingoutput <- sqrt(traininginput)

trainingdata <- cbind(traininginput, trainingoutput)
colnames(trainingdata) <- c("Input", "Output")

net.sqrt <- neuralnet(Output ~ Input, trainingdata, algorithm = "backprop", hidden = 10, threshold = .01, learningrate = .01, )

print(net.sqrt)
plot(net.sqrt)

testdata <- as.data.frame((1:10)^2) #產生測試資料
net.results <- compute(net.sqrt, testdata) #用compute()來預測模型結果

cleanoutput <- cbind(testdata, sqrt(testdata), as.data.frame(net.results$net.result))
colnames(cleanoutput) <- c("Input", "Expected Output", "Neural Net Output")
print(cleanoutput)


install.packages("DMwR")
library(DMwR)
regr.eval(cleanoutput[, 'Expected Output'],  #求平均絕對誤差MAE、均方根誤差RMSE
          cleanoutput[, 'Neural Net Output'], stats = c("mae", "rmse"))
#資料挖礦與大數據分析==========

#CART

library(MASS)
library(rpart) #CART決策樹
data("Pima.tr")
summary(Pima.tr)
set.seed(1111)
cart <- rpart(type ~., Pima.tr, control = rpart.control(cp = 0)) #cp是複雜係數alpha
summary(cart)
par(xpd = TRUE)
plot(cart)
text(cart)

cart_prune <- prune(cart, cp = .03)
par(xpd = TRUE)
plot(cart_prune)
text(cart)

pre <- predict(cart, Pima.te, type = "class")
confusion_matrix <- table(Type = Pima.te$type, Predict = pre)
confusion_matrix
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
accuracy


#C5.0

install.packages("C50")
library(C50) #C5.0決策樹
library(MASS)
data("Pima.tr")
C50_tree <- C5.0(type ~., Pima.tr, control = C5.0Control(noGlobalPruning = T)) #F則會修剪樹
summary(C50_tree)
plot(C50_tree)


#CHAID

install.packages("CHAID", repos="http://R-Forge.R-project.org", type="source")
library(CHAID)
data("Pima.tr")
data("Pima.te")
Pima <- rbind(Pima.tr, Pima.te)
level_name <- {} #將資料離散化，CHAID只能做離散
for(i in 1:7){
        Pima[, i] = cut(Pima[, i], breaks = 3, ordered_result = T, include.lowest = T)
        level_name <- rbind(level_name, levels(Pima[, i]))
}
level_name <- data.frame(level_name)
rownames(level_name) <- colnames(Pima[1:7])
colnames(level_name) <- paste("L", 1:3, sep = "")
level_name

Pima.tr <- Pima[1:200, ]
Pima.te <- Pima[201:nrow(Pima), ]
set.seed(1111)
CHAID_tree <- chaid(type ~., Pima.tr)
CHAID_tree
plot(CHAID_tree)


#利用R語言打通大數據的經脈==========

##資料準備==========

rpart(formula, data, weights, subset, na.action = na.rpart, method, model = FALSE, x = FALSE, y = TRUE, parms, control, cost, ...) #method參數用於選擇決策樹的類型，包含anova、poisson、class、exp四種，不進行設定的預設情況下，R會自行猜測，例如當y為因數型變數時，預設為class型，其中，anova對應的是回歸樹，class則是分類樹
rpart.control(minsplit = 20, minbucket = round(minsplit/3), cp = .01, maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, xval = 10, surrogatestyle = 0, maxdepth = 30) #minsplit為每個節點最少需要幾個樣本，minbucket為每個葉節點最少需要幾個樣本，cp為複雜度參數，maxdepth控制樹的高度


library(mvpart)
data(car.test.frame)
head(car.test.frame)
str(car.test.frame)
summary(car.test.frame)

Group_Mileage <- matrix(0, 60, 1)
Group_Mileage[which(car.test.frame$Mileage >= 26)] = "A" #連續變數轉為類別變數
Group_Mileage[which(car.test.frame$Mileage <= 22)] = "C"
Group_Mileage[which(Group_Mileage == 0)] = "B"

car.test.frame$Group_Mileage = Group_Mileage
car.test.frame[1:10, c(4, 9)]

a <- round(1/4 * sum(car.test.frame$Group_Mileage == "A"))
b <- round(1/4 * sum(car.test.frame$Group_Mileage == "B"))
c <- round(1/4 * sum(car.test.frame$Group_Mileage == "C"))
a; b; c

library(sampling)
sub <- strata(car.test.frame, stratanames = "Group_Mileage", 
              size = c(c, b, a), method = "srswor") #用strata()進行分層抽樣
sub

Train_Car <- car.test.frame[-sub$ID_unit, ]
Test_Car <- car.test.frame[sub$ID_unit, ]
nrow(Train_Car)
nrow(Test_Car)


##CART應用==========

###建立迴歸樹

library(rpart)
formula_Car_Reg <- Mileage ~ Price + Country + Reliability + Type + Weight + Disp. + HP
rp_Car_Reg <- rpart(formula_Car_Reg, Train_Car, method = "anova")
print(rp_Car_Reg)

printcp(rp_Car_Reg) #匯出回歸樹的CP表格用printcp()

summary(rp_Car_Reg)

rp_Car_Reg1 <- rpart(formula_Car_Reg, Train_Car, method = "anova", minsplit = 10) #anova為迴歸樹
print(rp_Car_Reg1)
printcp(rp_Car_Reg1)

rp_Car_Reg2 <- rpart(formula_Car_Reg, Train_Car, method = "anova", cp = .1)
print(rp_Car_Reg2)
printcp(rp_Car_Reg2)

rp_Car_Reg3 <- prune.rpart(rp_Car_Reg, cp = .1) #用prune.rpart達到相同較果
print(rp_Car_Reg3)
printcp(rp_Car_Reg3)

rp_Car_Reg4 <- rpart(formula_Car_Reg, Train_Car, method = "anova", maxdepth = 1) #也可用maxdepth控制樹的大小（事前修剪）
print(rp_Car_Reg4)
printcp(rp_Car_Reg4)


install.packages("rpart.plot") #繪製決策樹
library(rpart.plot)
rp_Car_Plot <- rpart(formula_Car_Reg, Train_Car, method = "anova", minsplit = 10)
print(rp_Car_Plot)
rpart.plot(rp_Car_Plot)

rpart.plot(rp_Car_Plot, type = 4) #更改決策樹類型
rpart.plot(rp_Car_Plot, type = 4, branch = 1) #branch
rpart.plot(rp_Car_Plot, type = 4, fallen.leaves = TRUE) #fallen.leaves

install.packages("maptree")
library(maptree)
draw.tree(rp_Car_Plot, col = rep(1, 7), nodeinfo = TRUE) #另一種決策樹製圖工具

plot(rp_Car_Plot, uniform = TRUE, main = "plot: Regression Tree") #plot()也可以畫決策樹
text(rp_Car_Plot, use.n = TRUE, all = TRUE)

post(rp_Car_Plot, file = "") #post()也可以畫決策樹


###建立分類樹

formula_Car_Cla <- Group_Mileage ~ Price + Country + Reliability + Type + Weight + Disp. + HP
rp_Car_Cla <- rpart(formula_Car_Cla, Train_Car, method = "class", minsplit = 5)
print(rp_Car_Cla)

rpart.plot(rp_Car_Cla, type = 4, fallen.leaves = TRUE) #繪製分類樹


###進行預測

pre_Car_Cla <- predict(rp_Car_Cla, Test_Car, type = "class")
pre_Car_Cla
p <- sum(as.numeric(pre_Car_Cla != Test_Car$Group_Mileage)) / nrow(Test_Car)
p
table(Test_Car$Group_Mileage, pre_Car_Cla)


##C4.5應用==========

install.packages("RWeka")
library(RWeka)
names(Train_Car) <- c("Price", "Country", "Reliability", "Mileage", "Type", "Weight", "Disp.", "HP", "Oil_Consumption")
formula <- Oil_Consumption ~ Price + Country + Reliability + Type + Weight + Disp. + HP
Train_Car$Oil_Consumption <- as.factor(Train_Car$Oil_Consumption) #C5.0的輸出變數資料型態必須是因素，若不是可用factor()函數來轉換
C45_0 <- J48(formula, Train_Car)
C45_0

summary(C45_0)

C45_1 <- J48(formula, Train_Car, control = Weka_control(M = 3))
C45_1
summary(C45_1)

plot(C45_1) #可能要先裝partykit套件


#實用R語言於資料分析==========

library(rpart) #CART
data(iris)
np <- ceiling(.1 * nrow(iris)) #抽10%當作測試資料
np
test.index <- sample(1:nrow(iris), np)
iris.testdata <- iris[test.index, ]
iris.traindata <- iris[-test.index, ]

iris.tree <- rpart(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, method = "class", data = iris.traindata)
iris.tree
summary(iris.tree)

plot(iris.tree)
text(iris.tree)

species.traindata <- iris$Species[-test.index]
train.predict <- factor(predict(iris.tree, iris.traindata, type = 'class'), levels = levels(species.traindata))
table.traindata <- table(species.traindata, train.predict)
table.traindata
correct.traindata <- sum(diag(table.traindata)) / sum(table.traindata) * 100
correct.traindata #訓練資料正確率

species.testdata <- iris$Species[test.index]
test.predict <- factor(predict(iris.tree, iris.testdata, type = 'class'), levels = levels(species.testdata))
table.testdata <- table(species.testdata, test.predict)
table.testdata
correct.testdata <- sum(diag(table.testdata)) / sum(table.testdata) * 100
correct.testdata


library(rpart)
data(iris)
np <- ceiling(.1 * nrow(iris))
np
test.index <- sample(1:nrow(iris), np)
iris.testdata <- iris[test.index, ]
iris.traindata <- iris[-test.index, ]

iris.tree <- rpart(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, method = "class", data = iris.traindata, control = rpart.control(minsplit = 5, cp = .0001, maxdepth = 30))
species.traindata <- iris$Species[-test.index]
train.predict <- factor(predict(iris.tree, iris.traindata, type = "class"), levels = levels(species.traindata))

table.traindata <- table(species.traindata, train.predict)
table.traindata
correct.traindata <- sum(diag(table.traindata)) / sum(table.traindata) * 100
correct.traindata

species.testdata <- iris$Species[test.index]
test.predict <- factor(predict(iris.tree, iris.testdata, type = "class"), levels = levels(species.testdata))
table.testdata <- table(species.testdata, test.predict)
table.testdata
correct.testdata <- sum(diag(table.testdata)) / sum(table.testdata) * 100
correct.testdata


library(C50) #載入C50套件
data(iris)
np <- ceiling(.1 * nrow(iris))
np

test.index <- sample(1:nrow(iris), np)
iris.test <- iris[test.index, ]
iris.train <- iris[-test.index, ]


c <- C5.0Control(subset = FALSE, bands = 0, winnow = FALSE, noGlobalPruning = FALSE, CF = .25, minCases = 2, fuzzyThreshold = FALSE, sample = 0, seed = sample.int(4096, size = 1) - 1L, earlyStopping = TRUE)
iris_treeModel <- C5.0(x = iris.train[, -5], y = iris.train$Species, control = c)
summary(iris_treeModel)

test.output <- predict(iris_treeModel, iris.test[, -5], type = "class")
n <- length(test.output)
number = 0
for(i in 1:n) {
        if(test.output[i] == iris.test[i, 5]) {
                number = number + 1
        }
}
test.accuracy = number / n * 100
test.accuracy

iris.train$Species = factor(iris.train$Species) #C5.0的輸出變數資料型態必須是因素，若不是可用factor()函數來轉換


library(C50)
library(stringr)
data(iris)
c <- C5.0Control(subset = FALSE, bands = 0, winnow = FALSE, noGlobalPruning = FALSE, CF = .25, minCases = 2, fuzzyThreshold = FALSE, sample = .9, seed = sample.int(4096, size = 1) - 1L, earlyStopping = TRUE, label = "Species")
iris_treeModel <- C5.0(x = iris[, -5], y = iris$Species, control = c)
summary(iris_treeModel)

x <- str_locate_all(iris_treeModel$output, "%)") #用str_locate_all()取得iris_treeModel這個array的output元素中%的位置
y <- substr(iris_treeModel$output, x[[1]][2]-4, x[[1]][2]-1) #用substr()取得%前1~前4位置的文字
test.error <- as.numeric(y) #用as.numeric()轉換成數字
test.correct <- 100 - test.error #扣掉錯誤率得到正確率
test.correct
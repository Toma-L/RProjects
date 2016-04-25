#R軟體資料分析基礎與應用==================================================

#19正規化和壓縮方法==================================================

acs <- read.table("http://jaredlander.com/data/acs_ny.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
head(acs)
testFrame <- data.frame(First = sample(1:10, 20, replace = TRUE),
                        Second = sample(1:20, 20, replace = TRUE),
                        Third = sample(1:10, 20, replace = TRUE),
                        Fourth = factor(rep(c("Alice", "Bob", "Charlie", "David"),
                                            5)),
                        Fifth = ordered(rep(c("Edward", "Frank", "Georgia", "Hank", "Isaac"), 4)),
                        Sixth = rep(c("a", "b"), 10), stringsAsFactors = FALSE)
head(testFrame)
head(model.matrix(First ~ Second + Fourth + Fifth, testFrame)) #model.matrix()

#Fifth是ordered factor，level間有大小關係
#Fourth的level少了一個

require(useful)
head(build.x(First ~ Second + Fourth + Fifth, testFrame, contrasts = FALSE)) #build.x()
head(build.x(First ~ Second + Fourth + Fifth, testFrame, 
             contrasts = c(Fourth = FALSE, Fifth = TRUE))) #只對Fourth使用所有level

acs$Income <- with(acs, FamilyIncome >= 150000)
head(acs)
acsX <- build.x(Income ~ NumBedrooms + NumChildren + NumPeople + NumRooms + NumUnits + NumVehicles + NumWorkers + OwnRent + YearBuilt + ElectricBill + FoodStamp + HeatingFuel + Insurance + Language - 1, 
                data = acs, contrasts = FALSE) #建立預測函數矩陣
class(acsX)
dim(acsX)
topleft(acsX, c = 6)
topright(acsX, c = 6)
acsY <- build.y(Income ~ NumBedrooms + NumChildren + NumPeople + NumRooms + 
                        NumUnits + NumVehicles + NumWorkers + OwnRent + YearBuilt + 
                        ElectricBill + FoodStamp + HeatingFuel + Insurance + Language - 1, data = acs) #建立反應變數
head(acsY)
tail(acsY)


require(glmnet)
set.seed(1863561)
acsCV1 <- cv.glmnet(x = acsX, y = acsY, 
                    family = "binomial", nfold = 5) #cv.glmnet()可以自動交叉驗證的值，預設alpha為1（LASSO）

acsCV1$lambda.min
acsCV1$lambda.1se
plot(acsCV1)

coef(acsCV1, s = "lambda.1se") #點代表沒被選中的變數，LASSO會把高度相關的變數排除掉
plot(acsCV1$glmnet.fit, xvar = "lambda")
abline(v = log(c(acsCV1$lambda.min, acsCV1$lambda.1se)), lty = 2)


set.seed(71623)
acsCV2 <- cv.glmnet(x = acsX, y = acsY, 
                    family = "binomial", nfold = 5, alpha = 0) #alpha = 0建立脊迴歸模型
acsCV2$lambda.min
acsCV2$lambda.1se
coef(acsCV2, s = "lambda.1se") #每個變數都會被保留，只是會被壓縮接近0
plot(acsCV2)
plot(acsCV2$glmnet.fit, xvar = "lambda")
abline(v = log(c(acsCV2$lambda.min, acsCV2$lambda.1se)), lty = 2)

set.seed(2834673)
theFolds <- sample(rep(x = 1:5, length.out = nrow(acsX))) #建立層別，要觀測值每次執行都落在同一層
alphas <- seq(from = .5, to = 1, by = .05) #尋找最佳的alpha值要加一層交叉驗證，只考慮>0.5的alpha，因為傾向LASSO好過Ridge Reg.
set.seed(5127151)
cl <- makeCluster(2) #啟動叢集

library(doParallel) #進行平行化運算
registerDoParallel(cl)
before <- Sys.time()
acsDouble <- foreach(i = 1:length(alphas), 
                     .errorhandling = "remove", #若發生錯誤，該迭代跳過
                     .inorder = FALSE, #整合的先後次序不重要
                     .multicombine = TRUE, #能夠同時接受好幾個引數
                     .export = c("acsX", "acsY", "alphas", "theFolds"), #將幾個變數通過.export載入foreach environment
                     .packages = "glmnet") %dopar% #每個worker都載入glmnet，%dopar%讓foreach以平行運算的方式執行
{
        print(alphas[i])
        cv.glmnet(x = acsX, y = acsY, family = "binomial", nfolds = 5,
                  foldid = theFolds, alpha = alphas[i])
}
after <- Sys.time()
stopCluster(cl) #確保完成後將叢集終止
after - before
sapply(acsDouble, class) #檢測該list是一個有11個cv.glmnet的物件列表

extractGlmnetInfo <- function(object) {
        lambdaMin <- object$lambda.min 
        lambda1se <- object$lambda.1se #找出備選中的lambda
        
        whichMin <- which(object$lambda == lambdaMin) #找出lambda落在路徑何處
        which1se <- which(object$lambda == lambda1se)
        data.frame(lambda.min = lambdaMin, error.min = object$cvm[whichMin],
                   lambda.1se = lambda1se, error1se = object$cvm[which1se]) #建立data.frame含有備選中的lambda和相關錯誤訊息
}

alphaInfo <- Reduce(rbind, lapply(acsDouble, extractGlmnetInfo)) #整合到一個data.frame
alphaInfo2 <- plyr::ldply(acsDouble, extractGlmnetInfo) #也可以用ldply
identical(alphaInfo, alphaInfo2)

alphaInfo$Alpha <- alphas
alphaInfo

require(reshape2)
require(stringr)
alphaMelt <- melt(alphaInfo, id.vars = "Alpha", value.name = "Value", variable.name = "Measure")
alphaMelt$Type <- str_extract(string = alphaMelt$Measure, pattern = "(min)|(1se)")
alphaMelt$Measure <- str_replace(string = alphaMelt$Measure, pattern = "//,(min|1se)", replacement = "")
alphaCast <- dcast(alphaMelt, Alpha + Type ~ Measure, value.var = "Value")
together <- function (x) {
        for(i in 1:dim(x)[1]) {
                if(i %% 2 == 1) {
                        x[i, 3] = x[i, 4]
                }
        }
        print(x)
}
alphaCast <- together(alphaCast)
alphaCast <- alphaCast[, -4]
names(alphaCast)[3] <- "error"

ggplot(alphaCast, aes(x = Alpha, y = error)) + geom_line(aes(group = Type)) + facet_wrap(~Type, scales = "free_y", ncol = 1) + geom_point(aes(size = lambda))


#Exploratory Data Analysis==================================================

#Principal Components Analysis and Singular Value Decomposition

set.seed(12345)
par(mar = rep(.2, 4))
dataMatrix <- matrix(rnorm(400), nrow = 40)
image(1:10, 1:40, t(dataMatrix)[, nrow(dataMatrix):1])
par(mar = rep(.2, 4))
heatmap(dataMatrix)

set.seed(678910)
for(i in 1:40) {
        #flip a coin
        coinFlip <- rbinom(1, size = 1, prob = .5)
        #if coin is heads add a common pattern to that row
        if(coinFlip) {
                dataMatrix[i, ] <- dataMatrix[i, ] + rep(c(0, 3), each = 5)
        }
}
par(mar = rep(.2, 4))
image(1:10, 1:40, t(dataMatrix)[, nrow(dataMatrix):1]) #增加了pattern

par(mar = rep(.2, 4))
heatmap(dataMatrix)

hh <- hclust(dist(dataMatrix))
dataMatrixOrdered <- dataMatrix[hh$order, ]
par(mfrow = c(1, 3))
image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1])
plot(rowMeans(dataMatrixOrdered), 40:1, , xlab = "Row Mean", ylab = "Row", pch = 19)
plot(colMeans(dataMatrixOrdered), xlab = "Column", ylab = "Column Mean", pch = 19)


svd1 <- svd(scale(dataMatrixOrdered))
par(mfrow = c(1, 3))
image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1])
plot(svd1$u[, 1], 40:1, , xlab = "Row", ylab = "First left singular vector", pch = 19)
plot(svd1$v[, 1], xlab = "Column", ylab = "First right singular vector", pch = 19)

par(mfrow = c(1, 2))
plot(svd1$d, xlab = "Column", ylab = "Singular value", pch = 19)
plot(svd1$d^2/sum(svd1$d^2), xlab = "Column", ylab = "Prop. of variance explained", pch = 19)

svd1 <- svd(scale(dataMatrixOrdered))
pca1 <- prcomp(dataMatrixOrdered, scale = TRUE)
plot(pca1$rotation[, 1], svd1$v[, 1], pch = 19, xlab = "Principal Component 1", ylab = "Rigth Sigular Vector 1")
abline(c(0, 1))

constantMatrix <- dataMatrixOrdered * 0
for(i in 1:dim(dataMatrixOrdered)[1]) {constantMatrix[i, ] <- rep(c(0, 1), each = 5)}
svd1 <- svd(constantMatrix)
par(mfrow = c(1, 3))
image(t(constantMatrix)[, nrow(constantMatrix):1])
plot(svd1$d, xlab = "Column", ylab = "Singular value", pch = 19)
plot(svd1$d^2/sum(svd1$d^2), xlab = "Column", ylab = "Prop. of variance explained", pch = 19)

set.seed(678910)
for(i in 1:40){
        coinFlip1 <- rbinom(1, size = 1, prob = .5)
        coinFlip2 <- rbinom(1, size = 1, prob = .5)
        if(coinFlip1) {
                dataMatrix[i, ] <- dataMatrix[i, ] + rep(c(0, 5), each = 5)
        }
        if(coinFlip2) {
                dataMatrix[i, ] <- dataMatrix[i, ] + rep(c(0, 5), 5)
        }
}
hh <- hclust(dist(dataMatrix))
dataMatrixOrdered <- dataMatrix[hh$order, ]

svd2 <- svd(scale(dataMatrixOrdered))
par(mfrow = c(1, 3))
image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1])
plot(rep(c(0, 1), each = 5), pch = 19, xlab = "Column", ylab = "Pattern 1")
plot(rep(c(0, 1), 5), pch = 19, xlab = "Column", ylab = "Pattern 2")

svd2 <- svd(scale(dataMatrixOrdered))
par(mfrow = c(1, 3))
image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1])
plot(svd2$v[, 1], pch = 19, xlab = "Column", ylab = "First right sigular vector")
plot(svd2$v[, 2], pch = 19, xlab = "Column", ylab = "Second right sigular vector")

svd1 <- svd(scale(dataMatrixOrdered))
par(mfrow = c(1, 2))
plot(svd1$d, xlab = "Column", ylab = "Singular value", pch = 19)
plot(svd1$d^2/sum(svd1$d^2), xlab = "Column", ylab = "Percent of variance explained", pch = 19)

dataMatrix2 <- dataMatrixOrdered
dataMatrix2[sample(1:100, size = 40, replace = FALSE)] <- NA
svd1 <- svd(scale(dataMatrix2)) #ERROR


install.packages("impute")
library(impute)
dataMatrix2 <- dataMatrixOrdered
dataMatrix2[sample(1:100, size = 40, replace = FALSE)] <- NA
dataMatrix2 <- impute.knn(dataMatrix2)$data
svd1 <- svd(scale(dataMatrixOrdered))
svd2 <- svd(scale(dataMatrix2))
par(mfrow = c(1, 2))
plot(svd1$v[, 1], pch = 19)
plot(svd2$v[, 1], pch = 19)


load("/Users/thomas/Downloads/face.rda")
image(t(faceData)[, nrow(faceData):1])
svd1 <- svd(scale(faceData))
plot(svd1$d^2/sum(svd1$d^2), pch = 19, xlab = "Singular vector", ylab = "Variance explained")

svd1 <- svd(scale(faceData))
approx1 <- svd1$u[, 1] %*% t(svd1$v[, 1]) * svd1$d[1] #%*%為矩陣相乘
approx5 <- svd1$u[, 1:5] %*% diag(svd1$d[1:5]) %*% t(svd1$v[, 1:5])
approx10 <- svd1$u[, 1:10] %*% diag(svd1$d[1:10]) %*% t(svd1$v[, 1:10])

quartz()
par(mfrow = c(1, 4))
image(t(approx1)[, nrow(approx1):1], main = "(a)")
image(t(approx5)[, nrow(approx5):1], main = "(b)")
image(t(approx10)[, nrow(approx10):1], main = "(c)")
image(t(faceData)[, nrow(faceData):1], main = "(d)")


#Practical Machine Learning==================================================

#Preprocessing with Principal Components Analysis (PCA)==================================================

library(caret)
library(kernlab)
data(spam)

inTrain <- createDataPartition(y = spam$type, p = .75, list = FALSE)
training <- spam[inTrain, ]
testing <- spam[-inTrain, ]
M <- abs(cor(training[, -58]))
diag(M) <- 0
which(M > .8, arr.ind = TRUE) #挑出cor > .8的

names(spam)[c(34, 32)]
plot(spam[, 34], spam[, 32])

X <- .71 * training$num415 + .71 * training$num857
Y <- .71 * training$num415 - .71 * training$num857
plot(X, Y)


#SVD法：解釋強、干擾大、壓抑小
#PCA法：解釋弱、干擾小、壓抑大

smallSpam <- spam[, c(34, 32)]
prComp <- prcomp(smallSpam)
plot(prComp$x[, 1], prComp$x[, 2])

prComp$rotation #PC1: .7081 * num415 + .7061 * num857

typeColor <- ((spam$type == "spam") * 1 + 1)
prComp <- prcomp(log10(spam[, -58] + 1))
plot(prComp$x[, 1], prComp$x[, 2], col = typeColor, xlab = "PC1", ylab = "PC2")


preProc <- preProcess(log10(spam[, -58] + 1), method = "pca", pcaComp = 2) #用preProcess()做PCA
spamPC <- predict(preProc, log10(spam[, -58] + 1))
plot(spamPC[, 1], spamPC[, 2], col = typeColor)

preProc <- preProcess(log10(training[, -58] + 1), method = "pca", pcaComp = 2)
trainPC <- predict(preProc, log10(training[, -58] + 1))
modelFit <- train(training$type ~., method = "glm", data = trainPC)

testPC <- predict(preProc, log10(testing[, -58] + 1))
confusionMatrix(testing$type, predict(modelFit, testPC))

modelFit <- train(training$type ~., method = "glm", preProcess = "pca", data = training)
confusionMatrix(testing$type, predict(modelFit, testing))

#主成份分析可能會使結果更難解釋
#最好先轉換資料（logs/Box Cox）
#可以把圖畫出來確認問題在哪
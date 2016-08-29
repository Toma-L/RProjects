# R Case Studies =====


## 4 - 決策樹與隨機森林 =====

### 4.1 party套件 =====

str(iris)
set.seed(1234)
ind <- sample(2, nrow(iris), replace = TRUE, prob = c(0.7, 0.3)) # 建立index
trainData <- iris[ind == 1, ] # 分出訓練組
testData <- iris[ind == 2, ]

library(party)
# detach(package:mvtnorm)
# unloadNamespace("mvtnorm")

myFormula <- Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width
iris_ctree <- ctree(myFormula, data = trainData)
table(predict(iris_ctree), trainData$Species)

print(iris_ctree) # 輸出規則

plot(iris_ctree)
plot(iris_ctree, type = "simple")


testPred <- predict(iris_ctree, newdata = testData)
table(testPred, testData$Species)


### 4.2 rpart套件 =====

# install.packages("mboost")
library(mboost)
data("bodyfat", package="TH.data")
dim(bodyfat)
attributes(bodyfat)
bodyfat[1:5, ]

set.seed(1234)
ind <- sample(2, nrow(bodyfat), replace = TRUE, prob = c(.7, .3))
bodyfat.train <- bodyfat[ind == 1, ]
bodyfat.test <- bodyfat[ind == 2, ]
# install.packages("rpart")
library(rpart)
myFormula <- DEXfat ~ age + waistcirc + hipcirc + elbowbreadth + kneebreadth

bodyfat_rpart <- rpart(myFormula, data = bodyfat.train, 
                       control = rpart.control(minsplit = 10))
attributes(bodyfat_rpart)
print(bodyfat_rpart$cptable)

print(bodyfat_rpart)
plot(bodyfat_rpart)
text(bodyfat_rpart, use.n = TRUE)


opt <- which.min(bodyfat_rpart$cptable[, "xerror"]) # 選擇最小誤差的決策樹
cp <- bodyfat_rpart$cptable[opt, "CP"]
bodyfat_prune <- prune(bodyfat_rpart, cp = cp)
print(bodyfat_prune)

plot(bodyfat_prune)
text(bodyfat_prune, use.n = TRUE)


DEXfat_pred <- predict(bodyfat_prune, newdata = bodyfat.test)
xlim <- range(bodyfat$DEXfat)
plot(DEXfat_pred ~ DEXfat, data = bodyfat.test, xlab = "Observed", ylab = "Predicted", ylim = xlim, xlim = xlim)
abline(a = 0, b = 1)


### 4.3 randomForest套件 =====

ind <- sample(2, nrow(iris), replace = TRUE, prob = c(.7, .3))
trainData <- iris[ind == 1, ]
testData <- iris[ind == 2, ]

library(randomForest)
rf <- randomForest(Species ~ ., data = trainData, ntree = 100, proximity = TRUE)
table(predict(rf), trainData$Species)
print(rf)
attributes(rf)
plot(rf)
importance(rf) # 變數的重要性
varImpPlot(rf) # 變數的重要性（圖示）

irisPred <- predict(rf, newdata = testData)
table(irisPred, testData$Species)

plot(margin(rf, testData$Species)) # 查看預測結果


## 10 - 文本挖掘 =====

### Twitter文本檢索 =====

library(RCurl)
library(twitteR)

consumerKey <- "GlSaLonfRtsFE3KOtcjTfSwMZ"
consumerSecret <- "yaTIBuqWR3JMyiwk0E1UQzCUbJoYgHDVzexxaJhTH2Sszc3vUj"
accessToken <- "1177664006-ka3vw49zEQk7fgCLenP4yRDlsV0DTzlEFCjltT8"
accessTokenSecret <- "1fRWUNjjGNDzFqKxInRS7wSww2gTdn5v8wjyvUiFumHwZ"

setup_twitter_oauth(consumerKey, consumerSecret, accessToken, accessTokenSecret)

rdmTweets <- userTimeline("rdatamining", n = 200)
(nDocs <- length(rdmTweets))

rdmTweets[11:15]
rdmTweets[[11]]$getText()
# ?cat()

for(i in 11:15) {
        cat(paste("[[", i, "]]", sep = ""))
        writeLines(strwrap(rdmTweets[[i]]$getText(), width = 73))
}

### 轉換文本 =====

df <- do.call("rbind", lapply(rdmTweets, as.data.frame))
dim(df)

library(tm)
myCorpus <- Corpus(VectorSource(df$text))# 轉換為語料庫

myCorpus <- tm_map(myCorpus,
                   content_transformer(function(x) iconv(x, to = 'UTF-8-MAC', sub = 'byte')),
                   mc.cores = 1)

inspect(myCorpus)

myCorpus <- tm_map(myCorpus, tolower, lazy = TRUE) # 轉為小寫
myCorpus <- tm_map(myCorpus, removePunctuation) # 移除標點
myCorpus <- tm_map(myCorpus, removeNumbers)
removeURL <- function(x) gsub("http[[:alnum:]]*", "", x)
myCorpus <- tm_map(myCorpus, removeURL) # 移除超連結
myStopwords <- c(stopwords('english'), "available", "via") # 設定停用字
myStopwords <- setdiff(myStopwords, c("r", "big"))
myCorpus <- tm_map(myCorpus, removeWords, myStopwords)

getTransformations() # 查看所有可用的轉換


### 提取詞幹 =====

myCorpusCopy <- myCorpus
myCorpus <- tm_map(myCorpus, stemDocument, lazy = TRUE) # 提取詞幹
inspect(myCorpus[11:15])

for(i in 11:15) {
        cat(paste("[[", i, "]] ", sep = ""))
        writeLines(strwrap(myCorpus[[i]], width = 73))
}

### 有問題
myCorpus <- tm_map(myCorpus, stemCompletion, dictionary = myCorpusCopy, lazy = TRUE)
inspect(myCorpus[11:15])

miningCases <- tm_map(myCorpusCopy, grep, pattern = "\\<mining")
sum(unlist(miningCases))
minerCases <- tm_map(myCorpusCopy, grep, pattern = "\\<miners")
sum(unlist(minerCases))

myCorpus <- tm_map(myCorpus, content_transformer(gsub), pattern = "miners", replacement = "mining") # 用mining取代miner


### 建立詞項 - 文檔矩陣 =====

# myCorpus <- tm_map(myCorpus, PlainTextDocument)
# myCorpus <- tm_map(myCorpus, content_transformer)

myTdm <- TermDocumentMatrix(myCorpus, control = list(wordLengths = c(1, Inf)))
myTdm # 此矩陣可用於分群、分類及關聯規則分析

idx <- which(dimnames(myTdm)$Terms == "r")
inspect(myTdm[idx+(0:5), 101:110])

rownames(myTdm) # 詞項列表


### 頻繁詞項與關聯 =====


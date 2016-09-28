#Getting and Cleaning Data==================================================

#Using data.table

library(data.table)
DF <- data.frame(x = rnorm(9), y = rep(c("a", "b", "c"), each = 3), z = rnorm(9))
head(DF, 3)
DT <- data.table(x = rnorm(9), y = rep(c("a", "b", "c"), each = 3), z = rnorm(9))
head(DT, 3)
tables() #查看現在有哪些data.table
DT[2, ]
DT[DT$y == "a", ]
DT[c(2, 3)] #都是row
DT[, c(2, 3)] #執行c(2, 3)
{
        x = 1
        y = 2
}
k <- {print(10); 5} #10，會執行前面的
print(k) #5
DT[, list(mean(x), sum(z))]
DT[, table(y)]
DT[, w := z^2] #使用:=指定新的變數存在原本的資料裡
DT2 <- DT
DT[, y := 2] #y全變2
DT
head(DT, 3)
head(DT2, 3) #改變DT也會改變到DT2
DT[, m := {tmp <- (x+z); log2(tmp + 5)}]
head(DT)

DT[, a := x > 0]
head(DT)

DT[, b := mean(x + w), by = a] #把資料分成a為TRUE和a為FALSE，各自計算mean(x+w)，當成其b的值

set.seed(123)
DT <- data.table(x = sample(letters[1:3], 1E5, TRUE))
DT[, .N, by = x] #用.N計算次數

DT <- data.table(x = rep(c("a", "b", "c"), each = 100), y = rnorm(300))
setkey(DT, x)
DT['a']

DT1 <- data.table(x = c('a', 'a', 'b', 'dt1'), y = 1:4)
DT2 <- data.table(x = c('a', 'b', 'dt2'), z = 5:7)
setkey(DT1, x)
setkey(DT2, x)
merge(DT1, DT2) #只留下都有的

big_df <- data.frame(x = rnorm(1E6), y = rnorm(1E6))
file <- tempfile()
write.table(big_df, file = file, row.names = FALSE, col.names = TRUE, sep"\t", quote = FALSE)
system.time(fread(file)) #fread()類似read.table()
system.time(read.table(file, header = TRUE, sep = "\t"))

#Subsetting and sorting

set.seed(13435)
X <- data.frame("var1" = sample(1:5), "var2" = sample(6:10), "var3" = sample(11:15))
X <- X[sample(1:5), ]
X$var2[c(1, 3)] <- NA
X

X[, 1]
X[, "var1"]
X[1:2, "var2"]
X[(X$var1 <= 3 & X$var3 > 11), ]
X[(X$var1 <= 3 | X$var3 > 15), ]
X[which(X$var2 > 8), ]
sort(X$var1)
sort(X$var1, decreasing = TRUE)
sort(X$var2, na.last = TRUE)
X[order(X$var1), ]
X[order(X$var1, X$var3), ]

library(plyr)
arrange(X, var1)
arrange(X, desc(var1))

X$var4 <- rnorm(5)
X

Y <- cbind(X, rnorm(5))
Y

#Summarizing data

if(!file.exists("./data")) {dir.create("./data")}
fileUrl <- "http://data.baltimorecity.gov/api/views/k5ry-ef3g/rows.csv?accessType = DOWNLOAD"
download.file(fileUrl, destfile = "./data/restaurants.csv")
restData <- read.csv("./data/restaurants.csv")
head(restData, n = 3)
tail(restData, n = 3)
summary(restData)
str(restData)

quantile(restData$councilDistrict, na.rm = TRUE)
table(restData$zipCode, useNA = "ifany") #如果有NA的話
table(restData$councilDistrict, restData$zipCode)

sum(is.na(restData$councilDistrict))
any(is.na(restData$councilDistrict))
all(restData$zipCode > 0)

colSums(is.na(restData))
all(colSums(is.na(restData)) == 0)

table(restData$zipCode %in% c("21212")) #%in%前面有多少符合後面
table(restData$zipCode %in% c("21212", "21213"))

restData[restData$zipCode %in% c("21212", "21213"), ]

data(UCBAdmissions)
DF <- as.data.frame(UCBAdmissions)
summary(DF)

xt <- xtabs(Freq ~ Gender + Admit, data = DF) #cross tabs

warpbreaks$replicate <- rep(1:9, len = 54)
xt <- xtabs(breaks ~., data = warpbreaks)
xt

ftable(xt)

fakeData <- rnorm(1e5)
object.size(fakeData) #資料大小

#Creating new variables

s1 <- seq(1, 10, by = 2); s1
s2 <- seq(1, 10, length = 3); s2
x <- c(1, 3, 8, 25, 100); seq(along = x)
restData$nearMe <- restData$neighborhood %in% c("Roland Park", "Homeland")
table(restData$nearMe)

restData$zipWrong <- ifelse(restData$zipCode < 0, TRUE, FALSE)
table(restData$zipWrong, restData$zipCode < 0)

restData$zipGroups <- cut(restData$zipCode, breaks = quantile(restData$zipCode)) #cut() # 可加入 label = () 參數
table(restData$zipGroups)
table(restData$zipGroups, restData$zipCode)

library(Hmisc)
restData$zipGroups <- cut2(restData$zipCode, g = 4) #easier cutting
table(restData$zipGroups)

restData$zcf <- factor(restData$zipCode)
restData$zcf[1:10]
class(restData$zcf)

yesno <- sample(c("yes", "no"), size = 10, replace = TRUE)
yesnofac <- factor(yesno, levels = c("yes", "no")) #levels為factor排序
relevel(yesnofac, ref = "yes")

as.numeric(yesnofac)

library(Hmisc)
restData$zipGroups <- cut2(restData$zipCode, g = 4)
table(restData$zipGroups)

library(Hmisc)
library(plyr)
restData2 <- mutate(restData, zipGroups = cut2(zipCode, g = 4)) #mutate()創造新變數
table(restData2$zipGroups)

#Reshaping data

library(reshape2)
head(mtcars)
mtcars$carname <- rownames(mtcars)
carMelt <- melt(mtcars, id = c("carname", "gear", "cyl"), measure.vars = c("mpg", "hp"))
head(carMelt, n = 3)
tail(carMelt, n = 3)

cylData <- dcast(carMelt, cyl ~ variable)
cylData

cylData <- dcast(carMelt, cyl ~ variable, mean) #對variable的值做mean
cylData

head(InsectSprays)
tapply(InsectSprays$count, InsectSprays$spray, sum)

spIns <- split(InsectSprays$count, InsectSprays$spray)
spIns
sprCount <- lapply(spIns, sum)
sprCount
unlist(sprCount)

sapply(spIns, sum)

ddply(InsectSprays, .(spray), summarize, sum = sum(count))

spraySums <- ddply(InsectSprays, .(spray), summarize, sum = ave(count, FUN = sum))
head(spraySums)

#Merging data

if(!file.exists("./data")) {dir.create("./data")}
fileUrl1 <- "http://dl.dropboxusercontent.com/u/7710864/data/reviews-apr29.csv"
fileUrl2 <- "http://dl.dropboxusercontent.com/u/7710864/data/solutions-apr29.csv"
download.file(fileUrl1, destfile = "./data/reviews.csv")
download.file(fileUrl2, destfile = "./data/solutions.csv")
reviews <- read.csv("./data/reviews.csv")
solutions <- read.csv("./data/solutions.csv")
head(reviews, 2)
head(solutions, 2)
names(reviews)
names(solutions)

mergedData <- merge(reviews, solutions, by.x = "solution_id", by.y = "id", all = TRUE)
head(mergedData)
intersect(names(solutions), names(reviews)) #intersect()找共同都有的變數

mergedData2 <- merge(reviews, solutions, all = TRUE) #不設定會自動merge名字一樣的變數
head(mergedData2)

df1 <- data.frame(id = sample(1:10), x = rnorm(10))
df2 <- data.frame(id = sample(1:10), y = rnorm(10))
arrange(join(df1, df2), id) #用join()合併

df1 <- data.frame(id = sample(1:10), x = rnorm(10))
df2 <- data.frame(id = sample(1:10), y = rnorm(10))
df3 <- data.frame(id = sample(1:10), z = rnorm(10))
dfList <- list(df1, df2, df3)
join_all(dfList) #很多個data.frame用join_all合併

#Editing text variables

cameraData <- read.csv("./data/cameras.csv")
names(cameraData)

splitNames <- strsplit(names(cameraData), "\\.")
splitNames[[5]]
splitNames[[6]]

mylist <- list(letters = c("A", "b", "c"), numbers = 1:3, matrix(1:25, ncol = 5))
head(mylist)
mylist[1]
mylist$letters
mylist[[1]]

splitNames[[6]][1]
firstElement <- function(x) {x[1]}
sapply(splitNames, firstElement)

names(reviews)
sub("_", "", names(reviews), ) #用""取代掉"_"，用sub()

testName <- "this_is_a_test"
sub("_", "", testName) #sub只會作用在第1個

gsub("_", "", testName) #gsub()作用在全部

grep("Alameda", cameraData$intersection) #grep()
table(grepl("Alameda", cameraData$intersection)) #grepl()變布林值

grep("Alameda", cameraData$intersection, value = TRUE) #直接給值，不要編號
grep("JeffStreet", cameraData$intersection)
length(grep("JeffStreet", cameraData$intersection))

library(stringr)
nchar("Jeffrey Leek") #計算字串長度
substr("Jeffrey Leek", 1, 7) #只取字串的某部分
paste("Jeffrey", "Leek")
paste0("Jeffrey", "Leek") #paste0()不會有空白
str_trim("Jeff       ") #清除空白

#Regular Expressions

## ^i think：以i think為首的
## morning$：以morning作結的
## [Bb][Uu][Ss][Hh]：不論大小寫
## ^[Ii] am
## ^[0-9][a-zA-Z]
## [^?.]$：[]放開頭代表排除

#Regular Expressions II

## 9.11：.代表什麼字都可以
## flood|fire：|表示「或」
## flood|earthquake|hurricane|coldfire
## ^[Gg]ood|[Bb]ad：bad或Bad就未必在句首
## ^([Gg]ood|[Bb]ad)：good和bad都只要句首
## [Gg]eorge([Ww]\.)? [Bb]ush：?代表「未必要」這個
## [Gg]eorge([Ww]\.)? [Bb]ush：加\代表.是純文字
## (.*)：.代表()裡面可以是任何字，*代表「不限字數」
## [0-9]+(.*)[0-9]+：+代表「至少要有一個」
## [Bb]ush( +[^ ]+ +){1,5}：{1,5}代表重複1~5次

## +([a-zA-Z]+) +\1 +
## ^s(.*)s
## ^s(.*?)s$

#Working with dates

d1 <- date()
d1
class(d1) #character

d2 <- Sys.Date()
d2
class(d2) #Date

format(d2, "%a %b %d")

x <- c("1jan1960", "2jan1960", "31mar1960", "30jul1960")
z <- as.Date(x, "%d%b%Y")
z
z[1] - z[2]
as.numeric(z[1] - z[2])

weekdays(d2)
months(d2)
julian(d2)

library(lubridate)
ymd("20140108")
mdy("08/04/2013")
dmy("03-04-2013")
ymd_hms("2011-08-03 10:15:03")
ymd_hms("2011-08-03 10:15:03", tz = "Pacific/Auckland")
?Sys.timezone

x <- dmy(c("1jan2013", "2jan2013", "31mar2013", "30jul2013"))
wday(x[1])
wday(x[1], label = TRUE)


#利用R語言打通大數據的經脈==========

#資料前置處理==========

##資料集載入==========

install.packages("lattice")
install.packages("MASS")
install.packages("nnet")
library(lattice)
library(MASS)
library(nnet)

install.packages("mice")
library(mice)
data(nhanes2)
nrow(nhanes2)
ncol(nhanes2)
summary(nhanes2)
head(nhanes2)


##資料清理==========

pay <- c(11, 19, 14, 22, 14, 28, 13, 81, 12, 43, 11, 16, 31, 16, 23, 42, 22, 26, 17, 22, 13, 27, 180, 
         16, 43, 82, 14, 11, 51, 76, 28, 66, 29, 14, 14, 65, 37, 16, 37, 35, 39, 27, 14, 17, 13, 38, 
         28, 40, 85, 32, 25, 26, 16, 12, 54, 40, 18, 27, 16, 14, 33, 29, 77, 50, 19, 34)
par(mfrow = c(2, 2))
hist(pay)
dotchart(pay)
boxplot(pay, horizontal = TRUE)
qqnorm(pay); qqline(pay)


###遺漏值處理==========

sum(is.na(nhanes2)) #計算有多少樣本有遺漏值
sum(complete.cases(nhanes2)) #計算完整樣本的數量
md.pattern(nhanes2) #觀察遺漏值的情況

#若要以chl為依變數，以age、hyp、bmi為自變數建立線性迴歸模型，因nhanes2有遺漏值，
#不能直接建構，因此要用mice()函數，用下列方式建模型

imp <- mice(nhanes2, m = 4) #產生4組完整資料
#多重插補法（Multivariate Imputation）透過變數間的關係對缺失資料進行預測，
#利用蒙地卡羅方法產生多個完整資料集，再對這些資料集分別進行分析。

fit <- with(imp, lm(chl ~ age + hyp + bmi))
pooled <- pool(fit) #整理4個迴歸模型
summary(pooled)

####刪除法==========
####插補法==========

sub <- which(is.na(nhanes2[, 4] == TRUE))
dataTR <- nhanes2[-sub, ]
dataTE <- nhanes2[sub, ]
dataTE[, 4] <- sample(dataTR[, 4], length(dataTE[, 4]), replace = TRUE) #隨機插補法
dataTE

sub <- which(is.na(nhanes2[, 4]) == TRUE)
dataTR <- nhanes2[-sub, ]
dataTE <- nhanes2[sub, ]
dataTE[, 4] <- mean(dataTR[, 4]) #平均值法
dataTE

sub <- which(is.na(nhanes2[, 4]) == TRUE)
dataTR <- nhanes2[-sub, ]
dataTE <- nhanes2[sub, ]
dataTE
lm <- lm(chl ~ age, data = dataTR)
nhanes2[sub, 4] <- round(predict(lm, dataTE)) #迴歸模型插補法
head(nhanes2)

#熱平台插補法是在非遺漏值資料集中找一個與遺漏值所在樣本相似的樣本
accept <- nhanes2[which(apply(is.na(nhanes2), 1, sum) != 0), ] #有遺漏值的樣本
donate <- nhanes2[which(apply(is.na(nhanes2), 1, sum) == 0), ] #無遺漏值樣本
accept[1, ]
donate[1, ]
accept[2, ]
sa <- donate[which(donate[, 1] == accept[2, 1] & donate[, 3] == accept[2, 3] & donate[, 4] == accept[2, 4]), ]
#尋找跟accept第2個樣本相似的樣本
accept[2, 2] = sa[1, 2]
accept[2, ]
#熱平台插補法在變數數量很多時很難找到需要的樣本

#冷平台插補法是按照某些變數將資料分層，在層中對遺漏值使用平均值插補法
level1 <- nhanes2[which(nhanes2[, 3] == "yes"), ] #按照變數hyp分層，挑出hyp值為yes的樣本
level1

level1[4, 4] <- mean(level1[1:3, 4])
level1


###雜訊資料處理==========

install.packages("outliers")
library(outliers)
set.seed(1)
s1 <- .Random.seed #產生100個標準常態亂數
y <- rnorm(100)
outlier(y)
outlier(y, opposite = TRUE)
dotchart(y)
dim(y) <- c(20, 5) #將y的資料重新劃分成20*5的矩陣
outlier(y) #求每列的離群最遠值
outlier(y, opposite = TRUE) #離群最遠值的相反值

set.seed(1)
s1 <- .Random.seed
y <- rnorm(10)
outlier(y, logical = TRUE)
plot(y)

#離群值也可以用分群方法檢測

#雜訊檢查後，常用分箱、迴歸、電腦檢查和人工方法檢查結合等方法「光滑」資料

#分箱方法透過對資料排序，利用資料「近鄰」來光滑有序數據的局部光滑方法

set.seed(1)
s1 <- .Random.seed
x <- rnorm(12)
x <- sort(x)
dim(x) <- c(3, 4)
x[1, ] <- apply(x, 1, mean)[1]
x[2, ] <- apply(x, 1, mean)[2]
x[3, ] <- apply(x, 1, mean)[3]
x


###資料不一致的處理==========

x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE, FALSE, FALSE, TRUE))
x

probs <- c(1:3/4)
rt.value <- c(0, 0, 0) #設定傳回值為3個數字
vapply(x, quantile, FUN.VALUE = rt.value, probs = probs)

probs <- c(1:4/4)
vapply(x, quantile, FUN.VALUE = rt.value, probs = probs) #ERROR，傳回值與要求格式不一

rt.value <- c(0, 0, 0, 0) #設定傳回值為4個數字
vapply(x, quantile, FUN.VALUE = rt.value, probs = probs)
rt.value <- c(0, 0, 0, "") #設定傳回值為3個數字、1個字串


##資料整合==========

x <- cbind(sample(c(1:50), 10), sample(c(1:50), 10))
chisq.test(x) #對矩陣x進行卡方檢定，檢查兩列是否相關（Reject H0，變數不相關）

x <- cbind(rnorm(10), rnorm(10))
cor(x) #相關係數，相關係數大的話，其中一個可被視為容錯而刪除
cov(x) #協方差，正值表示兩個屬性趨向一起改變，負值則相反方向

x <- cbind(sample(c(1:10), 10, replace = TRUE), rnorm(10), rnorm(10))
head(x)
y <- unique(x[, 1])
sub <- rep(0, length(y))
for(i in 1:length(y)) {
        sub[i] <- which(x[, 1] == y[i])[1]
} #去掉重複觀測值
x <- x[sub, ]
head(x)


##資料轉換==========

#光滑：分箱、迴歸、分群等
#屬性建構：由指定的屬性建構出新屬性
#聚集：整理資料
#規範化：把資料按比例縮放，實質落入一個特定的小區間

set.seed(1)
s1 <- .Random.seed #設定亂數種子
a <- rnorm(5) #產生亂數
b <- scale(a) #對亂數標準化
b
a


#離散化：數值屬性的原始值用區間標籤或概念化標籤取代，將定量資料轉成定性資料

a <- c(.7063422, .7533599, .6675749, .6100253, .9341495, .6069284, .3462011)
n <- length(a)
la <- rep(0, n)
la[which(a > .5)] = 1
la #0~1的連續性資料，轉為0或1的離散型資料

#由額定資料產生概念分層：屬性，資料泛化可以視為資料合併

city <- c(6, 7, 6, 2, 2, 6, 2, 1, 5, 7, 2, 1, 1, 6, 1, 3, 8, 8, 1, 1)
province <- rep(0, 20)
province[which(city > 4)] <- 1
province


##資料精簡

#AIC準則（赤池資訊準則）：透過最佳模型選擇變數（AIC越小越好）
#LASSO：透過限制條件選擇變數
#分類樹、隨機森林：透過對分類效果的影響大小篩選屬性
#小波轉換、主成份分析：把原資料轉換或投影到較小的空間

x <- matrix(rnorm(100 * 20), 100, 20) #產生20行常態亂數
y <- rnorm(100) #產生常態亂數作為依變數
library(glmnet)
fit1 <- glmnet(x, y) #廣義線性回歸，自變數為分組的，預設為LASSO
b <- coef(fit1, s = .01) #s代表lambda值，lambda越小，約束越寬，篩選出的變數越多
b #b為變數係數，有值的代表有選入模型

predict(fit1, newx = x[1:10, ], s = c(.01, .005)) #lambda為0.01和0.005情況下的預測值

#數值精簡指用較小的資料表示形式取代原資料。如參數方法中使用模型估算資料，就可以只儲存模型參數代替實際資料，如回歸模型和對數線性模型
#非參數方法可以使用長條圖、分群、抽樣和資料立方體聚集等方法


#R軟體資料分析基礎與應用==========

#cbind和rbind資料合併==========

sport <- c("Hockey", "Baseball", "Football")
league <- c("NHL", "MLB", "NFL")
trophy <- c("Stanley Cup", "Commissioner's Trophy", "Vince Lombardi Trophy")
trophies1 <- cbind(sport, league, trophy)
trophies2 <- data.frame(sport = c("Basketball", "Golf"), league = c("NBA", "PGA"), trophy = c("Larry O'Brien Championship Trophy", "Wanamaker Trophy"), stringsAsFactors = FALSE)
trophies <- rbind(trophies1, trophies2)

cbind(Sport = sport, Association = league, Prize = trophy)


#資料連結==========

download.file(url = "http://jaredlander.com/data/US_Foreign_Aid.zip", destfile = "ForeignAid.zio")
unzip("ForeignAid.zip", exdir = "data12")

require(stringr)
theFiles <- dir("data12/", pattern = "\\.csv") #取得檔案的列表
for(a in theFiles) {
        nameToUse <- str_sub(string = a, start = 12, end = 18)
        temp <- read.table(file = file.path("data12", a), header = TRUE, #用file.path()來指定資料夾和檔名
                           sep = ",", stringsAsFactors = FALSE)
        assign(x = nameToUse, value = temp)
}


##用merge合併兩個data.frame==========

Aid90s00s <- merge(x = Aid_90s, y = Aid_00s, by.x = c("Country.Name", "Program.Name"),
                   by.y = c("Country.Name", "Program.Name")) #merge()可指定不同名稱的變數進行連結，但速度比其他方法慢
head(Aid90s00s)


##用plyr join合併data.frame==========

require(plyr)
Aid90s00sJoin <- join(x = Aid_90s, y = Aid_00s, 
                      by = c("Country.Name", "Program.Name")) #join()速度比較快，但變數名稱要一樣
head(Aid90s00sJoin)

frameNames <- str_sub(string = theFiles, start = 12, end = 18)
frameList <- vector("list", length(frameNames)) #建立空的list
names(frameList) <- frameNames
for(a in frameNames) {
        frameList[[a]] <- eval(parse(text = a)) #parse解析字元，轉換為變數
}

head(frameList[[1]]) #data.frame都在list裡了

head(frameList[["Aid_00s"]])
head(frameList[[5]])
head(frameList[["Aid_60s"]])

allAid <- Reduce(function(...){ #用Reduce()函數把所有元素連結在一起，加快運行速度
        join(..., by = c("Country.Name", "Program.Name"))
}, frameList) #Reduce的原理是把1跟2先加起來，再把3加進來
dim(allAid)
require(useful)
corner(allAid, c = 15)
bottomleft(allAid, c = 15)


##data.table中的資料合併==========

require(data.table)
dt90 <- data.table(Aid_90s, key = c("Country.Name", "Program.Name"))
dt00 <- data.table(Aid_00s, key = c("Country.Name", "Program.Name"))
dt0090 <- dt90[dt00] #90在左，00在右，該指令是左連結，連結時需要用到關鍵詞，在建立data.table時已指定


#用reshape2套件置換行、列資料==========

##melt==========

head(Aid_00s)
require(reshape2)
melt00 <- melt(Aid_00s, id.vars = c("Country.Name", "Program.Name"), variable.name = "Year", value.name = "Dollars")
tail(melt00, 10)

require(scales)
melt00$Year <- as.numeric(str_sub(melt00$Year, start = 3, 6)) #把Year中的FY去掉，轉成numeric
meltAgg <- aggregate(Dollars ~ Program.Name + Year, data = melt00, sum, na.rm = TRUE)
meltAgg$Program.Name <- str_sub(meltAgg$Program.Name, start = 1, end = 10)
ggplot(meltAgg, aes(x = Year, y = Dollars)) + 
        geom_line(aes(group = Program.Name)) +
        facet_wrap(~ Program.Name) +
        scale_x_continuous(breaks = seq(from = 2000, to = 2009, by = 2)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0)) +
        scale_y_continuous(labels = multiple_format(extra = dollar, multiple = "B"))


##dcast==========

cast00 <- dcast(melt00, Country.Name + Program.Name ~ Year, #~左邊放要保留的，右邊放新指定的
                value.var = "Dollars")
head(cast00)


#字串處理==========

##用paste建立字串==========

paste("Hello", "Jared", "and others")
paste("Hello", "Jared", "and others", sep = "/")
paste(c("Hello", "Hey", "Howdy"), c("Jared", "Bob", "David")) #paste()也允許向量化運算
paste("Hello", c("Jared", "Bob", "David"))
paste("Hello", c("Jared", "Bob", "David"), c("Goodbye", "Seeya"))
vectorOfText <- c("Hello", "Everyone", "out there", ".")
paste(vectorOfText, collapse = " ") #collapse參數將文字折疊


##用sprintf建立含有變數的字串==========

person <- "Jared"
partySize <- "eight"
waitTime <- 25

paste("Hello ", person, ", your party of ", partySize, 
      " will be seated in ", waitTime, " minutes.", sep = " ") #很麻煩

sprintf("Hello %s, your party of %s will be seated in %s minutes", person, partySize, waitTime)
sprintf("Hello %s, your party of %s will be seated in %s minutes", 
        c("Jared", "Bob"), c("eight", 16, "four", 10), waitTime)


##抽取文字==========

require(XML)
theURL <- "http://www.loc.gov/rr/print/list/057_chron.html"
presidents <- readHTMLTable(theURL, which = 3, as.data.frame = TRUE, 
                            skip.rows = 1, header = TRUE, stringsAsFactors = FALSE)

head(presidents)
tail(presidents$YEAR)
presidents <- presidents[1:64, ]
require(stringr)
yearList <- str_split(string = presidents$YEAR, pattern = "-") #拆開字串
head(yearList)

yearMatrix <- data.frame(Reduce(rbind, yearList)) #合併成一個matrix
head(yearMatrix)
names(yearMatrix) <- c("Start", "Stop")
presidents <- cbind(presidents, yearMatrix)
presidents$Start <- as.numeric(as.character(presidents$Start))
presidents$Stop <- as.numeric(as.character(presidents$Stop))
head(presidents)
tail(presidents)

str_sub(string = presidents$PRESIDENT, start = 1, end = 3)
str_sub(string = presidents$PRESIDENT, start = 4, end = 8)
presidents[str_sub(string = presidents$Start, start = 4, end = 4) == 1, c("YEAR", "PRESIDENT", "Start", "Stop")]


##正規表示法==========

johnPos <- str_detect(string = presidents$PRESIDENT, pattern = "John") #不知道字會出現在什麼地方用str_detect()
presidents[johnPos, c("YEAR", "PRESIDENT", "Start", "Stop")]

badSearch <- str_detect(presidents$PRESIDENT, "john")
goodSearch <- str_detect(presidents$PRESIDENT, ignore.case("John")) #若要忽略大小寫，ignore.case
sum(badSearch)
sum(goodSearch)

con <- url("http://jaredlander.com/data/warTimes.rdata")
load(con) #用load()載入
close(con) #用close()關閉連結

head(warTimes, 10)

warTimes[str_detect(string = warTimes, pattern = "-")]
theTimes <- str_split(string = warTimes, pattern = "(ACAEA)|-", n = 2) #維基百科用的分隔符號為ACAEA
head(theTimes)

which(str_detect(string = warTimes, pattern = "-"))
theTimes[[147]]
theTimes[[150]]

theStart <- sapply(theTimes, FUN = function(x) x[1])
head(theStart)

theStart <- str_trim(theStart) #移除不必要的空格
head(theStart)

str_extract(string = theStart, pattern = "January") #抽出含有January字串的地方，用str_extract()
theStart[str_detect(string = theStart, pattern = "January")] #找出含有January的元素並整項回傳，用str_detect()

head(str_extract(string = theStart, "[0-9][0-9][0-9][0-9]"), 20) #找四個數字一起出現的情況
head(str_extract(string = theStart, "[0-9]{4}"), 20) #搜尋由4位數組成的數字
head(str_extract(string = theStart, "\\d{4}"), 20) #\\d代表任何整數
str_extract(string = theStart, "\\d{1,3}") #搜尋任何一個出現過1次、2次或3次的數字
head(str_extract(string = theStart, pattern = "^\\d{4}"), 30) #前端用^後端用$
head(str_extract(string = theStart, pattern = "\\d{4}$"), 30)
head(str_extract(string = theStart, pattern = "^\\d{4}$"), 30)
head(str_replace(string = theStart, pattern = "\\d", replacement = "x"), 30) #將第一個數字取代為x
head(str_replace_all(string = theStart, pattern = "\\d", replacement = "x"), 30) #將所有數字取代為x
head(str_replace_all(string = theStart, pattern = "nndf1, 4g", replacement = "x"), 30) #取代任何由1位數到4位數組成的數字為x

commands <- c("<a href = index.html>The Link is here</a>", "<b>This is bold text</b>")
str_replace(string = commands, pattern = "<.+?>(.+?)<.+>", replacement = "\\1")
#.代表搜尋任何文字
#+代表至少搜尋1次
#?代表非貪婪搜尋，採用貪婪搜尋若字串中多組符合匹配結果，可能會抓到額外內容
#\\1是回溯引用的符號，將搜尋的規則取代掉，1代表用搜尋的第一組文字做取代
#完成整詳見?regex


# DataCamp =====

# 資料框整理技巧 =====

library(tidyr)
straw_hat_long_df <- gather(straw_hat_wide_df, key = var, value = name, age, height)
straw_hat_wide_df <- spread(straw_hat_long_df, key = cate, value = int)

straw_hat_df_devil_fruit <- merge(straw_hat_df, straw_hat_devil_fruit, all.x = TRUE) # 保留左邊所有row


# Cleaning Data in R =====

library(dplyr)
glimpse(bmi)

bmi_long <- gather(bmi, year, bmi_val, -Country)
bmi_wide <- spread(bmi_long, year, bmi_val)

bmi_cc_clean <- separate(bmi_cc, col = Country_ISO, into = c("Country", "ISO"), sep = "/")
bmi_cc <- unite(bmi_cc_clean, Country_ISO, Country, ISO, sep = "-")
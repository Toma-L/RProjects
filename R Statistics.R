#R軟體統計應用分析實務==============================

#資料的處理與轉換==============================

##匯入與儲存資料檔==============================

dsc <- read.csv("sc_1.csv", header = T)
head(dsc, 5)
colnames(dsc)
names(dsc)
lapply(dsc, length)
sapply(dsc, length)
attach(dsc)
dsc[[6]]
dsc[6]
dsc[1:10, "esco"]
dsc[c(1, 25, 49, 73), "esco"]

attach(dsc)
search()
dsc$scav <- round((msco + esco) / 2, 1)
head(dsc)


##資料檔的分割==============================

cut(msco, breaks = 5, labels = NULL)
cut(msco, breaks = 5, labels = c(1:5))
cut(msco, breaks = 5, labels = c("E", "D", "C", "B", "A"))

table(cut(msco, breaks = 5, labels = c(1:5)))
table(cut(msco, breaks = 5, labels = c("E", "D", "C", "B", "A")))


##分割檔案函數==============================

split(esco, cla)
split(dsc, cla)
gro.cla <- split(dsc, cla)
head(gro.cla[[1]])
head(gro.cla[[2]])
head(gro.cla[[3]])

gro.disc <- split(dsc, disc)
head(gro.disc[[1]], 3)
head(gro.disc[[2]], 3)
head(gro.disc[[3]], 3)

mean(gro.disc[[1]]$esco)
mean(gro.disc[[2]]$esco)
mean(gro.disc[[3]]$esco)

gro.disc <- unstack(data.frame(esco, disc)) #與split()有類似功能
print(gro.disc[[1]]) #unstack(data.frame(計量變數, 因子變數))
print(gro.disc[[2]])
print(gro.disc[[3]])

mean(gro.disc[[1]])
mean(gro.disc[[2]])
mean(gro.disc[[3]])

tapply(esco, factor(disc), mean)


##資料檔重新編碼==============================

library(car)
dsc$sex <- recode(sex, "1 = 'Male'; 2 = 'Female'")
dsc$disc <- recode(disc, "1 = '民主'; 2 = '權威'; 3 = '放任'")
head(dsc, 3)
msco_g = recode(msco, "0:59 = 1; 60:69 = 2; 70:79 = 3; 80:89 = 4; 90:100 = 5")
table(msco_g)

dsc$g.msc <- recode(msco, "0:59 = 1; 60:69 = 2; 70:79 = 3; 80:89 = 4; 90:100 = 5")
head(dsc)

table(msco_g = recode(msco, "lo:59 = 1; 60:69 = 2; 70:79 = 3; 80:89 = 4; 90:hi = 5")) #最低lo，最高hi
table(msc_g <- recode(msco, "lo:59 = '不及格'; 60:hi = '及格'"))

table(msc_g <- recode(msco, "lo:59 = '不及格'; else = '及格'")) #可用else
table(pedu)
table
table(npedu = recode(pedu, "1 = 1; 2 = 2; 3 = 3; 4:5 = 4"))

dsc$n.pedu <- recode(pedu, "1 = '國中小'; 2 = '高中職'; 3 = '專科'; 4:5 = '大學以上'")
head(dsc)
write.csv(dsc, "dsc_reco.csv", row.names = F)
list.files()


##資料檔的排序與等級化==============================

a <- c(2, 13, 15, 10, NA, NA, 9, 4)
sort(a, decreasing = T)
sort(a, decreasing = F)
sort(a, decreasing = F, na.last = F)
sort(a, decreasing = F, na.last = T)
sort(a, decreasing = T, na.last = T)
sort(a, decreasing = T, na.last = F)

sort(esco, decreasing = T)
esco_s <- sort(esco, decreasing = T)
esco_s[30]
esco_s <- sort(esco, decreasing = F)
esco_s[30]
esco_s[27:30]

dsc$g.esco <- recode(esco, "lo:54 = 3; 55:74 = 2; 75:hi = 1")
print(dsc)

a <- c(4, 13, 15, 10, NA, NA, 9, 4, 13, 13, NA)
sa <- sort(a, decreasing = T, na.last = T)
rank(sa, ties.method = "average", na.last = NA)
rank(sa, ties.method = "min", na.last = NA)
rank(sa, ties.method = "max", na.last = T)

sa <- sort(a, decreasing = T, na.last = NA)
rank((sa <- sort(a, decreasing = T, na.last = NA)), ties.method = "random")

rank(msco)
dsc$rmsco <- (length(msco) + 1 - rank(msco, ties.method = "average"))

score <- c(45, 23, 65, 98, 72, 84, 69, 54, 61)
r.score <- rank(score)
print(r.score)
rl.score <- length(score) + 1 - r.score
print(rl.score)


##資料框架之資料型態的轉換==============================

dms <- read.csv("ms1.csv", header = T)
print(dms)

library(reshape)
dms_1 <- melt(dms, id.vars = c("id", "sex"), measure.vars = c("s1", "s2", "s3"))
dms_1
names(dms_1)[1:4] = c("編號", "性別", "次數", "測量值")
dms_1

dms_1$"次數" <- as.numeric(sub("s", "", dms_1$"次數"))
dms_1

st <- c("stu1", "stu2", "stu3", "stu4", "stu5") #sbu(原字串, 新字串, 向量名稱或變數)為取代
stu_1 <- sub("stu", "900", st)
stu_1

stu_2 <- as.numeric(stu_1)
class(stu_1)
class(stu_2)

wel <- c("data", "analysis")
nchar(wel) #計算文字串的字元數

bay_a <- paste("Bayesian", "data", "analysis")
print(bay_a)
bay_a <- paste("Bayesian", "data", "analysis", "002", sep = "-")
print(bay_a)

dabr <- c("daa", "analysis")
substr(dabr, 1, 3) #substr()抽取子字串


#次數分配==============================

dsc <- read.csv("dsc_reco.csv", header = T, fileEncoding = "BIG-5") #csv內含中文時，fileEncoding = "BIG-5"！！！
head(dsc)

attach(dsc)
table(sex)
table(cla)
table(disc)
table(disc, exclude = "放任") #可用exclude排除不要的次數統計
table(sex, disc, dnn = c("性別", "管教")) #dnn參數界定變數名稱
addmargins(table(sex, disc, dnn = c("性別", "管教"))) #加入邊際加總次數
rowSums(table(sex, disc, dnn = c("性別", "管教")))
colSums(table(sex, disc, dnn = c("性別", "管教")))
tsd <- table(sex, disc, dnn = c("性別", "管教"))
round(prop.table(tsd, 1), 3) #參數1表示計算橫列次數百分比
round(prop.table(tsd, 2), 3) #參數2表示計算直行次數百分比

tsd <- table(sex, disc, dnn = c("性別", "管教"))
cro.table <- round(prop.table(tsd, 1), 3) * 100
print(cro.table)

prop.table(margin.table(tsd, 2)) #margin.table()只會加總，不會算百分比
prop.table(margin.table(tsd, 1))
table(sex, disc, cla, dnn = c("性別", "管教", "班級"))
addmargins(table(sex, disc, cla, dnn = c("性別", "管教", "班級")))

sdc <- addmargins(table(sex, disc, cla, dnn = c("性別", "管教", "班級")))
round(prop.table(sdc), 3)


library(Hmisc)
describe(moti) #Hmisc套件的describe()也有類似table()的功能

describe(dsc[, 2:4])

describe(dsc[10])
describe(dsc[, 10])
describe(n.pedu) #3種方法都一樣

describe(dsc[5:6])


##xtabs()函數的應用==============================

xtabs(~ sex + disc, data = dsc)
xtabs(~ sex + disc, data = dsc, subset = cla == 1)
xtabs(~ sex + disc, data = dsc, subset = cla == 2)
xtabs(~ sex + disc, data = dsc, subset = cla == 3)
xtabs(~ sex + disc, data = dsc, subset = (cla = 1))
xtabs(~ sex + disc, data = dsc, subset = cla %in% c(1, 2)) #前數值 %in% 數值範圍
xtabs(~ sex + disc, data = dsc, subset = cla < 3)

table(moti)
barplot(table(moti))
barplot(table(moti), main = "學習動機長條圖", xlab = "學習動機", ylab = "次數", col = "green", border = "blue", density = 40)
labels = c("很低", "低", "普通", "高", "很高")
barplot(table(moti), main = "學習動機長條圖", xlab = "學習動機", ylab = "次數", names = labels, col = "green", border = "blue", density = 20)

f.col <- c(0, 1)[as.numeric(table(moti) > 25) + 1] #超過25的畫黑色
barplot(table(moti), col = f.col)

library(lattice)
f.col <- c(0, 1)[as.numeric(table(moti) >= 20) + 1]
barchart(table(moti), col = f.col, ylab = "moti")

plot(table(moti), type = "b", col = 1, ylim = c(0, 35)) #折線圖
text(table(moti), labels = c(15, 18, 30, 20, 17), pos = 1) #text為低階繪圖函數，要依附在高階繪圖函數中，pos為原始圖示點位移的位置量

varx = c(1:5) #繪製完整折線圖
vary = table(moti)
plot(vary, type = "b", col = 1, cex = 2, lwd = 2, ylim = c(0, 36), xaxt = "n")
lines(varx, vary, lwd = 2)
text(varx, vary, labels = vary, pos = 1)
axis(1, at = c(1:5), labels = c("很低", "低", "普通", "高", "很高"))


pie(table(moti))
pie(table(moti), main = "學習動機圓餅圖", border = "blue", col = c("white", "red", "green", "yellow", "black"), labels = c("很低", "低", "普通", "高", "很高"))
nlab = c("很低", "低", "普通", "高", "很高")
col_c = c("white", "red", "green", "yellow", "black")
pie(table(moti), main = "學習動機圓餅圖", border = "blue", label = nlab, col = col_c)


lab_1 <- c("很低", "低", "普通", "高", "很高")
col_c <- c(0:4)
per.moti <- (table(moti) * 100) / length(moti)
lab_2 <- paste(lab_1, "(", per.moti, "%)") #加入百分比
pie(table(moti), border = "blue", label = lab_2, col = col_c)


install.packages("plotrix") #準備畫3D圓餅圖
library(plotrix)
lab_1 <- c("very low", "low", "normal", "high", "very high")
col_c <- c("white", "gray80", "gray60", "gray50", "black")
per.moti <- (table(moti) * 100) / length(moti)
lab_2 <- paste(lab_1, "(", per.moti, "%)")
pie3D(table(moti), labels = lab_2, col = col_c, labelcex = .8, explode = .04) #labelcex標記文字的大小，explode立體扇形分開的程度


#描述性統計量==============================

library(fBasics)
attach(dsc)
head(dsc)
round(basicStats(msco), 2)
quantile(msco, .50)
quantile(msco, .25)
quantile(msco, .75)

round(basicStats(dsc[, 5:7]), 2)
round(basicStats(dsc[5:7]), 2)

fivenum(msco)
fivenum(esco)
fivenum(moti)

summary(msco)
summary(esco)
summary(moti)

summary(dsc[, 5:7])
names(dsc)

sapply(dsc, mean) #資料必須是計量變數
round(sapply(dsc, sd), 2)
dsc_part <- dsc[c(4:7)]
head(dsc_part)

sapply(dsc_part, range)
lapply(dsc_part, mean) #和sapply()函數的功能類似，只是回傳的屬性與格式不同
lapply(dsc_part, range)


##tapply()函數與aggregate()函數的應用==============================

round(tapply(msco, factor(cla), mean), 2) #tapply可以求出間斷變數各水準類別在計量變數的平均值
round(tapply(msco, disc, mean), 2)
round(tapply(msco, sex, mean), 2)
round(tapply(msco, factor(disc, levels = c("民主", "權威", "放任")), mean), 2)
round(tapply(msco, disc, levels = c("民主", "權威", "放任"), mean), 2)

round(tapply(msco, factor(cla), sd), 2)
round(tapply(msco, factor(disc, levels = c("民主", "權威", "放任")), sd), 2)
round(tapply(msco, factor(disc), sd), 2)
round(tapply(msco, factor(sex), sd), 2)
round(tapply(msco, sex, sd), 2)
dsc$sex <- recode(sex, "'Male' = '男生'; 'Female' = '女生'")
round(tapply(msco, factor(sex, levels = c("男生", "女生")), sd), 2)
round(tapply(msco, factor(cla, levels = c(1, 2, 3)), min), 2)
round(tapply(msco, factor(cla, levels = c(1, 2, 3)), max), 2)
round(tapply(msco, factor(cla, levels = c(1, 2, 3)), length), 2)

round(tapply((msco + esco) / 2, factor(cla, levels = c(1, 2, 3)), mean), 2)
round(tapply((msco + esco) / 2, factor(cla), mean), 2) #引述可以為變數間的數學運算式
tapply(msco, n.pedu, function(x) c(sum(x), mean(x), sd(x), length(x)))
esco.disc <- tapply(esco, disc, function(x) c(sum(x), mean(x), sd(x), length(x)))
esco.disc[1]

cm1 <- msco[which(cla == 1)]
cm2 <- msco[which(cla == 2)]
cm3 <- msco[which(cla == 3)]
mean(cm1)
mean(cm2)
mean(cm3)

mean(msco[which(cla == 1)])
mean(msco[which(cla == 2)])
mean(msco[which(cla == 3)])
sd(msco[which(cla == 1)])
sd(msco[which(cla == 2)])
sd(msco[which(cla == 3)])

table(moti[which(cla == 1)])
table(moti[which(cla == 2)])
table(moti[which(cla == 3)])

by(dsc$msco, sex, summary)
by(dsc$esco, disc, summary)

aggregate(dsc[, c(5:6)], by = list(dsc[, 3]), mean)
aggregate(dsc[, c(5:6)], by = list(dsc[, "sex"]), mean)
aggregate(dsc[, c(5:6)], by = list(dsc[, 8]), mean)
aggregate(dsc[, c(5:6)], by = list(dsc[, "disc"]), mean)
aggregate(dsc[, c(5:6)], by = list(dsc[, "disc"]), sd)
aggregate(dsc[, c(5:6)], by = list(dsc[, 8]), sd)
aggregate(dsc[, c(5:6)], by = list(dsc[, 8]), function(x) c(mean(x), sd(x), length(x)))
aggregate(dsc[, c(5:6)], by = list(dsc[, "sex"]), function(x) c(mean(x), sd(x), length(x)))
aggregate(dsc[, c(5:6)], by = list(dsc[, 2], dsc[, 3]), function(x) c(mean(x), sd(x)))
aggregate(dsc[, c(5:6)], by = list(dsc[, "disc"], dsc[, "sex"]), function(x) c(mean(x), sd(x)))


##標準分數==============================

dsc_c1 <- dsc[dsc$cla == "1", ]
head(dsc_c1)
dsc_c2 <- dsc[dsc$cla == "2", ]
head(dsc_c2)
dsc_c3 <- dsc[dsc$cla == "3", ]
head(dsc_c3)

dsc_c1 <- dsc[dsc$cla == "1" & dsc$sex == "女生", ]
attach(dsc)
dsc_12 <- dsc[cla == 1 & sex == "女生", ]
head(dsc_12)
dsc_3m <- dsc[cla == 3 & msco > 90, ]
head(dsc_3m)


dsc_c1 <- dsc[dsc$cla == "1", ]
dsc_c1$zms1 <- (dsc_c1$msco - mean(dsc_c1$msco)) / sd(dsc_c1$msco)
dsc_c1$zms1 <- round(dsc_c1$zms1, 3)
head(dsc_c1)

dsc_c1 <- dsc[dsc$cla == "1", ]
zms2 <- round(scale(dsc_c1$msco), 3) #用scale()求Z分數
head(zms2)

dsc_c1$tms1 <- dsc_c1$zms1 * 10 + 50
head(dsc_c1)
round(sd(dsc_c1$zms1), 2)
round(sd(dsc_c1$tms1), 2)
round(mean(dsc_c1$zms1), 2)
round(mean(dsc_c1$tms1), 2)

dsc <- read.csv("sc_1.csv", header = TRUE)
attach(dsc)
dsc_c1 <- dsc[dsc$cla == 1, ]
dsc_c2 <- dsc[dsc$cla == 2, ]
dsc_c3 <- dsc[dsc$cla == 3, ]

dsc_c1$zms1 <- (dsc_c1$msco - mean(dsc_c1$msco)) / sd(dsc_c1$msco)
dsc_c1$zms1 <- round(dsc_c1$zms1, 3)

dsc_c2$zms1 <- (dsc_c2$msco - mean(dsc_c2$msco)) / sd(dsc_c2$msco)
dsc_c2$zms1 <- round(dsc_c2$zms1, 3)

dsc_c3$zms1 <- (dsc_c3$msco - mean(dsc_c3$msco)) / sd(dsc_c3$msco)
dsc_c3$zms1 <- round(dsc_c3$zms1, 3)

dsc_zm <- rbind(dsc_c1, dsc_c2, dsc_c3)
dsc_zm$tms1 <- dsc_zm$zms1 * 10 + 50

dsc_c1 <- dsc[dsc$cla == 1, ]
dsc_c2 <- dsc[dsc$cla == 2, ]
dsc_c3 <- dsc[dsc$cla == 3, ]

dsc_c1$zms1 <- round(scale(dsc_c1$msco), 3)
dsc_c2$zms1 <- round(scale(dsc_c2$msco), 3)
dsc_c3$zms1 <- round(scale(dsc_c3$msco), 3)
dsc_zm <- rbind(dsc_c1, dsc_c2, dsc_c3)
dsc_zm$tms1 <- dsc_zm$zms1 * 10 + 50
print.data.frame(dsc_zm)

dsc <- read.csv("sc_1.csv", header = TRUE)
attach(dsc)
for(i in 1:3) {
        dsc$zms[cla == i] <- round(scale(msco[cla == i]), 2)
        dsc$tms[cla == i] <- round(scale(msco[cla == i]), 2) * 10 + 50
}
print.data.frame(dsc)

dsc <- read.csv("sc_1.csv", header = TRUE)
attach(dsc)
gronum <- 3
for(i in 1:gronum) {
        dsc$zms[cla == i] <- round(scale(msco[cla == i]), 2)
        dsc$tms[cla == i] <- round(scale(msco[cla == i]), 2) * 10 + 50
}
print.data.frame(dsc)

dscore <- read.csv("tenscore.csv", header = TRUE)
attach(dscore)
names(dscore)
rm(score)
tapply(score, factor(class), mean)
round(tapply(score, factor(class), sd), 2)

gronum <- 10
for(i in 1:gronum) {
        dscore$zscore[class == i] <- round(scale(score[class == i]), 3)
        dscore$tscore[class == i] <- round(scale(score[class == i]), 3) * 10 + 50
}
print.data.frame(dscore)

attach(dscore)
round(tapply(zscore, factor(class), mean), 3)
round(tapply(zscore, factor(class), sd), 3)

round(tapply(tscore, factor(class), mean), 2)
round(tapply(tscore, factor(class), sd), 2)


##計量變數的圖形==============================

###直方圖==============================

hist(x, 
     breaks = , #數值向量分割點
     freq = NULL, #TRUE表示次數
     probability = !freq, #TRUE表示機率
     include.lowest = TRUE, 
     right = TRUE, #右端封閉左端開放
     density = NULL, #陰影密度
     angle = 45, 
     col = NULL, #長條的顏色
     border = NULL, #長條邊框顏色
     main = "圖標題名稱", 
     xlim = range(breaks), 
     ylim = NULL, 
     xlab = xname, 
     ylab, 
     axes = TRUE, #繪製軸線
     plot = TRUE, #繪製直方圖，否則只呈現區間次數
     labels = FALSE) #邏輯或文字設定

hist(msco, col = "green", border = "blue", density = 60)
hist(msco, breaks = 10, col = "green", border = "blue", density = 60)
hist(msco, breaks = 10, prob = TRUE)
lines(density(msco), col = "black", lwd = 2.0) #寬度2的黑線

str(hist(esco, labels = TRUE, col = 3, border = 1, density = 60, main = NULL, ylim = c(0, 25))) #str()函數可以增列各長條對應的數值（次數），ylim界定y軸數值的上下限
str(hist(moti, labels = TRUE, density = 20, main = "", ylim = c(0, 35)))


###盒形圖==============================

boxplot(esco, col = "gray") #horizontal = T可繪製水平盒形圖，超過1.5倍盒長為偏離值用圈圈表示，超過3倍為極端值用星星表示
boxplot(dsc[, 5:6], col = c("gray", "green"))
boxplot(dsc[, 5:6], col = c(8, 3))
boxplot(dsc[, 5:6], horizontal = TRUE, col = c("gray", "green"))


library(lattice)
bwplot(esco) #用lattice套件的bwplot()函數繪製水平盒形圖
library(lattice)
bwplot(esco ~ factor(sex))


###常態機率圖==============================

qqnorm(msco, ylab = "m_score") #若樣本為常態分配，理論上會分佈在常態累積機率直線圖上，但只是初步檢核而已
qqline(msco)

qqnorm(esco, ylab = "score")
qqline(esco)

install.packages("normtest")
library(normtest) #進行常態性分佈之峰態與偏態之檢定

skewness.norm.test(msco) #常態性偏態檢定，偏態係數是否等於0的檢定
skewness.norm.test(esco) #T值>0為右偏（正偏），<0為左偏（負偏）

kurtosis.norm.test(msco) #常態性峰態檢定，峰態係數是否等於0的檢定，等於0為常態峰
kurtosis.norm.test(esco) #T值>0為高狹峰


install.packages("nortest")
library(nortest) #進行常態性檢定

ad.test(msco) #觀察值個數最少為7，Anderson-Darling
cvm.test(msco) #觀察值個數最少為7，Cramer-von Mises

lillie.test(msco) #觀察值個數最少為4，Lilliefors(Kolmogorov-Smirnov)

pearson.test(msco) #Pearson chi-square

sf.test(msco) #觀察值個數在5至5000之間，Shapiro-Francia
shapiro.test(msco) #R基本套件中提供的函數

ad.test(esco)
cvm.test(esco)
lillie.test(esco)
pearson.test(esco) #母數統計法中，許多統計分析對於計量變數違反常態性假設時也有很高的強韌性，所以進行常態性檢定時，可以將顯著水準訂嚴格一點
sf.test(esco) #若有些達到顯著水準，有些未達，可綜合判斷，若達顯著水準的統計量較多，可做出拒絕虛無假設的結論


#卡方檢定==============================


##適合度考驗==============================

it01 <- c("非常同意" = 30, "同意" = 70, "無意見" = 32, "不同意" = 28, "非常不同意" = 40)
chisq.test(it01) #卡方適合度檢驗，變數只有一個，檢定變數各水準類別實際觀察次數與理論期望次數的差異

exp <- 40
chi.sta <- (30 - exp)^2 / exp + (70 - exp)^2 / exp + (32 - exp)^2 / exp + (28 - exp)^2 / exp + (40 - exp)^2 / exp #卡方統計量的求法
print(chi.sta)

chisq.test(it01, p = c(.15, .25, .15, .25, .20)) #加入選項機率

it01 <- c("非常同意" = 30, "同意" = 70, "無意見" = 32, "不同意" = 28, "非常不同意" = 40)
ratio <- c(.15, .25, .15, .25, .20)
chisq.test(it01, p = ratio)


chi.sta <- (30 - 30) ^ 2 / 30 + (70 - 50) ^ 2 / 50 + (32 - 30) ^ 2 / 30 + (28 - 50) ^ 2 / 50 + (40 - 40) ^ 2 / 40
print(chi.sta)


dsc <- read.csv("dsc_reco.csv", header = TRUE)
attach(dsc)

tabulate(factor(moti)) #各選項次數
prop.table(table(moti)) #各選項百分比

chisq.test(tabulate(factor(moti))) #看各選項百分比是否有顯著不同
chisq.test(table(moti))

f.moti <- table(moti)
chisq.test(f.moti)

chisq.test(moti) #WARNING! 一定要先轉為因子變數！！！不然結果不對


##百分比同質性檢定==============================

chisq.test(sex, moti) #檢定變數有兩個，探究兩個間斷變數之類別間的關係，探究不同J個群體在I個反應的百分比是否相同
chisq.test(sex, factor(moti)) #最好還是要界定為因子變數

chisq.test(factor(moti), sex) #順序不同沒關係

addmargins(table(sex, moti))

round(prop.table(table(sex, moti), 1), 3)

chisq.test(sex, disc)

addmargins(table(sex, disc))
round(prop.table(table(sex, disc), 1), 3)


##卡方獨立性檢定==============================

dchi <- read.csv("chi_sq.csv", header = TRUE)
tail(dchi)
attach(dchi)
length(clot)
library(car)
dchi$year <- recode(year, "1 = '一年級'; 2 = '二年級'; 3 = '三年級'")
dchi$clot <- recode(clot, "1 = '甲款式'; 2 = '乙款式'; 3 = '丙款式'")
dchi$hat <- recode(hat, "1 = '黑色'; 2 = '白色'; 3 = '紅色'")

head(dchi)
tail(dchi)

attach(dchi) #框架物件dchi的內容已被修改，要重新attach()！！！
table(clot)
round(prop.table(table(clot)), 3)
table(hat)
prop.table(table(hat))

chisq.test(table(clot)) #卡方適合度考驗
chisq.test(table(hat))

addmargins(table(year, clot))
round(prop.table(table(year, clot), 1), 3)

addmargins(table(year, hat))
round(prop.table(table(year, hat), 1), 3)

chisq.test(year, clot) #卡方獨立性檢定

sqrt(13.866 / (13.866 + 120))
chisq.test(year, hat)
chisq.test(hat, year) #對調結果相同

sqrt(15.718 / (15.718 + 120)) #方形列聯表，統計量可以用列聯係數(C)表示


##2X2列聯表分析==============================

twot <- read.csv("t_test_1.csv", header = T)
attach(twot)
addmargins(table(hom, pass))
fisher.test(hom, pass, alternative = "two.sided") #2X2列聯表採用Fisher's精確性檢定

addmargins(table(pass, sex))
fisher.test(pass, sex, alternative = "two.sided")


#單一樣本檢定==============================

##單一樣本t檢定==============================

o_test = read.csv("o_test_1.csv", header = T)
attach(o_test)
head(o_test)
length(num)
is.na(o_test)[2:4]
mean(reti)
mean(spti)
mean(wrsc)
for(i in 2:4)
        print(mean(o_test[[i]]))
for(i in 2:4)
        print(sd(o_test[[i]]))

t.test(reti, mu = 50, alternative = "two.sided") #雙尾
t.value <- (mean(reti) - 50) / (sd(reti) / sqrt(length(reti)))
print(t.value)

t.test(reti, mu = 50, alternative = "greater") #單尾右側

t.test(spti, mu = 55, alternative = "greater")
t.value <- (mean(spti) - 55) / (sd(spti) / sqrt(length(spti)))
print(t.value)

t.test(spti, mu = 55, alternative = "two.sided")

t.test(wrsc, mu = 66, alternative = "two.sided")
t.value <- (mean(wrsc) - 66) / (sd(wrsc) / sqrt(length(wrsc)))

t.test(wrsc, mu = 66, alternative = "less")

2 * pt(-2.038, 19) #雙尾檢定p值
pt(-2.038, 19) #單尾檢定p值

qt(.95, 19) #單尾檢定臨界t值
qt(.975, 19) #雙尾檢定臨界t值


##單一母體比例（百分比）檢定==============================

binom.test(x = 15, n = 40, p = .25, alternative = "two.sided") #40位通過15位的比例與母體比例.25是否有差別
binom.test(x = 15, n = 40, p = 1/4, alternative = "t")
binom.test(x = 15/40, n = 1, p = .25, alternative = "two.sided") #Error

prot.test(x = 15, n = 40, p = .25, alternative = "two.sided")
prop.test(x = 15/40, n = 1, p = .25, alternative = "two.sided") #Warning

binom.test(x = 15, n = 40, p = .25, alternative = "g")

binom.test(x = 7, n = 50, p = .25, alternative = "t")
prop.test(x = 7, n = 50, p = 1/4, alternative = "t")

binom.test(x = 7, n = 50, p = 1/4, alternative = "l")


#雙樣本檢定==============================

##不同性別學生每天閱讀時間差異比較==============================

twot <- read.csv("t_test_1.csv", header = T)
attach(twot)
head(twot)

twot1 <- data.frame(twot)
require(car)
twot1$sex <- recode(twot1$sex, "1 = 'boy'; 2 = 'girl'")
twot1$hom <- recode(twot1$hom, "1 = 'complete'; 2 = 'single'")
attach(twot1)
head(twot1, 5)

names(twot)
names(twot1)

var.test(reti ~ sex, alternative = "t") #變異數同質性檢定
bartlett.test(reti ~ sex) #兩個或多個群體之變異數同質性檢定也可以用bartlett.test()

library(car)
leveneTest(reti, sex)
leveneTest(reti, factor(sex)) #最好明確界定變數為因子類型

t.test(reti ~ sex, mu = 0, alternative = "t", var.equal = T) #進行獨立樣本t檢定
t.test(reti ~ sex, mu = 0, alternative = "t", var.equal = T, paired = F, conf.level = .95)

with(twot1, {t.test(reti ~ sex, mu = 0, alternative = "t", var.equal = T)})
with(twot, t.test(reti[sex == 1], reti[sex == 2]), var.equal = T) #也可以用群組factor level的邏輯判斷式

round(tapply(reti, sex, mean), 3)
round(tapply(reti, sex, sd), 3)

round(tapply(twot1$reti, twot1$sex, mean), 3)
round(tapply(twot1$reti, twot1$sex, sd), 3)

mg <- subset(twot, sex == 1)
fg <- subset(twot, sex == 2)
et <- matrix(NA, nrow = 2, ncol = 5)
colnames(et) <- c("Mean", "Sd", "N", "MeanDif", "SE")
rownames(et) <- c("Male", "Female")
j = 4; i = 1
et[i, 1] <- mean(mg[, j])
et[i, 2] <- sd(mg[, j])
et[i, 3] <- length(mg[, j])
et[i, 4] <- mean(mg[, j]) - mean(fg[, j])
et[i, 5] <- sd(mg[, j]) / (sqrt(length(mg[, j])))
et[i + 1, 1] <- mean(fg[, j])
et[i + 1, 2] <- sd(fg[, j])
et[i + 1, 3] <- length(fg[, j])
et[i + 1, 4] <- mean(mg[, j]) - mean(fg[, j])
et[i + 1, 5] <- sd(fg[, j]) / (sqrt(length(fg[, j])))
round(et, 3)

mean(reti)
mean(reti[sex == 1])
mean(reti[sex == 2])

mean(spti)
mean(spti[sex == 1])
mean(spti[sex == 2])

mean(wrsc)
mean(wrsc[sex == 1])
mean(wrsc[sex == 2])

tdes <- matrix(NA, nrow = 2, ncol = 5)
colnames(tdes) <- c("Mean", "Sd", "N", "MeanDif", "SE")
rownames(tdes) <- c("Male", "Female")
for(i in 1:2) {
        tdes[i, 1] <- round(mean(reti[sex == i]), 2)
        tdes[i, 2] <- round(sd(reti[sex == i]), 2)
        tdes[i, 3] <- length(reti[sex == i])
        tdes[i, 4] <- round(mean(reti[sex == 1]) - mean(reti[sex == 2]), 2)
        tdes[i, 5] <- round(sd(reti[sex == i]) / sqrt(length(reti[sex == i])), 2)
}
print(tdes)

for(j in 4:6) {
        for(i in 1:2) {
              tdes[i, 1] <- round(mean(twot[, j][sex == i]), 2)
              tdes[i, 2] <- round(mean(twot[, j][sex == i]), 2)
              tdes[i, 3] <- length(twot[, j][sex == i])
              tdes[i, 4] <- round(mean(twot[, j][sex == 1]), 2) - round(mean(twot[, j][sex == 2]), 2)
              tdes[i, 5] <- round(sd(twot[, j][sex == i]) / (sqrt(length(twot[, j][sex == i]))), 2)
        }
        print(tdes)
}

for(j in 4:6) {
        print(paste("dependent variable: ", names(twot[j]), " descriptive statistics"))
        for(i in 1:2) {
                tdes[i, 1] <- round(mean(twot[, j][sex == i]), 2)
                tdes[i, 2] <- round(sd(twot[, j][sex == i]), 2)
                tdes[i, 3] <- length(twot[, j][sex == i])
                tdes[i, 4] <- round(mean(twot[, j][sex == 1]), 2) - round(mean(twot[, j][sex == 2]), 2)
                tdes[i, 5] <- round(sd(twot[, j][sex == i]) / (sqrt(length(twot[, j][sex == i]))), 2)
        }
        print(tdes)
}


##不同性別學生每天運動時間差異比較==============================

var.test(spti ~ sex, data = twot, alternative = "t")
leveneTest(spti, factor(sex))

t.test(spti ~ sex, mu = 0, alternative = "t", var.equal = T)

var.test(reti ~ hom, alternative = "t")
var.test(spti ~ hom, alternative = "t")
var.test(wrsc ~ hom, alternative = "t")

for(i in 4:6) {
        print(paste("DEvar:", names(twot[i])))
        print(var.test(twot[[i]] ~ hom, alternative = "t"))
}

for(i in 4:6) {
        print(t.test(twot[[i]] ~ hom, alternative = "t", var.equal = T))
}


group <- rep(c("treat", "control"), each = 10)
score <- c(7, 6, 7, 9, 10, 8, 6, 4, 5, 3, 6, 5, 4, 7, 8, 2, 4, 5, 3, 3)
var.test(score ~ group, alternative = "t")
t.test(score ~ group, alternative = "t", var.equal = T)

t.test(score ~ group, alternative = "l", var.equal = T)


##相依樣本檢定==============================

t.test(typa, typb, mu = 0, alternative = "t", var.equal = T, paired = T, conf.level = .95)
t.test(typb, typc, mu = 0, alternative = "t", var.equal = T, paired = T, conf.level = .95)
t.test(typa, typc, mu = 0, alternative = "t", var.equal = T, paired = T, conf.level = .95)

print(mean(typa))
print(mean(typb))
print(mean(typc))

print(sd(typa))
print(sd(typb))
print(sd(typc))


##兩群體比例檢定==============================

pass <- c(8, 12) #通過人數
scho <- c(22, 26) #參加人數
prop.test(pass, scho)


nopass <- c(5, 15) #2間學校通過人數
class <- c(16, 20) #2間學校參加人數
prop.test(nopass, class) #2間學校通過率比較

pass <- c(20, 15, 18) #3間
school <- c(40, 50, 25)
prop.test(pass, school) #至少有一間的通過率跟其他間不同


##無母數檢定==============================

attach(twot)
kruskal.test(reti, hom, data = twot) #Mann-Whitney U
kruskal.test(wrsc, sex, data = twot) #群組間計量變數之分配的位置參數（location parameters）是否相同，虛無假設為兩群體的等級和沒有顯著不同
wilcox.test(reti ~ hom, mu = 0, paired = FALSE, alternative = "t") #兩變數的中位數是否相同（Wilcoxon signed rank）

wilcox.test(wrsc ~ sex, mu = 0, paired = FALSE, alternative = "t")

wilcox.test(typa, typb, paired = TRUE, alternative = "t")
wilcox.test(typa, typc, paired = TRUE, alternative = "t")

ks.test(reti, spti) #檢定兩計量變數是否在相同的連續型機率分配（Kolmogorov-Smirnov）

locationTest(typa, typb, method = "kw2") #kw2表示界定檢定的量數為兩變數的中位數
locationTest(typa, typc, method = "kw2") #檢定兩變數是否有相同的中位數或相同的集中趨勢
locationTest(typa, typb, method = "t") #改為t表示檢定的是兩變數的平均數的差異

locationTest(typa, typc, method = "t")


#相關分析==============================

dco <- read.csv("cor_1.csv", header = T, fileEncoding = "BIG-5") #csv內含中文時，fileEncoding = "BIG-5"！！！)
head(dco)

is.na(dco[, 5:8])
all(is.na(dco[, 5:8]) == FALSE)
head(dco[, 5:8])
tail(dsc[, 5:8])


##Pearson積差相關==============================

attach(dco)
cor.test(manx, moti, alternative = "t", method = "pearson")

cov(manx, moti)
sd(manx)
sd(moti)
cov(manx, moti) / (sd(manx) * sd(moti))

library(fBasics)
correlationTest(manx, moti, "pearson")
pearsonTest(manx, moti) #fBasics套件中的pearsonTest()也可以進行相關分析，單尾左側檢定的p值等於1減掉雙尾檢定的p值除以2


cor.test(manx, matt, alternative = "t", method = "pearson")
correlationTest(manx, matt, "pearson")

summary(pearsonTest(manx, moti)) #pearsonTest()函數分類為fHTEST，模式型態為S4


cor.test(manx, msco, alternative = "t", method = "pearson") #0.7以上為高度，0.4以下為低度，之間稱為中度
cor.test(moti, matt, alternative = "t", method = "pearson")
cor.test(moti, msco, alternative = "t", method = "pearson")
cor.test(matt, msco, alternative = "t", method = "pearson")


dco_1 <- dco[-c(1:4)]
names(dco_1)
head(dco_1, 3)

dco_2 <- dco[c(5:8)]
head(dco_2, 3)
round(cor(dco_1), 2) #cor()求3個以上計量變數間的相關矩陣
cor(dco) #Error，若變數中的因子變數水準數值標記為文字群組，會出現錯誤

attach(dco)
names(dco)
cor4 <- array(, c(4, 4)) #空的陣列
nu <- 4
for(i in 1:4) 
        for(j in 1:4) 
                cor4[i, j] <- round(cor(dco[, i + nu], dco[, j + nu]) , 3)
                colnames(cor4) <- names(dco[(nu + 1):(nu + 4)])
                rownames(cor4) <- names(dco[(nu + 1):(nu + 4)])
print(cor4)


cor4 <- array(NA, c(4, 4))
nu <- 4
for(i in 1:4) 
        for(j in 1:4)
                cor4[i, j] <- round(cor(dco[, i + nu], dco[, j + nu]), 3)
                colnames(cor4) <- paste(colnames(dco[-c(1:4)]))
                rownames(cor4) <- paste(colnames(dco[-c(1:4)]))
print(cor4)


cor4 <- array(, c(4, 4))
nu = 4
for(i in 1:4) 
        for(j in 1:4)
                if(i < j) {
                        cor4[i, j] <- round(cor(dco[, i + nu], dco[, j + nu])^2, 3)
                } else {
                        cor4[i, j] <- round(cor(dco[, i + nu], dco[, j + nu]), 3)
                }
                colnames(cor4) <- names(dco[(nu + 1):(nu + 4)])
                rownames(cor4) <- names(dco[(nu + 1):(nu + 4)])
print(cor4) #右上角為決定係數，左下角為相關係數


iq <- c(100, 105, 107, 108, 109, 104, 110, 111, 106, 112, 114, 116)
score <- c(60, 65, 72, 71, 66, 74, 78, 81, 86, 88, 74, 80)
require(fBasics)
correlationTest(iq, score, "pearson")


##等級相關==============================

cor.test(matt, msco, alternative = "t", method = "kendall") #變數尺度為次序變數，兩變數間的相關參數估計值為等級相關係數，使用Kendall等級相關統計量
cor.test(matt, msco, alternative = "t", method = "spearman") #Warning的意思是兩變數中，變數內的值有相同者，而無法估算精確的顯著性，雖無法精確估算p值，但p值仍有很高的強韌性（偏誤值很小）

spearmanTest(matt, msco) #fBasics套件中的函數spearmanTest()進行等級相關分析
kendallTest(x = matt, y = msco) #fBasics套件中的函數kendallTest()進行等級相關分析


X <- c(50, 70, 60, 40, 30, 20, 10)
Y <- c(26, 67, 33, 45, 34, 23, 11)

cor.test(X, Y, alternative = "t", method = "spearman")
cor.test(X, Y, alternative = "t", method = "kendall")

tempx <- c(50, 70, 60, 40, 30, 20, 10)
varx <- length(tempx) + 1 - rank(tempx) #數值向量轉換等級
print(varx)
tempy <- c(26, 67, 33, 45, 34, 23, 11)
vary <- length(tempy) + 1 - rank(tempy)
print(vary)

cor.test(varx, vary, alternative = "t", method = "spearman")
cor.test(varx, vary, alternative = "t", method = "kendall") #使用原始測量值或轉換等級，兩者得出的結果是相同的


##兩個變數間散佈圖的繪製==============================

plot(1:7, 1:7, cex = 1:7, pch = 0:6, col = 1:7)

plot(1:7, 1:7, cex = 1:7, pch = 0:6, col = 1:7, xlim = c(0, 9), ylim = c(0, 9))
text(1:7, 1:7, labels = paste(0:6), cex = 1:7, col = 1:7, xlim = c(0, 9), ylim = c(0, 9))
plot(matt, msco, col = 4, type = "p") #p點l線c只有線b點與線o強化點與線h垂直線s階梯線n不畫圖

plot(matt, msco, type = "n")
grid() #增加散佈圖的格線 #lty = 1 "solid"實線; = 2 "dashed"短線條; = 3 "dotted"點線條; = 4 "dotdash"點短線; = 5 "longdash"長短線; = 6 "twodash"雙短線
points(matt, msco, pch = 1, lwd = 2) #lwd決定符號線條粗細

plot(matt, msco, pch = home)
grid(lty = 2)

plot(matt, msco, pch = as.integer(home))
grid(lty = 2)

plot(matt, msco, pch = home, lwd = 2)
par(family = "Hannotate TC Regular") #字體查詢：Finder -> 應用程式 -> 字體簿
labnames <- c("完整家庭", "單親家庭", "隔代教養")
grid(lty = 2)
legend(20, 50, labnames, pch = c(1:3))

plot(manx, matt, col = 4)
plot(manx, moti, col = 4, type = "b")

qplot(x, y, data, size = colour = , xlim =, ylim =, main =, xalb =, ylab =, asp = NA) #ggplot2套件中的函數
library(ggplot2)
qplot(matt, msco, size = 4, colour = 2)

library(ggplot2)
qplot(matt, manx, size = 4)

varnum <- c(5:8)
m_cor <- cor(dco[varnum])
library(ellipse)
plotcorr(m_cor, col = "black") #繪製矩陣中所有配對變數的相關圖 #越接近圓形表示零相關

plotcorr(m_cor, col = "black", diag = T, type = "lower") #type = "lower"只畫左下角，diag = T表示要畫對角線的圖示


#單因子變異數分析==============================

##整體檢定==============================

###匯入資料檔==============================

dco <- read.csv("cor_1.csv", header = TRUE, fileEncoding = "BIG-5")
attach(dco)
tail(dco)

dco_new <- data.frame(dco)
attach(dco_new)

library(car)
dco_new$year <- recode(year, "1 = '一年級'; 2 = '二年級'; 3 = '三年級'; 4 = '四年級'")
dco_new$home <- recode(home, "1 = '完整家庭'; 2 = '單親家庭'; 3 = '隔代教養'")
head(dco_new)
tail(dco_new)

round(tapply(manx, list(factor(year)), mean), 2)
round(tapply(manx, list(factor(year)), sd), 2)

objects() #查看目前物件

with(dco_new, {
        print(round(tapply(manx, list(year), mean), 2))
        print(round(tapply(manx, list(year), sd), 2))
})


###aov()函數的應用==============================

ymanx <- aov(manx ~ factor(year))
ymanx
summary(ymanx)

ymoti <- aov(moti ~ factor(year))
summary(ymoti)

model.tables(ymanx, type = "effects")
model.tables(ymanx, type = "means")


###anova()函數的應用==============================

anova(ymanx)
coef(ymanx)
effects(ymanx)
residuals(ymanx)
summary(ymanx)
fitted.values(ymanx)
fitted.values(ymanx)

library(xtable) #得到可以放網頁的程式碼
print(xtable(summary(ymanx)), type = "html")
ymanx_t <- xtable(ymanx)
print(ymanx_t, type = "html")
print(xtable(ymanx), type = "html")


##事後比較==============================

###使用基本套件函數==============================

TukeyHSD(ymanx)

new.ymanx <- aov(dco_new$manx ~ dco_new$year)
summary(new.ymanx)
TukeyHSD(new.ymanx)

install.packages("asbio")
library(asbio)

pairw.anova(manx, factor(year), conf.level = .95, method = "tukey") #會輸出決策欄
pairw.anova(manx, factor(year), conf.level = .95, method = "lsd") #LSD法（最小顯著差異法）
pairw.anova(manx, factor(year), conf.level = .95, method = "scheffe") #Scheffe法（雪費法）

dunnettCI(manx, factor(year), conf.level = .95, control = 1) #Dunn法
dunnettCI(manx, factor(year), conf.level = .95, control = 2)

scheffeCI(manx, factor(year), conf.level = .95) #Scheffe法
bonfCI(manx, factor(year), conf.level = .95) #Bonferroni法


##不同地區在數學態度的差異比較==============================

loatt <- aov(matt ~ area)
summary(loatt)
library(asbio)
scheffeCI(matt, factor(area), conf.level = .95)
pairw.anova(matt, factor(area), conf.level = .95, method = "tukey")


##R編輯器命令稿在ANOVA的應用==============================

colnames(dco_new)
attach(dco_new)
for(i in 5:8) {
        ano.model <- aov(dco_new[[i]] ~ dco_new[[2]])
        print(summary(ano.model))
        print(TukeyHSD(ano.model))
}

attach(dco_new)
library(asbio)
for(i in 2:4) {
        for(j in 5:8) {
                ano.model <- aov(dco_new[[j]] ~ dco_new[[i]])
                print(summary(ano.model))
                print(scheffeCI(dco_new[[j]], factor(dco_new[[i]])))
        }
}

for(i in 2:4) {
        for(j in 5:8) {
                print(paste("varDE:", colnames(dco_new[j]), "---varID:", colnames(dco_new[i])))
                print(round(tapply(dco_new[[j]], list(factor(dco_new[[i]])), mean), 2))
                print(round(tapply(dco_new[[j]], list(factor(dco_new[[i]])), sd), 2))
        }
}
colnames(dco_new[3])
colnames(dco_new[5])


##單因子相依樣本變異數分析==============================

###量表內向度的差異比較==============================

ddim <- read.csv("dimen_anova.csv", header = TRUE, fileEncoding = "BIG-5")
attact(ddim)
ddim

temp <- read.csv("dim4.csv", header = TRUE, fileEncoding = "BIG-5")
temp

library(reshape)
temp_1 <- melt(temp, id.vars = "subject") #用melt()將資料重新排列
temp_1

names(temp_1) <- c("blocks", "dimen", "score")
print.data.frame(temp_1)

ddim_md <- aov(score ~ factor(blocks) + factor(dimen), data = ddim)
summary(ddim_md)

attach(ddim)
blocks <- as.factor(blocks)
dimen <- as.factor(dimen)
ddim_md_1 <- aov(score ~ blocks + dimen, data = ddim)
print(ddim_md_1)
summary(ddim_md_1)

TukeyHSD(ddim_md)

library(asbio)
scheffeCI(score, factor(dimen), conf.level = .95)

install.packages("laercio")
library(laercio)

LDuncan(ddim_md, "factor(dimen)", conf.level = .95) #用laercio套件的LDuncan()函數進行事後比較
LTukey(ddim_md, "factor(dimen)", conf.level = .95)


###重複量數設計==============================

blo <- read.csv("o_block.csv", header = TRUE, fileEncoding = "BIG-5")
attach(blo)
blo
round(tapply(anxi, list(factor(times)), mean), 2)
round(tapply(anxi, list(factor(blocks)), sd), 2)
round(tapply(anxi, list(factor(blocks)), mean), 2)

o_blo <- aov(anxi ~ factor(blocks) + factor(times), data = blo)
summary(o_blo)

TukeyHSD(o_blo)

library(asbio)
scheffeCI(anxi, factor(times), conf.level = .95)


#二因子變異數分析==============================

#典型相關==============================

#迴歸分析==============================

reg <- read.csv("reg_1.csv", header = TRUE, fileEncoding = "BIG-5")
attach(reg)
names(reg)
head(reg, 3)
tail(reg, 3)
length(rownames(reg)); nrow(reg) #一樣
length(colnames(reg))


##簡單迴歸分析==============================

###進行簡單迴歸分析假定的檢核==============================

reg_m <- lm(acad ~ moti, data = reg)

library(car)
ncvTest(reg_m) #檢定「等分散性」（變異數齊一性）
durbinWatsonTest(reg_m) #檢定殘差「獨立性」，虛無假設為殘差間無自我相關（誤差項彼此獨立） #檢定結果：違反
shapiro.test(residuals(reg_m)) #檢定殘差「常態性」 #檢定結果：違反
plot(reg_m, which = 1)
plot(reg_m, which = 2) #常態Q-Q圖資料多數落在對角線上表示殘差服從常態分配
plot(reg_m, which = 3)


###進行迴歸分析

reg_m <- lm(acad ~ moti, data = reg)
summary(reg_m)

reg_m <- aov(acad ~ moti, data = reg)
summary(reg_m)
deviance(reg_m) #精確估算SSE（誤差變異）
round((r_sq <- 14862 / (14862 + 14785)), 3) #R-squared
names(reg_m) #模型物件包含的次函數
round(coef(reg_m), 3) #迴歸模型的估計值係數

plot(acad ~ moti, type = "p", cex = 1.5)
abline(reg_m, lwd = 3, col = 4, lty = 1) #繪製迴歸線
grid(nx = 10, ny = 20)


##複迴歸==============================

reg_d <- reg[-c(1, 2)]
names(reg_d)
round(cor(reg_d), 2) #cor()求資料集的相關矩陣

library(ellipse)
plotcorr(cor(reg_d), col = 1)

library(corrplot)
corm <- cor(reg[, 3:9])
corrplot(corm, method = "number", col = "black", cl.pos = "n")
corrplot(corm, method = "ellipse", col = "black")
corrplot(corm, method = "circle", col = "black", type = "lower") #type = "lower"只有下三角矩陣
corrplot(corm, method = "square", type = "lower")
corrplot(corm, method = "pie", col = "blue", type = "lower")


###強迫進入法==============================

reg_m <- lm(acad ~ hope + clas + atte + peer + moti + tact, data = reg)
summary(reg_m)

install.packages("psych")
library(psych)
setCor(9, c(3:8), temp, std = TRUE)


###逐步選取法==============================

library(MASS)
stepAIC(reg_m, direction = "forward", k = 2, steps = 1000) #前向選取法（forward method）
summary(stepAIC(reg_m, direction = "forward", k = 2, steps = 1000)) 

stepAIC(reg_m, direction = "backward", k = 2, steps = 1000) #後向選取法（backward method）

stepAIC(reg_m, direction = "both", k = log(120))
stepAIC(reg_m, direction = "both", k = log(nrow(reg)))


###基本套件中的step()函數應用

summary(step(reg_m), method = "forward", k = 2) #forward method
summary(step(reg_m), method = "backward", k = log(nrow(reg))) #backward


##模式比較選取函數

temp <- read.csv("reg_1.csv", header = TRUE, fileEncoding = "BIG-5")
attach(temp)

names(temp)
xvar <- as.matrix(temp[-c(1, 2, 9)])
yvar <- temp$acad

library(leaps)
regm <- summary(regsubsets(xvar, yvar, method = "forward", nbest = 1))
summary(regm)

attach(regm)
round(cbind(which, rsq, rss, cp, bic), 2) #which每個模型中的元素 #rsq為R平方 #rss為殘差均方和 #cp為Mallows' Cp值

regm <- summary(regsubsets(xvar, yvar, method = "backward", nbest = 2)) #改挑兩個最佳模式
attach(regm)
round(cbind(which, rsq, rss, cp, bic), 2)


##虛擬變數的轉換

temp <- read.csv("dummy_1.csv", header = TRUE, fileEncoding = "BIG-5")
attach(temp)

temp$d1_scale <- ifelse(scale == 1, 1, 0)
temp$d2_scale <- ifelse(scale == 2, 1, 0) #參照組為scale = 3

temp$d1_post <- ifelse(post == 1, 1, 0)
temp$d2_post <- ifelse(post == 2, 1, 0)
temp$d3_post <- ifelse(post == 3, 1, 0) #參照組為post = 4

print.data.frame(temp)

attach(temp) #內容更動過了所以要重新attach

reg.model <- lm(press ~ d1_scale + d2_scale + d1_post + d2_post + d3_post)
summary(reg.model)

library(car) #也可以用car套件的recode()增列虛擬變項
temp$d1_sc <- recode(scale, "1 = 1; else = 0")
temp$d3_sc <- recode(scale, "3 = 1; else = 0") #參照組為scale = 2
temp$d1_po <- recode(post, "1 = 1; else = 0")
temp$d3_po <- recode(post, "3 = 1; else = 0")
temp$d4_po <- recode(post, "4 = 1; else = 0")
print.data.frame(temp)


##二次曲線迴歸

temp <- read.csv("collinear.csv", header = TRUE, fileEncoding = "BIG-5")
attach(temp)
names(temp)


###簡單相關分析

cor.test(temp$anxiety, temp$score)
cor.test(temp[, 8], temp[, 9])


###線性迴歸分析

reg.lm <- lm(score ~ anxiety)
summary(reg.lm)


##散佈圖的繪製

par(mfrow = c(1, 1))
plot(temp[, 8], temp[, 9], pch = 16, cex = 2)
grid(nx = 12, ny = 10) #nx、ny界定網格線的數目

plot(temp[, 8], temp[, 9], pch = 15:16, col = 1:2, cex = 2)
grid(nx = 12, ny = 10)
legend("topright", legend = levels(factor(temp$sex)), pch = 15:16, col = 1:2)


##二次曲線迴歸分析

reg.m <- lm(score ~ anxiety + I(anxiety^2))
summary(reg.m)


plot(temp[, 8], temp[, 9], pch = 16, cex = 2)
grid(nx = 12, ny = 10)
reg.m <- lm(score ~ anxiety + I(anxiety^2))
score.p <- predict(reg.m)
lines(spline(anxiety, score.p), lwd = 3, col = 2)
plot(temp[, 8], temp[, 9], type = "n")
lines(spline(anxiety, score.p), lwd = 3, col = 2)


##多元共線性（collinear/multicollinearity）

###共線性範例

m.reg <- lm(acad ~ invo + like + atte + moti)
summary(m.reg) #跡象為：迴歸係數與相關係數正負相反、R sqr高但顯著的變數很少、不合理的估計值參數

library(DAAG)
vif(m.reg) #大於10就有問題，vif為變異數膨脹因素（variance inflation factor）
round(1/vif(m.reg), 3) #vif與tolerance互為倒數

cormatrix <- as.matrix(round(cor(temp[, 3:7]), 2)) #cor()求相關矩陣
cormatrix

cormatrix[lower.tri(cormatrix, diag = TRUE)] #選取相關矩陣左下角的部分
cormatrix[lower.tri(cormatrix, diag = FALSE)] #不要對角線的部分
which(cormatrix[lower.tri(cormatrix, diag = FALSE)] >= .900)


###主成份迴歸分析

pca <- prcomp(temp[, 3:6]) #共線性的解決方式：主成份分析、脊迴歸、排除變數
summary(pca) #前兩個主成份解釋的變異量達98.1%
pca$x #樣本觀察值在四個主成份的分數
round(cor(pca$x), 3) #主成份間的相關矩陣

pca.reg <- lm(acad ~ pca$x[, 1:2]) #只用前兩個主成份做迴歸
summary(pca.reg)


pc <- princomp(temp[, 3:6]) #也可以用princomp()函數
summary(pc)
pc$score #樣本觀察值在各主成份的分數，用score()
round(cor(pc$score), 3)

pc.reg <- lm(acad ~ pc$score[, 1:2])
summary(pc.reg)


#徑路分析

temp <- read.csv("reg_1.csv", header = TRUE, fileEncoding = "BIG-5")
names(temp)

library(psych) #psych套件中的mediate()函數進行徑路分析（path analysis）
mediate(y = 9, x = 7, m = 8, data = temp, std = T)

mediate(y = 9, x = 5, m = 7:8, data = temp, std = TRUE)


library(psych)
setCor(y = 9, x = c(7, 8), data = temp, std = TRUE) #setCor()適用兩個以上自變數之徑路圖

setCor(y = c(7, 9), x = c(3, 4, 5, 8), data = temp, std = TRUE)







#基礎統計分析 R程式在社會科學之應用==============================
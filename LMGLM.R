#R軟體資料分析基礎與應用==================================================

#17廣義線性模型==================================================

##羅吉斯迴歸（Logistic Regression）==================================================

acs <- read.table("http://jaredlander.com/data/acs_ny.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
acs$Income <- with(acs, FamilyIncome >= 150000)
require(ggplot2)
require(useful)
ggplot(acs, aes(x = FamilyIncome)) + geom_density(fill = "grey", color = "grey") + geom_vline(xintercept = 150000) + scale_x_continuous(label = multiple.dollar, limits = c(0, 1000000))

head(acs)

income1 <- glm(Income ~ HouseCosts + NumWorkers + OwnRent + NumBedrooms + FamilyType, 
               data = acs, family = binomial(link = "logit"))

summary(income1)

invlogit <- function(x){
        1/(1 + exp(-x))
}
invlogit(income1$coefficients) #解讀羅吉斯迴歸中的係數要先用羅吉反函數對係數作轉換！


##泊松迴歸模型==================================================

#適合運用在計數資料上

ggplot(acs, aes(x = NumChildren)) + geom_histogram(binwidth = 1) #家庭小孩數量直方圖

children1 <- glm(NumChildren ~ FamilyIncome + FamilyType + OwnRent, data = acs, 
                 family = poisson(link = "log"))
summary(children1)

install.packages("coefplot")
library(coefplot)
coefplot(children1)

#泊松迴歸可能存在過度離散（overdispersion）的問題

z <- (acs$NumChildren - children1$fitted.values) / sqrt(children1$fitted.values)
sum(z ^ 2) / children1$df.residual #計算過度離散率（OD），>= 2表示過度離散，或是< 2但p值為1

pchisq(sum(z ^ 2), children1$df.residual) #過度離散p值

#過度離散問題顯著，用準泊松分佈（quasi poisson）或負二項分佈（negative binomial）重新建模

children2 <- glm(NumChildren ~ FamilyIncome + FamilyType + OwnRent, data = acs, 
                 family = quasipoisson(link = "log"))
multiplot(children1, children2)


##其他廣義線性模型==================================================

#glm還支援伽瑪(Gamma)、反高斯（inverse gaussian）、準二項（quasibinomial）回歸
#對他們使用不同的連結函數（link functions）
#Gamma：inverse、identity、log
#Poisson：log、identity、sqrt
#inverse gaussian：1/mu^2、inverse、identity、log

#多項回歸模型：可以用來對幾種類別進行分類，可以執行好幾個羅吉斯回歸得到相同結果，或使用nnet套件的polr()函數或multinom()函數


##倖存分析==================================================

require(survival)
head(bladder) #stop（事件發生或病人離開研究的時間）、event（在該時間是否有發生事件）

bladder[100:105, ]

survObject <- with(bladder[100:105, ], Surv(stop, event)) #Surv()建立反應變數
survObject
survObject[, 1:2] #首兩列有事件發生，發生時間為12，最後兩列沒事件發生，但可能在該時間後發生，因此該資料的發生時間被設限了（censored）


#倖存分析最常使用的模型為Cox比例風險模型（Proportional Hazard Model）

cox1 <- coxph(Surv(stop, event) ~ rx + number + size + enum, data = bladder)
summary(cox1)

#倖存曲線顯示「在某個時間點有多少比例的受試者存活」

plot(survfit(cox1), xlab = "Days", ylab= "Survival Rate", conf.int = TRUE)


#此處的rx變數是病人接受治療或安慰劑的指標，可傳遞到strata機贓料分成兩群來分析，產生兩條倖存曲線

cox2 <- coxph(Surv(stop, event) ~ strata(rx) + number + size + enum, data = bladder)
summary(cox2)

plot(survfit(cox2), xlab = "Days", ylab = "Survival Rate", onf.int = TRUE, col = 1:2)
legend("bottomleft", legend = c(1, 2), lty = 1, col = 1:2, text.col = 1:2, title = "rx")

cox.zph(cox1) #檢測比例風險模型的假設
cox.zph(cox2)


#Andersen-Gill分析和倖存分析相似，但處理的是區間資料，而且可以處理多個事件
#例如不只可以處理一間急診室是否有人求診，還能計算出急診室求診個數
#同樣用coxph()，要加一個附加變數到Surv，且必須根據用來識別資料的欄位（id）對資料分群

head(bladder2)
ag1 <- coxph(Surv(start, stop, event) ~ rx + number + size + enum + cluster(id), data = bladder2)

ag2 <- coxph(Surv(start, stop, event) ~ strata(rx) + number + size + enum + cluster(id), data = bladder2)

plot(survfit(ag1), conf.int = TRUE)
plot(survfit(ag2), conf.int = TRUE, col = 1:2)
legend("topright", legend = c(1, 2), lty = 1, col = 1:2, text.col = 1:2, title = "rx")
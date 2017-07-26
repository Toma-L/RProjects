#資料挖礦與大數據分析==================================================

#關聯規則==================================================

## support（支持度）：前提項目X與結果項目Y一起出現的機率
## confidence（信賴度）：前提項目X發生的情況下，結果項目Y發生的條件機率
## lift（增益）：信賴度與結果項目Y單獨發生的機率大小比值


## Apriori：廣度 + 水平
## Partition Apriori：廣度 + 垂直，資料庫分區段
## DHP：廣度 + 水平，利用雜湊表
## MSApriori：廣度 + 水平，挖掘低頻率但重要事件的關聯規則
## FP-growth：深度 + 水平，改善Apriori無法處理大量資料的缺點


library(arules)
library(arulesViz)

data("IncomeESL")
class(IncomeESL)
IncomeESL %>% head
class(IncomeESL[, 4])
IncomeESL <- IncomeESL[complete.cases(IncomeESL), ] # 刪除遺漏值
dim(IncomeESL)
Income <- as(IncomeESL, "transactions") # 換成可以進行關聯分析的transactions物件，每個屬性值轉化為單一item
inspect(Income[1:10])
sort(itemFrequency(Income), decreasing = TRUE)
itemFrequencyPlot(Income, support = .2, cex.names = .8)

rules <- apriori(Income, parameter = list(support = .1, confidence = .6)) # 支持度門檻.1，信賴度門檻.6
summary(rules)
plot(rules, measure = c("confidence", "lift"), shading = "support")
plot(rules, method = "grouped") # 圓圈大小：支持度，顏色深淺：增益
rulesOwn <- subset(rules, subset = rhs %in% "householder status=own" & lift > 1) # 想觀察特定族群
inspect(head(sort(rulesOwn, by = "support"), n = 5)) # inspect()


data("IncomeESL")
IncomeESL <- IncomeESL[complete.cases(IncomeESL), ]
IncomeESL[["income"]] <- factor((as.numeric(IncomeESL[["income"]]) > 6) + 1, 
                                levels = 1:2, labels = c("$40-", "$40+")) # 對資料重新編碼

Income <- as(IncomeESL, "transactions")

rules <- apriori(Income, parameter = list(support = .2, confidence = .6))
rulesIncome <- subset(rules, subset = rhs %in% "income=$40+" & lift > 1)

inspect(sort(rulesIncome, by = "confidence"))


## 關聯規則的應用

### 分析顧客行為
### 進行市場區隔與選擇目標顧客
### 改進賣場陳設與實行目標行銷
### 組合搭售商品
### 發掘詐欺行為
### 流失顧客分析

# 並非全部符合篩選指標的關聯規則都可以應用，必須同時經過領域知識的推論與評估

# 不平衡資料可以先將比例較大的資料先刪除
# 避免關聯規則太冗長，可以僅選取變數組合長度較短的


# 利用R語言打通大數據的經脈==================================================

# 06連結分析==================================================

# install.packages("arules")
library(arules)

# Apriori效率較低
# Eclat執行效率有所提升
# FP-Growth高效最佳化演算法

# Apriori參數預設值
## support = .1
## confidence = .8
## maxlen = 10
## minlen = 1
## target = "rules"/"frequent itemsets"（輸出連結規則/頻繁項集）

## appearance：對先決條件X（lhs）和連結結果Y（rhs）實際包含哪些項進行限制
## control：控制函數效能，對項集進行升冪（sort = 1）或降冪（sort = -1）排序，是否向使用者報告處理程序（verbose = TREU/FALSE）


library(arules)
data("Groceries")
summary(Groceries)
inspect(Groceries[1:10])  #inspect()

rules0 <- apriori(Groceries, parameter = list(support = .001, confidence = .5))
rules0 # 5668筆規則
inspect(rules0[1:10])

# 先設定得很低，再加強support或confidence來調整，設定值較高容易遺失有用資訊

rules1 <- apriori(Groceries, parameter = list(support = .005, confidence = .5))
rules1 # 120筆

rules2 <- apriori(Groceries, parameter = list(support = .005, confidence = .60))
rules2 # 22筆

rules3 <- apriori(Groceries, parameter = list(support = .005, confidence = .64))
rules3 # 4筆
inspect(rules3)


rules.sorted_sup <- sort(rules0, by = "support") # 透過support控制
inspect(rules.sorted_sup[1:5])

rules.sorted_con <- sort(rules0, by = "confidence") # 透過confidence控制
inspect(rules.sorted_con[1:5])

rules.sorted_lift <- sort(rules0, by = "lift") # 透過lift控制
inspect(rules.sorted_lift[1:5])

# 用lift來篩選連結規則是最可靠的指標，結論也常常是最有用的


# 只想知道芥末（mustard）的強連結商品，且只要2個商品連結
rules4 <- apriori(Groceries, parameter = list(maxlen = 2, supp = 0.001, conf = .1),
                  appearance = list(rhs = "mustard", default = "lhs"))
inspect(rules4)

itemsets_apr <- apriori(Groceries, parameter = list(supp = .001, target = "frequent itemsets"), 
                        control = list(sort = -1))
itemsets_apr
inspect(itemsets_apr[1:5]) # 觀察銷量最高的商品（supp = 0.001, sort = -1）


# 用eclat()來取得最適合進行bundle銷售的商品（eclat()無法產生關聯規則）
itemsets_ecl <- eclat(Groceries, parameter = list(minlen = 1, maxlen = 3, supp = .001, target = "frequent itemsets"),
                      control = list(sort = -1))
itemsets_ecl
inspect(itemsets_ecl[1:5]) # 觀察前5個頻繁項目集（eclat(target = "frequent itemsets)）
# 頻繁項集只和support設定有關，confidence值不影響


library(arulesViz)
rules5 <- apriori(Groceries, parameter = list(support = .002, confidence = .5))
plot(rules5) # 顏色深淺為lift值高低

plot(rules5, measure = c("support", "lift"), shading = "confidence") # 改成由confidence決定顏色

plot(rules5, interactive = TRUE) # 設定互動參數interactive

plot(rules5, shading = "order", control = list(main = "Two key plot")) # shading = "order"點的顏色深淺代表連結規則中有多少樣商品

plot(rules5, method = "grouped") # lift是顏色深淺，support是尺寸大小

plot(rules5[1:50], method = "matrix", measure = "lift")
plot(rules5[1:50], method = "matrix3D", measure = "lift")
plot(rules5[1:50], method ="paracoord")
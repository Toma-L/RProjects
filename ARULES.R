#資料挖礦與大數據分析==================================================

#關聯規則==================================================

##support（支持度）：前提項目X與結果項目Y一起出現的機率
##confidence（信賴度）：前提項目X發生的晴望下，結果項目Y發生的條件機率
##lift（增益）：信賴度與結果項目Y單獨發生的機率大小比值


##Apriori：廣度 + 水平
##Partition Apriori：廣度 + 垂直，資料庫分區段
##DHP：廣度 + 水平，利用雜湊表
##MSApriori：廣度 + 水平，挖掘低頻率但重要事件的關聯規則
##FP-growth：深度 + 水平，改善Apriori無法處理大量資料的缺點


library(arules)
library(arulesViz)

data("IncomeESL")
IncomeESL <- IncomeESL[complete.cases(IncomeESL), ] #刪除遺漏值
dim(IncomeESL)

Income <- as(IncomeESL, "transactions") #換成可以進行關聯分析的transactions物件，每個屬性值轉化為單一item
sort(itemFrequency(Income), decreasing = TRUE)
itemFrequencyPlot(Income, support = .2, cex.names = .8)

rules <- apriori(Income, parameter = list(support = .1, confidence = .6)) #支持度門檻.1，信賴度門檻.6
summary(rules)
plot(rules, measure = c("confidence", "lift"), shading = "support")
plot(rules, method = "grouped") #圓圈大小：支持度，顏色深淺：增益
rulesOwn <- subset(rules, subset = rhs %in% "householder status=own" & lift > 1) #想觀察特定族群
inspect(head(sort(rulesOwn, by = "support"), n = 5)) #inspect()


data("IncomeESL")
IncomeESL <- IncomeESL[complete.cases(IncomeESL), ]
IncomeESL[["income"]] <- factor((as.numeric(IncomeESL[["income"]]) > 6) + 1, 
                                levels = 1:2, labels = c("$40-", "$40+")) #對資料重新編碼

Income <- as(IncomeESL, "transactions")

rules <- apriori(Income, parameter = list(support = .2, confidence = .6))
rulesIncome <- subset(rules, subset = rhs %in% "income=$40+" & lift > 1)

inspect(sort(rulesIncome, by = "confidence"))


##關聯規則的應用

###分析顧客行為
###進行市場區隔與選擇目標顧客
###改進賣場陳設與實行目標行銷
###組合搭售商品
###發掘詐欺行為
###流失顧客分析

#並非全部符合篩選指標的關聯規則都可以應用，必須同時經過領域知識的推論與評估

#不平衡資料可以先將比例較大的資料先刪除
#避免關聯規則太冗長，可以僅選取變數組合長度較短的
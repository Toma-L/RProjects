# Data Visualization in R

# Thomas JH Lin / Akihiro Hayashi / 임소홍


# R Graphics Cookbook ==================================================

library(ggplot2)

# 1. R 基礎 =====

# 2. 快速探索數據 =====

## 2.1 散佈圖 geom_point() =====

qplot(mtcars$wt, mtcars$mpg)
qplot(wt, mpg, data = mtcars)
ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()
str(mtcars)



## 2.2 折線圖 geom_line() =====

plot(pressure$temperature, pressure$pressure, type = "l") # 先plot()再加東西
points(pressure$temperature, pressure$pressure) # 記得是point"s"
lines(pressure$temperature, pressure$pressure / 2, col = "red")
points(pressure$temperature, pressure$pressure / 2, col = "red")

qplot(pressure$temperature, pressure$pressure, geom = "line") # geom = 
qplot(temperature, pressure, data = pressure, geom = "line") # 方法2
ggplot(pressure, aes(x = temperature, y = pressure)) + geom_line() #方法3

qplot(temperature, pressure, data = pressure, geom = c("line", "point")) # 一次寫完



## 2.3 長條圖 geom_bar() =====

barplot(BOD$demand, names.arg = BOD$Time)
table(mtcars$cyl)
barplot(table(mtcars$cyl))

# qplot(BOD$Time, BOD$demand, geom = "bar", stat = "identity") # Error
# qplot(factor(BOD$Time), BOD$demand, geom = "bar", stat = "identity") # Error

qplot(mtcars$cyl)
class(mtcars$cyl)
qplot(factor(mtcars$cyl))
# qplot(Time, demand, data = BOD, geom = "bar", stat = "identity") # Error
ggplot(BOD, aes(x = Time, y = demand)) + geom_bar(stat = "identity")
ggplot(BOD, aes(x = factor(Time), y = demand)) + geom_bar(stat = "identity")

qplot(factor(cyl), data = mtcars)
ggplot(mtcars, aes(x = factor(cyl))) + geom_bar()



## 2.4 直方圖 geom_histogram(binwidth) =====

hist(mtcars$mpg)
hist(mtcars$mpg, breaks = 10)

qplot(mtcars$mpg)
qplot(mpg, data = mtcars, binwidth = 4) # binwidth
ggplot(mtcars, aes(x = mpg)) + geom_histogram(binwidth = 4)



## 2.5 箱型圖 geom_boxplot(); interaction() =====

plot(ToothGrowth$supp, ToothGrowth$len) # x軸變數為factor，所以自動畫出boxplot
class(ToothGrowth$supp) # factor
boxplot(len ~ supp, data = ToothGrowth)
boxplot(len ~ supp + dose, data = ToothGrowth)

qplot(ToothGrowth$supp, ToothGrowth$len, geom = "boxplot") # 方法2
ggplot(ToothGrowth, aes(x = supp, y = len)) + geom_boxplot() # 方法3

# 使用 interaction() 創造更多分組
unique(interaction(ToothGrowth$supp, ToothGrowth$dose))
qplot(interaction(ToothGrowth$supp, ToothGrowth$dose), ToothGrowth$len, geom = "boxplot")
qplot(interaction(supp, dose), len, data = ToothGrowth, geom = "boxplot") # another way
ggplot(ToothGrowth, aes(x = interaction(supp, dose), y = len)) + geom_boxplot() # the same



## 2.6 函數圖形 curve(); stat_function(fun, geom) =====

curve(x ^ 3 - 5 * x, from = -4, to = 4)

myfun <- function(xvar) {
        1 / (1 + exp(-xvar + 10))
}
curve(myfun(x), from = 0, to = 20)
curve(1 - myfun(x), add = TRUE, col = "red")


# qplot(c(0, 20), fun = myfun, stat = "function", geom = "line") # Error
ggplot(data.frame(x = c(0, 20)), aes(x = x)) + 
        stat_function(fun = myfun, geom = "line")



# 3. 長條圖 =====

## 3.1 簡單長條圖 geom_bar() =====

library(gcookbook)
ggplot(pg_mean, aes(x = group, y = weight)) + geom_bar(stat = "identity")

BOD
str(BOD)
ggplot(BOD, aes(x = Time, y = demand)) + geom_bar(stat = "identity")
ggplot(BOD, aes(x = factor(Time), y = demand)) + 
        geom_bar(stat = "identity") # 轉為離散型變量

ggplot(pg_mean, aes(x = group, y = weight)) + 
        geom_bar(stat = "identity", fill = "lightblue", colour = "black")



## 3.2 簇狀條形圖 position = "dodge" =====

cabbage_exp
class(cabbage_exp$Cultivar); class(cabbage_exp$Date)
ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) +  # fill要用factor變數
        geom_bar(stat = "identity", position = "dodge") # 水平排列

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(position = "dodge", stat = "identity", colour = "black") +
        scale_fill_brewer(palette = "Pastel1") # scale_fill_brewer(palette = "")

ce <- cabbage_exp[1:5, ]
ce

ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(position = "dodge", stat = "identity", colour = "black") +
        scale_fill_brewer(palette = "Pastel1") # 有缺失的組合，鄰近的長條會自動擴充



## 3.3 頻度長條圖 geom_bar() =====

ggplot(diamonds, aes(x = cut)) + geom_bar() # 離散型會得到bar chart
ggplot(diamonds, aes(x = carat)) + geom_bar() # 連續型會得到histogram



## 3.4 長條圖著色 fill, colour =====

upc <- subset(uspopchange, rank(Change) > 40)
upc

ggplot(upc, aes(x = Abb, y = Change, fill = Region)) + 
        geom_bar(stat = "identity")

ggplot(upc, aes(x = reorder(Abb, Change), y = Change, fill = Region)) +  # reorder()按Change由低到高排序Abb
        geom_bar(stat = "identity", colour = "black") +
        scale_fill_manual(values = c("#669933", "#FFCC66")) +  # scale_fill_manual()
        xlab("State")



## 3.5 正負長條圖 position = "identity" =====

csub <- subset(climate, Source == "Berkeley" & Year >= 1900)
csub$pos <- csub$Anomaly10y >= 0
csub

ggplot(csub, aes(x = Year, y = Anomaly10y, fill = pos)) + 
        geom_bar(stat = "identity", position = "identity")
# position = "identity" 可以避免系統畫負值而產生的警告訊息
# 但是暖色通常用於正值才對

ggplot(csub, aes(x = Year, y = Anomaly10y, fill = pos)) + 
        geom_bar(stat = "identity", position = "identity", colour = "black", size = .25) + 
        scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"), guide = FALSE)
# size 是調整邊線寬度用的，不是長條寬度！
# guide = FALSE 可以刪除圖例



## 3.6 調整長條寬度和間距 width; position_dodge() =====

ggplot(pg_mean, aes(x = group, y = weight)) + geom_bar(stat = "identity")

ggplot(pg_mean, aes(x = group, y = weight)) + 
        geom_bar(stat = "identity", width = .5)

ggplot(pg_mean, aes(x = group, y = weight)) + 
        geom_bar(stat = "identity", width = 1)

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", width = .5, position = "dodge")

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", width = .5, position = position_dodge(.7))
# position = position_dodge()

# 四式等價
geom_bar(position = "dodge")
geom_bar(width = .9, position = position_dodge())
geom_bar(position = position_dodge(.9))
geom_bar(width = .9, position = position_dodge(.9))



## 3.7 堆積長條圖 stat = "identity" =====

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity") # 就不加position = "dodge"就是了，只留stat = "identity"
cabbage_exp

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity") + 
        guides(fill = guide_legend(reverse = TRUE))
# guides() 調整圖例，調整圖例的填充色順序用 fill = guide_legend(reverse = TRUE)

library(plyr) # 為了使用 order = desc()
ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar, order = desc(Cultivar))) + 
        geom_bar(stat = "identity")
# 注意！！！調整堆疊順序用 order = desc() ---> NO
# NOOOOOOOOOOO! order = desc()不是正確的做法，應該要用order把樣本重新排序
# 他的堆疊方式是在資料集中先出現的先疊

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", colour = "black") + 
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_fill_brewer(palette = "Pastel1")



## 3.8 百分比堆積長條圖 ddply(, transform) =====

library(plyr)
# 注意！！！
ce <- ddply(cabbage_exp, "Date", transform,  # 以Date為切割變數進行transform()
            percent_weight = Weight / sum(Weight) * 100)
ce
ggplot(ce, aes(x = Date, y = percent_weight, fill = Cultivar)) + 
        geom_bar(stat = "identity")

cabbage_exp
ce <- ddply(cabbage_exp, "Date", transform,
            percent_weight = Weight / sum(Weight) * 100)

ggplot(ce, aes(x = Date, y = percent_weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", colour = "black") +
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_fill_brewer(palette = "Pastel1")



## 3.9 加上資料標籤 geom_text(aes(label = ), vjust, hjust) =====

ggplot(cabbage_exp, aes(x = interaction(Date, Cultivar), y = Weight)) +
        geom_bar(stat = "identity") + 
        geom_text(aes(label = Weight), vjust = 1.5, colour = "white") # 頂端下方

ggplot(cabbage_exp, aes(x = interaction(Date, Cultivar), y = Weight)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = Weight), vjust = -.2) # 頂端上方
# vjust 調整標籤在長條圖頂端的上或下

ggplot(cabbage_exp, aes(x = interaction(Date, Cultivar), y = Weight)) + 
        geom_bar(stat = "identity") +
        geom_text(aes(label = Weight), vjust = -.2) + 
        ylim(0, max(cabbage_exp$Weight) * 1.05) # 提高y軸上限

ggplot(cabbage_exp, aes(x = interaction(Date, Cultivar), y = Weight)) + 
        geom_bar(stat = "identity") + 
        geom_text(aes(y = Weight + .1, label = Weight)) # 注意！！！設定標籤的y軸位置

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label = Weight), vjust = 1.5, colour = "white", 
                  position = position_dodge(.9), size = 3) # 分類間距的default為0.9


# 直向堆積長條圖要加入標籤之前，要先對各組資料求出累積和！
library(plyr)
ce <- arrange(cabbage_exp, Date, Cultivar)
# 注意！！！stacked barchart如何加入資料值，可用cumsum()
ce <- ddply(ce, "Date", transform, label_y = cumsum(Weight))
ce

ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity") +  # stacked column
        geom_text(aes(y = label_y, label = Weight), vjust = 1.5, colour = "white")

ce <- arrange(cabbage_exp, Date, Cultivar)
ce <- ddply(ce, "Date", transform, label_y = cumsum(Weight) - .5 * Weight) # 把標籤放在長條圖正中間
ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) +
        geom_bar(stat = "identity") + 
        geom_text(aes(y = label_y, label = Weight), colour = "White")

ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", colour = "black") +
        geom_text(aes(y = label_y, label = paste(format(Weight, nsmall = 2), "kg")), size = 4) +
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_fill_brewer(palette = "Pastel1")



## 3.10 Cleveland點圖 geom_point() =====

library(gcookbook)
tophit <- tophitters2001[1:25, ]
ggplot(tophit, aes(x = avg, y = name)) + geom_point()

tophit[, c("name", "lg", "avg")]
# 注意！！！用reorder()讓name照avg大小排序（此函數一次只能按照一個變量來排序）
ggplot(tophit, aes(x = avg, y = reorder(name, avg))) + 
        geom_point(size = 3) + 
        theme_bw() + 
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_line(colour = "grey60", 
                                                linetype = "dashed"))

ggplot(tophit, aes(x = reorder(name, avg), y = avg)) + 
        geom_point(size = 3) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_line(colour = "grey60", 
                                                linetype = "dashed"))

# 注意！！！手動實現多重變量排序
nameorder <- tophit$name[order(tophit$lg, tophit$avg)]
tophit$name <- factor(tophit$name, levels = nameorder) # 步驟2

# 注意！！！geom_segment()的用法
ggplot(tophit, aes(x = avg, y = name)) + 
        geom_segment(aes(yend = name), xend = 0, colour = "grey50") + 
        geom_point(size = 3, aes(colour = lg)) + 
        scale_colour_brewer(palette = "Set1", limits = c("NL", "AL")) + 
        theme_bw() + 
        theme(panel.grid.major.y = element_blank(),
              legend.position = c(1, .55), # 圖例靠右置中上
              legend.justification = c(1, .5))

ggplot(tophit, aes(x = avg, y = name)) + 
        geom_segment(aes(yend = name), xend = 0, colour = "grey50") + 
        geom_point(size = 3, aes(colour = lg)) + 
        scale_colour_brewer(palette = "Set1", limits = c("NL", "AL"), guide = FALSE) + 
        theme_bw() + 
        theme(panel.grid.major.y = element_blank()) + 
        facet_grid(lg ~ ., scales = "free_y", space = "free_y") # 注意！！！free_y


# 4. 折線圖 =====

## 4.1 簡單折線圖 geom_line() =====

ggplot(BOD, aes(x = Time, y = demand)) + geom_line()
BOD

BOD1 <- BOD
BOD1$Time <- factor(BOD1$Time) # 連續型變量轉factor
# 注意！！！
ggplot(BOD1, aes(x = Time, y = demand, group = 1)) + geom_line() # 對於factor變量，要加入 group = 1 確認為同一組資料

# 拓寬y軸範圍
ggplot(BOD, aes(x = Time, y = demand)) + geom_line() + ylim(0, max(BOD$demand))
ggplot(BOD, aes(x = Time, y = demand)) + geom_line() + expand_limits(y = 0)



## 4.2 加上資料標記 geom_line() + geom_point() =====

ggplot(BOD, aes(x = Time, y = demand)) + geom_line() + geom_point()

library(gcookbook)
ggplot(worldpop, aes(x = Year, y = Population)) + 
        geom_line() + geom_point()
# 注意！！！scale_y_log10()
ggplot(worldpop, aes(x = Year, y = Population)) + 
        geom_line() + geom_point() + scale_y_log10()



## 4.3 多重折線圖 colour, linetype, shape =====

library(plyr)
tg <- ddply(ToothGrowth, c("supp", "dose"), summarise, length = mean(len))
ggplot(tg, aes(x = dose, y = length, colour = supp)) + geom_line()
ggplot(tg, aes(x = dose, y = length, linetype = supp)) + geom_line() # 改變線的形狀

# 注意！！！一定一定要加 group = supp 讓電腦知道資料是一組的
ggplot(tg, aes(x = factor(dose), y = length, colour = supp, group = supp)) + 
        geom_line()

# 錯誤示範
ggplot(tg, aes(x = dose, y = length)) + geom_line() # 沒有正確分組造成一個x對應不只一個點

ggplot(tg, aes(x = dose, y = length, shape = supp)) + geom_line() + 
        geom_point(size = 4)

ggplot(tg, aes(x = dose, y = length, fill = supp)) + geom_line() + 
        geom_point(size = 4, shape = 21)

# 注意！！！標記可能互相重疊，可以適當地左右移動 position = position_dodge()
ggplot(tg, aes(x = dose, y = length, shape = supp)) + 
        geom_line(position = position_dodge(.2)) + 
        geom_point(position = position_dodge(.2), size = 4)


## 4.4 修改線條樣式 geom(linetype, size, colour) =====

ggplot(BOD, aes(x = Time, y = demand)) + 
        geom_line(linetype = "dashed", size = 1, colour = "blue")

library(plyr)
tg <- ddply(ToothGrowth, c("supp", "dose"), summarise, length = mean(len))
ggplot(tg, aes(x = dose, y = length, colour = supp)) + 
        geom_line() + 
        scale_colour_brewer(palette = "Set1")

ggplot(tg, aes(x = dose, y = length, group = supp)) + 
        geom_line(linetype = "dashed") + 
        geom_point(shape = 22, size = 3, fill = "white")

ggplot(tg, aes(x = dose, y = length, group = supp)) + 
        geom_line(colour = "darkgreen", size = 1.5)

ggplot(tg, aes(x = dose, y = length, colour = supp)) +
        geom_line(linetype = "dashed") + 
        geom_point(shape = 22, size = 3, fill = "white")



## 4.5 修改資料標記 geom_point(size, shape, colour, fill) =====

ggplot(BOD, aes(x = Time, y = demand)) + 
        geom_line() + 
        geom_point(size = 4, shape = 22, colour = "darkred", fill = "pink")

ggplot(BOD, aes(x = Time, y = demand)) + 
        geom_line() + 
        geom_point(size = 4, shape = 21, fill = "white")

library(plyr)
tg <- ddply(ToothGrowth, c("supp", "dose"), summarise, length = mean(len))
pd <- position_dodge(.2)
ggplot(tg, aes(x = dose, y = length, fill = supp)) + 
        geom_line(position = pd) + 
        geom_point(shape = 21, size = 3, position = pd) +
        scale_fill_manual(values = c("black", "white"))



## 4.6 面積圖 geom_area(colour, fill, alpha) =====

sunspotyear <- data.frame(
        Year = as.numeric(time(sunspot.year)),
        Sunspots = as.numeric(sunspot.year))

ggplot(sunspotyear, aes(x = Year, y = Sunspots)) + geom_area()

ggplot(sunspotyear, aes(x = Year, y = Sunspots)) + 
        geom_area(colour = "black", fill = "blue", alpha = .2) # 調整填充色

ggplot(sunspotyear, aes(x = Year, y = Sunspots)) + 
        geom_area(fill = "blue", alpha = .2) +  # 注意！！！不設定colour就沒有底部框線
        geom_line()



## 4.7 堆積面積圖 fill =====

ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup)) + 
        geom_area()
str(uspopage)
# 注意！！！default 的堆積順序與圖標有時是相反的，可以用 breaks 參數調整
ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup)) + 
        geom_area(colour = "black", size = .2, alpha = .4) +
        scale_fill_brewer(palette = "Blues", breaks = rev(levels(uspopage$AgeGroup)))

library(plyr)
# 可以用 order = desc() 反轉堆積順序 ---> No! 只有圖例中的順序會變！！！
ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup, order = desc(AgeGroup))) +
        geom_area(colour = "black", size = .2, alpha = .4) +  # size是框線粗細
        scale_fill_brewer(palette = "Blues")

# 注意！！！geom_line(position = "stack")
ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup, order = desc(AgeGroup))) +
        geom_area(colour = NA, alpha = .4) +  # 不設定colour就不會有底線
        scale_fill_brewer(palette = "Blues") +
        geom_line(position = "stack", size = .2) # 堆積面積圖加線



## 4.8 百分比堆積面積圖 ddply() =====

library(gcookbook)
library(plyr)
uspopage_prop <- ddply(uspopage, "Year", transform, Percent = Thousands / sum(Thousands) * 100)
# ddply 將 uspopage 資料集，按照 Year 拆成多個獨立的 data frames

ggplot(uspopage_prop, aes(x = Year, y = Percent, fill = AgeGroup)) + 
        geom_area(colour = "black", size = .2, alpha = .4) + 
        scale_fill_brewer(palette = "Blues", breaks = rev(levels(uspopage$AgeGroup)))

# scale_fill_brewer(breaks = rev(levels())); ggplot(aes(order = desc()))
# 這兩個方法都只是改變圖例中的順序


## 4.9 增加信賴區間 geom_ribbon(aes(ymin, ymax)) =====

clim <- subset(climate, Source == "Berkeley", select = c("Year", "Anomaly10y", "Unc10y"))
clim
# 注意！！！geom_ribbon() 要先畫，才不會讓 geom_line() 糊掉！
ggplot(clim, aes(x = Year, y = Anomaly10y)) + 
        geom_ribbon(aes(ymin = Anomaly10y - Unc10y, ymax = Anomaly10y + Unc10y), alpha = .2) +
        geom_line()

ggplot(clim, aes(x = Year, y = Anomaly10y)) + 
        geom_line(aes(y = Anomaly10y - Unc10y), colour = "grey50", linetype = "dotted") +
        geom_line(aes(y = Anomaly10y + Unc10y), colour = "grey50", linetype = "dotted") +
        geom_line()


# 5. 散佈圖 =====

## 5.1 散佈圖 geom_point(shape, size) =====

library(gcookbook)
heightweight[, c("ageYear", "heightIn")]
ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point()
ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point(shape = 21)
ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point(size = 1.5)


## 5.2 修改點的樣式 colour, shape =====

library(gcookbook)
heightweight[, c("sex", "ageYear", "heightIn")]

ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + geom_point()
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex)) + geom_point()

ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, colour = sex)) + 
        geom_point()

ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, colour = sex)) + 
        geom_point() + 
        scale_shape_manual(values = c(1, 2)) +  # 人工選擇點的形狀
        scale_colour_brewer(palette = "Set1") # 人工選擇點的顏色



## 5.3 使用非內建的點形 shape =====

ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point(shape = 3)
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex)) + 
        geom_point(size = 3) + 
        scale_shape_manual(values = c(1, 4))

hw <- heightweight
# 注意！！！guide_legend(override.aes = list(shape = ))
hw$weightGroup <- cut(hw$weightLb, breaks = c(-Inf, 100, Inf), labels = c("< 100", ">= 100"))
ggplot(hw, aes(x = ageYear, y = heightIn, shape = sex, fill = weightGroup)) + 
        geom_point(size = 2.5) +
        scale_shape_manual(values = c(21, 24)) +
        scale_fill_manual(values = c(NA, "black"), 
                          guide = guide_legend(override.aes = list(shape = 21)))


## 5.4 將連續型變數映射到點的顏色或大小 colour, size, scale_size_area() =====

heightweight[, c("sex", "ageYear", "heightIn", "weightLb")]

ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = weightLb)) + geom_point()
ggplot(heightweight, aes(x = ageYear, y = heightIn, size = weightLb)) + geom_point()

# 人天生對顏色和大小的變化不太敏銳，因此當一個變數不需要太精密的解釋時才適合用

ggplot(heightweight, aes(x = ageYear, y = heightIn, fill = weightLb)) + 
        geom_point(shape = 21, size = 2.5) + 
        scale_fill_gradient(low = "black", high = "white") # 黑白漸層

# 注意！！！黑白漸層搭配離散圖標
ggplot(heightweight, aes(x = ageYear, y = heightIn, fill = weightLb)) + 
        geom_point(shape = 21, size = 2.5) +
        scale_fill_gradient(low = "black", high = "white", breaks = seq(70, 170, by = 20), 
                            guide = guide_legend())

# 注意！！！scale_size_area()在ggplot()就要設定size！！！
ggplot(heightweight, aes(x = ageYear, y = heightIn, size = weightLb, colour = sex)) + 
        geom_point(alpha = .5) +
        scale_size_area() +  # 使面積與資料值成正比
        scale_colour_brewer(palette = "Set1")

# 不同形狀的點很難比較面積大小，因此不要同時操作


## 5.5 處理資料點重疊問題 stat_bin2d(), stat_binhex(), position = "jitter"=====

# 圖形重疊（overplotting）的解決方案：
## 半透明的點
## 矩形資料分箱
## 六邊形資料分箱
## 箱型圖

sp <- ggplot(diamonds, aes(x = carat, y = price))
sp + geom_point()

sp + geom_point(alpha = .1) # 90%透明度
sp + geom_point(alpha = .01) # 99%透明度
# 注意！！！
sp + stat_bin2d() # 分別在x軸y軸分割30組，共900個箱子
sp + stat_bin2d(bins = 50) + 
        scale_fill_gradient(low = "lightblue", high = "red", limits = c(0, 6000)) # 2500箱

# install.packages("hexbin")
library(hexbin) # 使用六邊形箱子
sp + stat_binhex() + 
        scale_fill_gradient(low = "lightblue", high = "red", limits = c(0, 8000))
sp + stat_binhex() + 
        scale_fill_gradient(low = "lightblue", high = "red", breaks = c(0, 250, 500, 1000, 2000, 4000, 6000), limits = c(0, 6000))
# 範圍外會變成灰色箱子

# 當其中一軸或兩軸為離散型變數時，也會出現overplotting，可用 position_jitter() 增加隨機擾動
sp1 <- ggplot(ChickWeight, aes(x = Time, y = weight))
sp1 + geom_point()
sp1 + geom_point(position = "jitter") # 跟 position_jitter() 意思一樣

sp1 + geom_point(position = position_jitter(width = .5, height = 0)) # width 和 height 調整擾動值的精度

sp1 + geom_boxplot(aes(group = Time))


## 5.6 模型擬合線 stat_smooth(method) =====

library(gcookbook)
sp <- ggplot(heightweight, aes(x = ageYear, y = heightIn))
sp + geom_point() + stat_smooth(method = lm) # default 為95%信賴區間

sp + geom_point() + stat_smooth(method = lm, level = .95) # 調整 level 可調整信賴區間

sp + geom_point(colour = "grey60") + 
        stat_smooth(method = lm, se = FALSE, colour = "black") # se = FALSE就不會畫出信賴區間

# stat_smooth 的 default 是 loess曲線
sp + geom_point(colour = "grey60") + stat_smooth()
sp + geom_point(colour = "grey60") + stat_smooth(method = loess) # the same

library(MASS)
b <- biopsy
b$classn[b$class == "benign"] <- 0
b$classn[b$class == "malignant"] <- 1
b
# 注意！！！
ggplot(b, aes(x = V1, y = classn)) + 
        geom_point(position = position_jitter(width = .3, height = .06), alpha = .4, shape = 21, size = 1.5) + 
        stat_smooth(method = glm, method.args = list(family = "binomial"))

sps <- ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + 
        geom_point() + 
        scale_colour_brewer(palette = "Set1")

sps + stat_smooth() # stat_smooth() 的範圍限定在預測資料對應的範圍內 # 預設為loess

sps + stat_smooth(method = lm, se = FALSE, fullrange = TRUE) # 可外推的模型要加入參數 fullrange = TRUE
sps + stat_smooth(method = lm, fullrange = TRUE) # 可外推的模型要加入參數 fullrange = TRUE



## 5.7 既有模型散佈圖加入擬合線 =====

model <- lm(heightIn ~ ageYear + I(ageYear ^ 2), heightweight)
model

xmin <- min(heightweight$ageYear)
xmax <- max(heightweight$ageYear)
predicted <- data.frame(ageYear = seq(xmin, xmax, length.out = 100))

predicted$heightIn <- predict(model, predicted)
predicted

sp <- ggplot(heightweight, aes(x = ageYear, y = heightIn)) + 
        geom_point(colour = "grey40")

sp + geom_line(data = predicted, size = 1)


predictvals <- function(model, xvar, yvar, xrange = NULL, samples = 100, ...) {
        if (is.null(xrange)) { # 如果xrange沒有輸入，則從模型對象中自動提取x軸範圍
                if (any(class(model) %in% c("lm", "glm"))) # 提取方法視模型而定
                        xrange <- range(model$model[[xvar]])
                else if (any(class(model) %in% "loess"))
                        xrange <- range(model$x)
        }
        newdata <- data.frame(x = seq(xrange[1], xrange[2], length.out = samples))
        names(newdata) <- xvar
        newdata[[yvar]] <- predict(model, newdata = newdata, ...)
        newdata
}

modlinear <- lm(heightIn ~ ageYear, heightweight)
modloess <- loess(heightIn ~ ageYear, heightweight)

lm_predicted <- predictvals(modlinear, "ageYear", "heightIn")
loess_predicted <- predictvals(modloess, "ageYear", "heightIn")

sp + geom_line(data = lm_predicted, colour = "red", size = .8) +
        geom_line(data = loess_predicted, colour = "blue", size = .8)


library(MASS)
b <- biopsy
b$classn[b$class == "benign"] <- 0
b$classn[b$class == "malignant"] <- 1

fitlogistic <- glm(classn ~ V1, b, family = binomial)
glm_predicted <- predictvals(fitlogistic, "V1", "classn", type = "response")
ggplot(b, aes(x = V1, y = classn)) + 
        geom_point(position = position_jitter(width = .3, height = .08), alpha = .4, shape = 21, size = 1.5) + 
        geom_line(data = glm_predicted, colour = "#1177FF", size = 1)



## 5.8 多模型擬合線 =====

make_model <- function(data) {
        lm(heightIn ~ ageYear, data)
}

library(gcookbook)
library(plyr)
models <- dlply(heightweight, "sex", .fun = make_model)
models

predvals <- ldply(models, .fun = predictvals, xvar = "ageYear", yvar = "heightIn")
predvals

ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + geom_point() + 
        geom_line(data = predvals)

predvals <- ldply(models, .fun = predictvals, xvar = "ageYear", yvar = "heightIn", xrange = range(heightweight$ageYear))
# 把兩組的x軸範圍調整成相同
ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + geom_point() + geom_line(data = predvals)



## 5.9 散佈圖加入模型係數 annotate("text", label = "", x, y) =====

library(gcookbook)
model <- lm(heightIn ~ ageYear, heightweight)
summary(model)

pred <- predictvals(model, "ageYear", "heightIn")
sp <- ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point() + geom_line(data = pred)
sp + annotate("text", label = "r^2 = 0.42", x = 16.5, y = 52) # 加入係數標籤
sp + annotate("text", label = "r^2 == 0.42", parse = TRUE, x = 16.5, y = 52) # parse = TRUE，可用R的方式表示數學符號

eqn <- as.character(as.expression(
        substitute(italic(y) == a + b * italic(x) * "," ~~ italic(r)^2 ~ "=" ~ r2, 
                   list(a = format(coef(model)[1], digits = 3),
                        b = format(coef(model)[2], digits = 3),
                        r2 = format(summary(model)$r.squared, digits = 2)
                        ))))
eqn
parse(text = eqn)

sp + annotate("text", label = eqn, parse = TRUE, x = Inf, y = -Inf, hjust = 1.1, vjust = -.5)
# 注意！！！x = Inf, y = -Inf 使公式置於右下角


# 5.10 散佈圖加入邊際地毯 geom_rug() =====

ggplot(faithful, aes(x = eruptions, y = waiting)) + 
        geom_point() + 
        geom_rug()
# marginal rugs 本質上是一個一維的散佈圖
ggplot(faithful, aes(x = eruptions, y = waiting)) + 
        geom_point() + 
        geom_rug(position = "jitter", size = .2)
# position = "jitter", size = .2 調整線寬及減輕重疊程度


## 5.11 散佈圖加標籤 geom_text(aes(label = , x, y)) =====

library(gcookbook)
subset(countries, Year == 2009 & healthexp > 2000)
sp <- ggplot(subset(countries, Year == 2009 & healthexp > 2000), 
             aes(x = healthexp, y = infmortality)) + geom_point()
sp + annotate("text", x = 4350, y = 5.4, label = "Canada") +
        annotate("text", x = 7400, y = 6.8, label = "USA")

sp + geom_text(aes(label = Name), size = 4) # geom_text() 可用factor或char類型的向量製作標籤

sp + geom_text(aes(label = Name), size = 4, vjust = 0)
sp + geom_text(aes(y = infmortality + .1, label = Name), size = 4, vjust = 0)

# 注意！！！左對齊 hjust = 0; 右對齊 hjust = 1
# 用這種方法，較長的標籤會有較大的移動，此時最好用 x 增減一個值來調整

sp + geom_text(aes(label = Name), size = 4, hjust = 0)
sp + geom_text(aes(x = healthexp + 100, label = Name), size = 4, hjust = 0)


# 注意！！！如果不想要全部加上標籤，可以複製一個新的標籤
cdat <- subset(countries, Year == 2009 & healthexp > 2000)
cdat$Name1 <- cdat$Name

idx <- cdat$Name1 %in% c("Canada", "Ireland", "United Kingdom", "United States", "New Zealand", "Iceland",
                         "Japan", "Luxembourg", "Netherlands", "Switzerland")
idx
cdat$Name1[!idx] <- NA
cdat

ggplot(cdat, aes(x = healthexp, y = infmortality)) + geom_point() + 
        geom_text(aes(x = healthexp + 100, label = Name1), size = 4, hjust = 0) + 
        xlim(2000, 10000)



## 5.12 氣泡圖 scale_size_area() =====

library(gcookbook)
cdat <- subset(countries, Year == 2009 & Name %in% c("Canada", "Ireland", "United Kingdom", 
                                                     "United States", "New Zealand", "Iceland","Japan", 
                                                     "Luxembourg", "Netherlands", "Switzerland"))
cdat
p <- ggplot(cdat, aes(x = healthexp, y = infmortality, size = GDP)) + 
        geom_point(shape = 21, colour = "black", fill = "cornsilk")
p
p + scale_size_area(max_size = 15) # 以GDP決定面積

# 當x軸y軸都是類別變數，氣泡圖可以用來表示變量值
hec <- HairEyeColor[, , "Male"] + HairEyeColor[, , "Female"] # HairEyeColor是個list
library(reshape2)
hec <- melt(hec, value.name = "count")
ggplot(hec, aes(x = Eye, y = Hair)) + 
        geom_point(aes(size = count), shape = 21, colour = "black", fill = "cornsilk") + 
        scale_size_area(max_size = 20, guide = FALSE) + 
        geom_text(aes(y = as.numeric(Hair) - sqrt(count) / 22, label = count), vjust = 1, colour = "grey60", size = 4)
# 注意！！！此處的y座標是計算得出



## 5.13 散佈圖矩陣 pairs() =====

library(gcookbook)
c2009 <- subset(countries, Year == 2009, select = c(Name, GDP, laborrate, healthexp, infmortality))
c2009
pairs(c2009[, 2:5])

# 自定義的面板函數
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y, use = "complete.obs"))
        txt <- format(c(r, 0.123456789), digits = digits)[1]
        txt <- paste(prefix, txt, sep = "")
        if(missing(cex.cor)) cex.cor <- 0.8 / strwidth(txt)
        text(.5, .5, txt, cex = cex.cor * (1 + r) / 2)
}

panel.hist <- function(x, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5))
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$counts
        y <- y / max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col = "white", ...)
}

pairs(c2009[, 2:5], 
      upper.panel = panel.cor,
      diag.panel = panel.hist,
      lower.panel = panel.smooth)

panel.lm <- function(x, y, col = par("col"), bg = NA, pch = par("pch"),
                     cex = 1, col.smooth = "black", ...) {
        points(x, y, pch = pch, col = col, bg = bg, cex = cex)
        abline(stats::lm(y ~ x), col = col.smooth, ...)
}

pairs(c2009[, 2:5], pch = ".", # 用小一點的點，方便辨認數據，或用 cex 參數
      upper.panel = panel.cor,
      diag.panel = panel.hist,
      lower.panel = panel.lm)



# 6. 敘述統計 ============================================================

## 6.1 簡單直方圖 geom_histogram(binwidth) =====

ggplot(faithful, aes(x = waiting)) + geom_histogram()

w <- faithful$waiting
ggplot(NULL, aes(x = w)) + geom_histogram()

ggplot(faithful, aes(x = waiting)) + 
        geom_histogram(binwidth = 5, 
                       fill = "white",  # 填滿
                       colour = "black") # 邊框顏色
# 注意！！！
binsize <- diff(range(faithful$waiting)) / 15 # 計算binwidth
ggplot(faithful, aes(x = waiting)) + 
        geom_histogram(binwidth = binsize, fill = "white", colour = "black")

h <- ggplot(faithful, aes(x = waiting)) # 儲存成變量以便重複利用

h + geom_histogram(binwidth = 8, fill = "white", colour = "black", origin = 31) # 設定組邊界
h + geom_histogram(binwidth = 8, fill = "white", colour = "black", origin = 35) # geom_bar(stat = "bin") 可得相同結果



## 6.2 多組資料直方圖 facet_grid() =====

library(MASS) # 取用數據
ggplot(birthwt, aes(x = bwt)) + 
        geom_histogram(fill = "white", colour = "black") + 
        facet_grid(smoke ~ . )

birthwt1 <- birthwt
birthwt1$smoke <- factor(birthwt1$smoke)
levels(birthwt1$smoke) # 標籤是0, 1不方便
# 注意！！！revalue()
library(plyr) # 為了使用 revalue()
birthwt1$smoke <- revalue(birthwt1$smoke, c("0" = "No Smoke", "1" = "Smoke"))

ggplot(birthwt1, aes(x = bwt)) + 
        geom_histogram(fill = "white", colour = "black") +
        facet_grid(smoke ~ .)

ggplot(birthwt, aes(x = bwt)) + 
        geom_histogram(fill = "white", colour = "black") + 
        facet_grid(race ~ .)

ggplot(birthwt, aes(x = bwt)) + 
        geom_histogram(fill = "white", colour = "black") + 
        facet_grid(race ~ ., scales = "free") # scales = "free" 可以單獨設定各y軸尺度

birthwt1$smoke <- factor(birthwt1$smoke)
# 注意！！！
ggplot(birthwt1, aes(x = bwt, fill = smoke)) + 
        geom_histogram(position = "identity",  # 使直方圖重疊
                       alpha = .4) # 設定透明度



## 6.3 密度曲線 geom_density(fill); geom_line(stat = "density", adjust) =====

ggplot(faithful, aes(x = waiting)) + geom_density()

ggplot(faithful, aes(x = waiting)) + 
        geom_line(stat = "density") +  # 注意！！！不想要邊線和底線可以用 geom_line()
        expand_limits(y = 0) # 擴大y軸範圍，包含0

w <- faithful$waiting
ggplot(NULL, aes(x = w)) + geom_density() # 可以傳遞向量作為參數

ggplot(faithful, aes(x = waiting)) +  # 帶寬不同，曲線的光滑程度不同
        geom_line(stat = "density", adjust = .25, colour = "red") +
        geom_line(stat = "density") + 
        geom_line(stat = "density", adjust = 2, colour = "blue")

ggplot(faithful, aes(x = waiting)) + 
        geom_density(fill = "blue", alpha = .2) + 
        xlim(35, 105) # 限制x軸範圍

ggplot(faithful, aes(x = waiting)) + 
        geom_density(fill = "blue", colour = NA, alpha = .2) + 
        geom_line(stat = "density") + 
        xlim(35, 105)

# 注意！！！同時畫density跟histogram的方法！！！
ggplot(faithful, aes(x = waiting, y = ..density..)) +  # y = ..density.. 可以縮小直方圖的尺度以配合密度曲線
        geom_histogram(fill = "cornsilk", colour = "grey60", size = .2) +
        geom_density() + 
        xlim(35, 105)



## 6.4 多組資料密度曲線 colour, facet_grid() =====

library(MASS)
birthwt1 <- birthwt
birthwt1$smoke <- factor(birthwt1$smoke) # 必須先轉化為factor
ggplot(birthwt1, aes(x = bwt, colour = smoke)) + 
        geom_density()

ggplot(birthwt1, aes(x = bwt, fill = smoke)) + 
        geom_density(alpha = .3)

ggplot(birthwt1, aes(x = bwt)) + 
        geom_density() + 
        facet_grid(smoke ~ .)

levels(birthwt1$smoke)
library(plyr)
birthwt1$smoke <- revalue(birthwt1$smoke, c("0" = "No Smoke", "1" = "Smoke"))
ggplot(birthwt1, aes(x = bwt)) + 
        geom_density() + 
        facet_grid(smoke ~ .)

ggplot(birthwt1, aes(x = bwt, y = ..density..)) + 
        geom_histogram(binwidth = 200, fill = "cornsilk", colour = "grey60", size = .2) +
        geom_density() + 
        facet_grid(smoke ~ .)



## 6.5 頻次多邊形 geom_freepoly() =====

# freqpoly 跟 histogram 非常類似，而核密度曲線只是一個估計

ggplot(faithful, aes(x = waiting)) + geom_freqpoly()
ggplot(faithful, aes(x = waiting)) + geom_freqpoly(binwidth = 4) # 控制組距
binsize <- diff(range(faithful$waiting)) / 15
ggplot(faithful, aes(x = waiting)) + geom_freqpoly(binwidth = binsize)


## 6.6 基本箱型圖（盒鬚圖） geom_boxplot() =====

library(MASS)
ggplot(birthwt, aes(x = factor(race), y = bwt)) + 
        geom_boxplot()
ggplot(birthwt, aes(x = factor(race), y = bwt)) + 
        geom_boxplot(width = .5) # 調整箱型圖寬度
ggplot(birthwt, aes(x = factor(race), y = bwt)) + 
        geom_boxplot(outlier.size = 1.5, outlier.shape = 21) # 注意！！！修改outlier的外觀
ggplot(birthwt, aes(x = 1, y = bwt)) + 
        geom_boxplot() + 
        scale_x_continuous(breaks = NULL) +  # 注意！！！只有1組，x軸不標記
        theme(axis.title.x = element_blank())


## 6.7 缺口箱型圖 notch =====

# 注意！！！notch
ggplot(birthwt, aes(x = factor(race), y = bwt)) + 
        geom_boxplot(notch = TRUE)
# 若各組缺口不重疊，表示各組中位數有差異


## 6.8 箱型圖加平均值 stat_summary(fun.y = "mean", geom = "point") =====

# 注意！！！
ggplot(birthwt, aes(x = factor(race), y = bwt)) + 
        geom_boxplot() + 
        stat_summary(fun.y = "mean", geom = "point", shape = 23, size = 3, fill = "white")


## 6.9 小提琴圖 geom_violin(trim, scale, adjust) =====

library(gcookbook)
p <- ggplot(heightweight, aes(x = sex, y = heightIn))
p + geom_violin()

p + geom_violin() +  # 同時畫小提琴圖跟箱型圖
        geom_boxplot(width = .1, fill = "black", outlier.colour = NA) +
        stat_summary(fun.y = median, geom = "point", fill = "white", shape = 21, size = 2.5)

p + geom_violin(trim = FALSE) # 注意！！！保留小提琴的尾部

# 注意！！！default情況是不同組的小提琴圖面積會相同，可以用scale = "count"調整
p + geom_violin(scale = "count")

p + geom_violin(adjust = 2) # adjust修改平滑程度
p + geom_violin(adjust = .5)


## 6.10 Wilkinson點圖 geom_dotplot(method, stackdir) =====

library(gcookbook)
countries2009 <- subset(countries, Year == 2009 & healthexp > 2000)

p <- ggplot(countries2009, aes(x = infmortality))
p + geom_dotplot()

p + geom_dotplot(binwidth = .25) + geom_rug() +  # geom_rug() 標示資料點的具體位置
        scale_y_continuous(breaks = NULL) +  # 移除y軸尺度
        theme(axis.title.y = element_blank()) # 移除座標軸標籤
# 數據在水平方向上並非均勻分佈，而是堆在他表示的資料點的中心位置

p + geom_dotplot(method = "histodot", binwidth = .25) +
        geom_rug() +
        scale_y_continuous(breaks = NULL) + 
        theme(axis.title.y = element_blank())
# 用像直方圖的固定間距分組算法

# 奇數數量與偶數數量保持一致的中心堆疊方式，設定stackdir = "center"或stackdir = "centerwhole"
p + geom_dotplot(binwidth = .25, stackdir = "center") +
        scale_y_continuous(breaks = NULL) + 
        theme(axis.title.y = element_blank())
p + geom_dotplot(binwidth = .25, stackdir = "centerwhole") +
        scale_y_continuous(breaks = NULL) + 
        theme(axis.title.y = element_blank())



## 6.11 多組資料點圖 geom_dotplot(binaxis) =====

library(gcookbook)
# sex為factor類型
# 注意！！！binaxis = "y"
ggplot(heightweight, aes(x = sex, y = heightIn)) + 
        geom_dotplot(binaxis = "y", binwidth = .5, stackdir = "center")

# 同時畫箱型圖跟點圖
ggplot(heightweight, aes(x = sex, y = heightIn)) + 
        geom_boxplot(outlier.colour = NA, width = .4) + 
        geom_dotplot(binaxis = "y", binwidth = .5, stackdir = "center", fill = NA)

# 注意！！！把箱型圖畫在點圖旁邊（移動整體y值）
ggplot(heightweight, aes(x = sex, y = heightIn)) + 
        geom_boxplot(aes(x = as.numeric(sex) + .2, group = sex), width = .25) + 
        geom_dotplot(aes(x = as.numeric(sex) - .2, group = sex), binaxis = "y", binwidth = .5, stackdir = "center") + 
        scale_x_continuous(breaks = 1:nlevels(heightweight$sex), labels = levels(heightweight$sex))
# nlevels(): 算某factor有幾個level


## 6.12 二維資料密度圖 geom_density2d() =====

p <- ggplot(faithful, aes(x = eruptions, y = waiting))
p + geom_point() + stat_density2d() # 此處等高線顏色相同

# 注意！！！..level..，把密度曲面的高度映射給等高線的顏色
p + stat_density2d(aes(colour = ..level..))

# 將密度估計映射給填充色 fill = ..density..
# 將密度估計映射給瓦片圖的透明度 alpha = ..density..
# 柵格: raster（更有效地進行渲染）
# 瓦片: tile
p + stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE)
p + geom_point() + 
        stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE)
p + stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE, h = c(.5, 5)) # h如同binwidth


# 7. 註解 =====

## 7.1 添加文本註解 annotate("text", x, y, label = ); geom_text(label = ) =====

p <- ggplot(faithful, aes(x = eruptions, y = waiting)) + geom_point()
p + annotate("text", x = 3, y = 48, label = "Group 1") + 
        annotate("text", x = 4.5, y = 66, label = "Group 2")

p + annotate("text", x = 3, y = 48, label = "Group 1", family = "serif", fontface = "italic", 
             colour = "darkred", size = 3) + 
        annotate("text", x = 4.5, y = 66, label = "Group 2", family = "serif", fontface = "italic", 
                 colour = "darkred", size = 3)

# 如果使用 geom_text() ，文本會在相同位置被遮蓋
p + annotate("text", x = 3, y = 48, label = "Group 1", alpha = .1) + 
        geom_text(x = 4.5, y = 66, label = "Group 2", alpha = .1)

p + annotate("text", x = -Inf, y = Inf, label = "Upper left", hjust = -.2, vjust = 2) + 
        annotate("text", x = mean(range(faithful$eruptions)), y = -Inf, vjust = -.4, label = "Bottom middle")
# 注意！！！Inf, -Inf 在繪圖區邊緣放置文本，hjust, vjust 調整文本相對於邊緣的位置



## 7.2 在註解中用數學方程式 annotate("text", label = ) =====

p <- ggplot(data.frame(x = c(-3, 3)), aes(x = x)) + stat_function(fun = dnorm)
p + annotate("text", x = 2, y = .3, parse = TRUE, label = "frac(1, sqrt(2 * pi)) * e ^ {-x ^ 2 / 2}")
# 加 parse = TRUE
# 內部引號閉合的每一部分都會被當作數學方程式的一個變量，要讓乘號 * 變可見要用 %*%

p + annotate("text", x = 0, y = .05, parse = TRUE, size = 4, label = "'Function: ' * y == frac(1, sqrt(2 * pi)) * e ^ {-x ^ 2 / 2}")
# ?plotmath 看更多數學方程
# ?demo(plotmath) 看數學方程的圖示


## 7.3 添加直線 geom_hline(yintercept); geom_vline(xintercept); geom_abline(intercept, slope) =====

library(gcookbook)
p <- ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + 
        geom_point()
p + geom_hline(yintercept = 60) + 
        geom_vline(xintercept = 14)
p + geom_abline(intercept = 37.4, slope = 1.75)

library(plyr)
hw_means <- ddply(heightweight, "sex", summarise, heightIn = mean(heightIn))
hw_means

p + geom_hline(aes(yintercept = heightIn, colour = sex), data = hw_means, linetype = "dashed", size = 1)


# 若 x 軸為離散型變數
pg <- ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_point()
pg + geom_vline(xintercept = 1)
pg + geom_vline(xintercept = which(levels(PlantGrowth$group) == "ctrl"))
# which(levels(PlantGrowth$group) == "ctrl") 其實就是1



## 7.4 添加線段和箭頭 annotate("segment", x, xend, y, yend, arrow = arrow()) =====

library(gcookbook)
p <- ggplot(subset(climate, Source == "Berkeley"), aes(x = Year, y = Anomaly10y)) + 
        geom_line()
# 注意！！！
p + annotate("segment", x = 1950, xend = 1980, y = -.25, yend = -.25)
library(grid)
p + annotate("segment", x = 1850, xend = 1820, y = -.8, yend = -.95, colour = "blue", size = 2, arrow = arrow()) + 
        annotate("segment", x = 1950, xend = 1980, y = -.25, yend = -.25, 
                 arrow = arrow(ends = "both", angle = 90, length = unit(.2, "cm")))
# 前者為30度箭頭，後者為90度雙箭頭
ggplot(subset(climate, Source == "Berkeley"), aes(x = Year, y = Anomaly10y)) + 
        geom_line() + 
        annotate("segment", x = 1850, xend = 1820, y = -.8, yend = -.95, colour = "blue", size = 2, arrow = arrow()) + 
        annotate("segment", x = 1950, xend = 1980, y = -.25, yend = -.25, arrow = arrow(ends = "both", angle = 90, length = unit(.2, "cm")))



## 7.5 添加矩形陰影 annotate("rect", xmin, xmax, ymin, ymax) =====

library(gcookbook)
p <- ggplot(subset(climate, Source == "Berkeley"), aes(x = Year, y = Anomaly10y)) + 
        geom_line()
p + annotate("rect", xmin = 1950, xmax = 1980, ymin = -1, ymax = 1, alpha = .1, fill = "blue")



## 7.6 高亮某一元素 scale_fill_manual(values) =====

# 注意！！！非常實用
pg <- PlantGrowth
pg$h1 <- "no"
pg$h1[pg$group == "trt2"] <- "yes" # 方法1：改值
ggplot(pg, aes(x = group, y = weight, fill = h1)) + geom_boxplot() + 
        scale_fill_manual(values = c("grey85", "#FFDDCC"), guide = FALSE)

# 方法2：逐一設定顏色
ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) + geom_boxplot() + 
        scale_fill_manual(values = c("grey85", "grey85", "#FFDDCC"), guide = FALSE)



## 7.7 添加誤差線 geom_errorbar(aes(ymin, ymax), position) =====

library(gcookbook)
ce <- subset(cabbage_exp, Cultivar == "c39")
ggplot(ce, aes(x = Date, y = Weight)) + 
        geom_bar(fill = "white", colour = "black", stat = "identity") + 
        geom_errorbar(aes(ymin = Weight - se, ymax = Weight + se), width = .2)

ggplot(ce, aes(x = Date, y = Weight)) + 
        geom_line(aes(group = 1)) + 
        geom_point(size = 4) + 
        geom_errorbar(aes(ymin = Weight - se, ymax = Weight + se), width = .2)

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(position = "dodge", stat = "identity") +  # 並列
        geom_errorbar(aes(ymin = Weight - se, ymax = Weight + se), 
                      position = "dodge", width = .2)

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(position = "dodge", stat = "identity") + 
        geom_errorbar(aes(ymin = Weight - se, ymax = Weight + se), 
                      position = position_dodge(.9), width = .2)
# position = "dodge"即為position = position_dodge()的縮寫版本

pd <- position_dodge(.3) # 注意！！！將圖形左右平移
ggplot(cabbage_exp, aes(x = Date, y = Weight, colour = Cultivar, group = Cultivar)) + 
        geom_errorbar(aes(ymin = Weight - se, ymax = Weight + se), width = .2, size = .25, colour = "black", position = pd) + 
        geom_line(position = pd) + 
        geom_point(position = pd, size = 2.5)



## 7.8 向獨立分面添加註解 geom_text(label = 變數) =====

p <- ggplot(mpg, aes(x = displ, y = hwy)) + 
        geom_point() + facet_grid(. ~ drv)
f_labels <- data.frame(drv = c("4", "f", "r"), label = c("4wd", "Front", "Rear"))
# 注意！！！
p + geom_text(x = 6, y = 40, aes(label = label), data = f_labels)
p + annotate("text", x = 6, y = 42, label = "label text")

lm_labels <- function(dat) {
        mod <- lm(hwy ~ displ, data = dat)
        formula <- sprintf("italic(y) == %.2f %+.2f * italic(x)",
                           round(coef(mod)[1], 2), round(coef(mod)[2], 2))
        r <- cor(dat$displ, dat$hwy)
        r2 <- sprintf("italic(R ^ 2) == %.2f", r ^ 2)
        data.frame(formula = formula, r2 = r2, stringsAsFactors = FALSE)
}

library(plyr)
labels <- ddply(mpg, "drv", lm_labels)
labels

p + geom_smooth(method = lm, se = FALSE) + 
        geom_text(x = 3, y = 40, aes(label = formula), data = labels, parse = TRUE, hjust = 0) + 
        geom_text(x = 3, y = 34, aes(label = r2), data = labels, parse = TRUE, hjust = 0) # 以x座標為起始

labels <- ddply(mpg, "drv", summarise, r2 = cor(displ, hwy) ^ 2)
labels$r2 <- sprintf("italic(R ^ 2) == %.2f", labels$r2)



# 8. 座標軸 =====

## 8.1 交換x軸和y軸 coord_flip(); scale_x_discrete(limits = rev(levels)) =====

ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot()
ggplot(PlantGrowth, aes(x = group, y = weight)) + 
        geom_boxplot() + coord_flip()

ggplot(PlantGrowth, aes(x = group, y = weight)) + 
        geom_boxplot() + 
        coord_flip() + 
        scale_x_discrete(limits = rev(levels(PlantGrowth$group))) # 反轉因子的排列順序



## 8.2 設置連續型座標軸的值域 scale_y_continuous(limits = c(), breaks); ylim() =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot()
p

p + ylim(0, max(PlantGrowth$weight))
ylim(0, 10)
scale_y_continuous(limits = c(0, 10)) # 兩式等價

p + ylim(0, 10) + scale_y_continuous(breaks = c(0, 5, 10))
p + scale_y_continuous(breaks = c(0, 5, 10)) + 
        ylim(0, 10) # 注意！！！同時使用會產生錯誤
p + scale_y_continuous(limits = c(0, 10), breaks = c(0, 5, 10)) # 想同時生效可加在參數

p + scale_y_continuous(limits = c(5, 6.5)) # 限制y的值域會使得部分資料被剪除
p + coord_cartesian(ylim = c(5, 6.5)) # 此法純粹縮放

p + expand_limits(y = 0) # 擴展值域



# 8.3 反轉一條連續型座標軸 scale_y_reverse() =====
# 注意！！！
ggplot(PlantGrowth, aes(x = group, y = weight)) + 
        geom_boxplot() + 
        scale_y_reverse() # y軸刻度反轉
#注意方法2！！！把值反著放也可以
ggplot(PlantGrowth, aes(x = group, y = weight)) + 
        geom_boxplot() + 
        ylim(6.5, 3.5)

ggplot(PlantGrowth, aes(x = group, y = weight)) + 
        geom_boxplot() + 
        scale_y_reverse(limits = c(8, 0))
# ylim()與scale_y_reverse()同樣不能配合，要設定在參數



## 8.4 修改類別型座標軸上項目的順序 scale_x_discrete(limits = c(, , )) =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot()
p + scale_x_discrete(limits = c("trt1", "ctrl", "trt2"))

p + scale_x_discrete(limits = c("ctrl", "trt1")) # 只展示部分類別
p + scale_x_discrete(limits = rev(levels(PlantGrowth$group))) # 反轉類別順序



## 8.5 設置x軸和y軸的縮放比例 coord_fixed() =====

library(gcookbook)
sp <- ggplot(marathon, aes(x = Half, y = Full)) + geom_point()
sp + coord_fixed() # 等長縮放座標軸

sp + coord_fixed() + 
        scale_y_continuous(breaks = seq(0, 420, 30)) + 
        scale_x_continuous(breaks = seq(0, 420, 30)) # 指定位置放刻度線

sp + coord_fixed(ratio = 1/2) +  # 讓座標軸延伸
        scale_y_continuous(breaks = seq(0, 420, 30)) + 
        scale_x_continuous(breaks = seq(0, 420, 15))



## 8.6 設置刻度線的位置 scale_y_continuous(breaks) =====

ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot()
ggplot(PlantGrowth, aes(x = group, y = weight)) + 
        geom_boxplot() + 
        scale_y_continuous(breaks = c(4, 4.25, 4.5, 5, 6, 8))
ggplot(PlantGrowth, aes(x = group, y = weight)) + 
        geom_boxplot() + 
        scale_x_discrete(limits = c("trt2", "ctrl"), breaks = "ctrl") # 用limits篩選組別



## 8.7 移除刻度線和標籤 theme(axis.text.y = element_blank(), axis.ticks = element_blank()) =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot()
p + theme(axis.text.y = element_blank()) # 移除刻度標籤

p + theme(axis.ticks = element_blank(), 
          axis.text.y = element_blank()) # 注意！！！移除兩軸刻度線
# axis.ticks 刻度線
# axis.text 刻度標籤
p + scale_y_continuous(breaks = NULL) # 移除刻度線、刻度標籤和網格線



## 8.8 修改刻度標籤的文本 scale_y_continuous(breaks, labels) =====

library(gcookbook)
hwp <- ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point()
hwp + scale_y_continuous(breaks = c(50, 56, 60, 66, 72), 
                         labels = c("Tiny", "Really\nshort", "Short", "Medium", "Tallish"))
# 注意！！！ \n 指的是換行！

# 注意！！！
footinch_formatter <- function(x) {
        foot <- floor(x / 12)
        inch <- x %% 12
        return(paste(foot, "'", inch, "\"", sep = "")) # \"加反斜線是為了區分引號字符與參數引號
}

footinch_formatter(56:64)
hwp + scale_y_continuous(labels = footinch_formatter)

timeHMS_formatter <- function(x) {
        h <- floor(x / 60)
        m <- floor(x %% 60)
        s <- round(60 * (x %% 1))
        lab <- sprintf("%02d:%02d:%02d", h, m, s)
        lab <- gsub("^00:", "", lab) # 開頭為00:則移除
        lab <- gsub("^0", "", lab) # 開頭為0則移除
        return(lab)
}

timeHMS_formatter(c(.33, 50, 51.25, 59.32, 60, 60.1, 130.23))

# 注意！！！
## library(scales)
## comma()：在數字加入逗號
## dollar()：將美元加入符號並四捨五入到最接近的美分
## percent()：乘以100並加上百分符號
## scientific()：科學記數法



## 8.9 修改刻度標籤的外觀 theme(axis.text.x = element_text()) =====

bp <- ggplot(PlantGrowth, aes(x = group, y = weight)) + 
        geom_boxplot() + 
        scale_x_discrete(breaks = c("ctrl", "trt1", "trt2"), labels = c("Control", " Treatment 1", "Treatment 2"))

bp + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) # 標籤轉90度
bp + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
bp + theme(axis.text.x = element_text(family = "Times", face = "italic", colour = "darkred", size = rel(0.9)))
# rel(0.9)指的是目前主題基礎字體大小的0.9倍



## 8.10 修改座標軸標籤的文本 scale_x_continuous(name); xlab() =====

library(gcookbook)
hwp <- ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + 
        geom_point()
hwp

hwp + xlab("Age in years") + ylab("Height in inches")

hwp + labs(x = "Age in years", y = "Height in inches")
hwp + scale_x_continuous(name = "Age in years")
hwp + scale_x_continuous(name = "Age\n(years)")



## 8.11 移除座標軸標籤 theme(axis.title.x = element_blank()) =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot()
p + theme(axis.title.x = element_blank()) # 隱藏x軸標籤

p + xlab("") # 另一種做法，但圖上還是會為x軸標籤留出空間



## 8.12 修改座標軸標籤的外觀 theme(axis.title.x = element_text()) =====

library(gcookbook)
hwp <- ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point()
hwp + theme(axis.title.x = element_text(face = "italic", colour = "darkred", size = 14))

hwp + ylab("Height\n(inches)") + 
        theme(axis.title.y = element_text(angle = 0, face = "italic", size = 14))

hwp + ylab("Height\n(inches)") + 
        theme(axis.title.y = element_text(angle = 90, face = "italic", 
                                          colour = "darkred", size = 14))



## 8.13 沿座標軸顯示直線 theme(panel.border = element_blank(), axis.line = element_line()) =====

library(gcookbook)
p <- ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point()
p
p + theme(axis.line = element_line(colour = "red")) # 看不出差在哪
p + theme_bw() + 
        theme(panel.border = element_blank(), 
              axis.line = element_line(colour = "black"))
p + theme_bw() + 
        theme(panel.border = element_blank(), 
              axis.line = element_line(colour = "black", size = 4, lineend = "square"))



## 8.14 使用對數座標軸 scale_x_log10() =====

library(MASS)
p <- ggplot(Animals, aes(x = body, y = brain, label = rownames(Animals))) + 
        geom_text(size = 3)
p
p + scale_x_log10() + scale_y_log10()
p + scale_x_log10(breaks = 10 ^ (-1:5)) + scale_y_log10(breaks = 10 ^ (0:3))

library(scales)
p + scale_x_log10(breaks = 10 ^ (-1:5), labels = trans_format("log10", math_format(10 ^ .x))) + 
        scale_y_log10(breaks = 10 ^ (0:3), labels = trans_format("log10", math_format(10 ^ .x))) # 用scale套件的trans_format()函數即可

# 也可以先轉換再畫圖
ggplot(Animals, aes(x = log10(body), y = log10(brain), label = rownames(Animals))) + geom_text(size = 3)
# 注意！！！難！！！
library(scales)
p + scale_x_continuous(trans = log_trans(),
                       breaks = trans_breaks("log", function(x) exp(x)),
                       labels = trans_format("log", math_format(e ^ .x))) + 
        scale_y_continuous(trans = log2_trans(),
                           breaks = trans_breaks("log2", function(x) 2 ^ x),
                           labels = trans_format("log2", math_format(2 ^ .x)))

library(gcookbook)
ggplot(aapl, aes(x = date, y = adj_price)) + geom_line()
# 對數轉換對金融資料很有用
ggplot(aapl, aes(x = date, y = adj_price)) + 
        geom_line() + 
        scale_y_log10(breaks = c(2, 10, 50, 250))


## 8.15 為對數座標軸添加刻度 labels = trans_format() =====

library(MASS)
library(scales)
ggplot(Animals, aes(x = body, y = brain, label = rownames(Animals))) + 
        geom_text(size = 3) + 
        annotation_logticks() + 
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x),
                      labels = trans_format("log10", math_format(10 ^ .x))) + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ x),
                      labels = trans_format("log10", math_format(10 ^ .x)))

ggplot(Animals, aes(x = body, y = brain, label = rownames(Animals))) + 
        geom_text(size = 3) + 
        annotation_logticks() + 
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ x),
                      labels = trans_format("log10", math_format(10 ^ .x)),
                      minor_breaks = log10(5) + -2:5) + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ x), 
                      labels = trans_format("log10", math_format(10 ^ .x)), 
                      minor_breaks = log10(5) + -1:3) + 
        coord_fixed() + 
        theme_bw()



## 8.16 繪製環狀圖形 geom_histogram(origin) + coord_polar(); scale_x_continuous(minor_breaks = ); %+% =====

library(gcookbook)
head(wind)
ggplot(wind, aes(x = DirCat, fill = SpeedCat)) + 
        geom_histogram(binwidth = 15, origin = -7.5) + 
        coord_polar() + 
        scale_x_continuous(limits = c(0, 360))
# 注意！！！minor_breaks
ggplot(wind, aes(x = DirCat, fill = SpeedCat)) + 
        geom_histogram(binwidth = 15, origin = -7.5, colour = "black", size = .25) + 
        guides(fill = guide_legend(reverse = TRUE)) + 
        coord_polar() + 
        scale_x_continuous(limits = c(0, 360), breaks = seq(0, 360, by = 45), 
                           minor_breaks = seq(0, 360, by = 15)) + 
        scale_fill_brewer()

ggplot(wind, aes(x = DirCat, fill = SpeedCat)) + 
        geom_histogram(binwidth = 15, origin = -7.5, colour = "black", size = .25) + 
        guides(fill = guide_legend(reverse = TRUE)) + 
        coord_polar(start = -45 * pi / 180) +  # 注意！！！設置圖形起始角度
        scale_x_continuous(limits = c(0, 360), breaks = seq(0, 360, by = 45), 
                           minor_breaks = seq(0, 360, by = 15)) + 
        scale_fill_brewer()


md <- data.frame(deaths = as.numeric(mdeaths), 
                 month = as.numeric(cycle(mdeaths)))
library(plyr)
md <- ddply(md, "month", summarise, deaths = mean(deaths))
md

p <- ggplot(md, aes(x = month, y = deaths)) + 
        geom_line() + 
        scale_x_continuous(breaks = 1:12)
p + coord_polar() # 注意！！！線幾何極座標圖
p + coord_polar() + ylim(0, max(md$deaths))

p + coord_polar() + ylim(0, max(md$deaths)) + 
        xlim(0, 12) # 注意！！！設置界限0來解決1跟12在同樣角度的問題


mdx <- md[md$month == 12, ]
mdx$month <- 0 # 注意！！！新增一筆0月資料，資料跟尾巴（12月）一樣，就可以將起點終點連起來
mdnew <- rbind(mdx, md)
# 注意！！！這招超神！！！
p %+% mdnew + coord_polar() + ylim(0, max(md$deaths)) # 舊圖設定 %+% 新資料 ---> 畫原本的圖，但使用新資料



## 8.17 在座標軸上使用日期 =====



## 8.18 在座標軸上使用相對時間 =====



# 9. 整體外觀 =====

## 9.1 設置圖形標題 ggtitle(); labs(title) =====

library(gcookbook)
p <- ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point()
p + ggtitle("Age and Height of Schoolchildren")
p + ggtitle("Age and Height\nof Schoolchildren")
p + labs(title = "Age and Height\nof Schoolchildren") # ggtitle("") 與 labs(title = "") 等價

p + ggtitle("Age and Height of Schoolchildren") + 
        theme(plot.title = element_text(vjust = -2.5)) # 似乎無效
p + annotate("text", x = mean(range(heightweight$ageYear)), y = Inf,
             label = "Age and Height of Schoolchildren", vjust = 1.5, size = 6)



## 9.2 修改文本外觀 theme(); annotate(); geom_text() =====

library(gcookbook)
p <- ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point()
p + theme(axis.title.x = element_text(size = 16, 
                                      lineheight = .9, 
                                      family = "Times", 
                                      face = "bold.italic", 
                                      colour = "red"))
p + ggtitle("Age and Height\nof Schoolchildren") + 
        theme(plot.title = element_text(
                size = rel(1.5), lineheight = .9, family = "Times", 
                face = "bold.italic", colour = "red")) # 針對主題元素

p + annotate("text", x = 15, y = 53, label = "Some text", size = 7, 
             family = "Times", fontface = "bold.italic", colour = "red")
p + geom_text(aes(label = weightLb), size = 4, family = "Times", colour = "red") # 針對文本幾何對象



## 9.3 使用主題 theme_grey(); theme_set() =====

library(gcookbook)
p <- ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point()
p + theme_grey() # default
p + theme_grey(base_size = 16, base_family = "Times") # 設定主題的基本字體和大小
theme_set(theme_bw()) # 設定default主題
p
theme_set(theme_grey()) # 設定重回theme_grey()



## 9.4 修改主題元素的外觀 theme(panel., axis., legend); facet_grid(strip.) =====

library(gcookbook)
p <- ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + geom_point()

p + theme(
        panel.grid.major = element_line(colour = "red"),
        panel.grid.minor = element_line(colour = "red", linetype = "dashed", size = .2),
        panel.background = element_rect(fill = "lightblue"),
        panel.border = element_rect(colour = "blue", fill = NA, size = 2))

p + ggtitle("Plot title here") + 
        theme(
                axis.title.x = element_text(colour = "red", size = 14), # 軸標題
                axis.text.x = element_text(colour = "blue"), # 刻度
                axis.title.y = element_text(colour = "red", size = 14, angle = 90),
                axis.text.y = element_text(colour = "blue"),
                plot.title = element_text(colour = "red", size = 20, face = "bold"))

p + theme(
        legend.background = element_rect(fill = "grey85", colour = "red", size = 1),
        legend.title = element_text(colour = "blue", face = "bold", size = 14),
        legend.text = element_text(colour = "red"),
        legend.key = element_rect(colour = "blue", size = .25)
)

p + facet_grid(sex ~ .) + theme(
        strip.background = element_rect(fill = "pink"),
        strip.text.y = element_text(size = 14, angle = -90, face = "bold")
)
# 注意！！！想套用現成的主題但又想要微調，微調的部分必須加在theme_bw()之後！
p + theme(axis.title.x = element_text(colour = "red")) + theme_bw() # 否則會被還原
p + theme_bw() + theme(axis.title.x = element_text(colour = "red", size = 12))



## 9.5 創建自定義主題 theme() =====
# 注意！！！
library(gcookbook)
mytheme <- theme_bw() + 
        theme(text = element_text(colour = "red"), 
              axis.title = element_text(size = rel(1.25)))
p <- ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point()
p + mytheme



## 9.6 隱藏網格線 theme(panel.grid.major/minor = element_blank()) =====

library(gcookbook)
p <- ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point()

p + theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

p + theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())

p + theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())



# 10. 圖例 =====

## 10.1 移除圖例 guides(fill = FALSE); scale_fill_discrete(guide = FALSE); theme(legend.position = "none") =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) + geom_boxplot()
p

p + guides(fill = FALSE) # 注意！！！移除某特定圖例
p + scale_fill_discrete(guide = FALSE)
# 移除標註fill的圖例
# 另外還有scale_colour_discrete(), scale_shape_discrete()
p + theme(legend.position = "none") # 移除「所有」圖例



## 10.2 修改圖例的位置 theme(legend.position = c(, ), legend.justification = c(, )) =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) + geom_boxplot() + 
        scale_fill_brewer(palette = "Pastel2")
p + theme(legend.position = "top")

p + theme(legend.position = c(1, 0), legend.justification = c(1, 0)) # 右下角
p + theme(legend.position = c(1, 1), legend.justification = c(1, 1)) # 右上角

p + theme(legend.position = c(.85, .2)) + 
        theme(legend.background = element_rect(fill = "white", colour = "black")) # 白底黑邊

p + theme(legend.position = c(.85, .2)) + 
        theme(legend.background = element_blank()) +  # 移除圖例背景
        theme(legend.key = element_blank()) # 移除圖例邊框



## 10.3 修改圖例項目的順序 scale_fill_discrete(limits) =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) + geom_boxplot()
p
# 注意！！！
p + scale_fill_discrete(limits = c("trt1", "trt2", "ctrl")) # 只能修改圖例的順序

p + scale_x_discrete(limits = c("trt1", "trt2", "ctrl")) + 
        scale_fill_discrete(limits = c("trt1", "trt2", "ctrl")) # 要改x軸順序要用scale_x_discrete()

p + scale_fill_grey(start = .5, end = 1, limits = c("trt1", "trt2", "ctrl")) # 用不同標度scale_fill_xxx()
p + scale_fill_brewer(palette = "Pastel2", limits = c("trt1", "trt2", "ctrl"))



## 10.4 反轉圖例項目的順序 guides(fill = guide_legend(reverse = TRUE); sale_fill_hue(guide_legend(reverse = TRUE)) =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) + geom_boxplot()
p
p + guides(fill = guide_legend(reverse = TRUE)) # x軸順序沒變
p + scale_fill_hue(guide_legend(reverse = TRUE)) # default情況下，scale_fill_discrete()和scale_fill_hue()等價



## 10.5 修改圖例標題 labs(fill = "", colour = "", size = "", shape = "");  =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) + geom_boxplot()
p
p + labs(fill = "Condition")
library(gcookbook)
# 注意！！！scale_size_continuous()
hw <- ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + 
        geom_point(aes(size = weightLb)) + 
        scale_size_continuous(range = c(1, 4))
hw
hw + labs(colour = "Male/Female", size = "Weight\n(pounds)")

hw1 <- ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, colour = sex)) + 
        geom_point()
hw1 + labs(shape = "Male/Female")
hw1 + labs(shape = "Male/Female", colour = "Male/Female") # 注意！！！sex同時映射到shape和colour，因此只有一個圖例



## 10.6 修改圖例標題的外觀 theme(legend.title = element_text()); guides(fill = guide_legend(title.theme = element_text())) =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) + 
        geom_boxplot()
p + theme(legend.title = element_text(face = "italic", 
                                      family = "Times", colour = "red", size = 14))
p + guides(fill = guide_legend(title.theme = element_text(face = "italic", family = "times", colour = "red", size = 14)))
# 麻煩又囉唆的方法2



## 10.7 移除圖例標題 guides(fill = guide_legend(title = NULL)) =====

ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) + geom_boxplot() + 
        guides(fill = guide_legend(title = NULL))
# scale_fill_hue(guide = guide_legend(title = NULL))



## 10.8 修改圖例標籤 scale_fill/shape/colour_discrete(labels) =====

library(gcookbook)
p <- ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) + geom_boxplot()
p + scale_fill_discrete(labels = c("Control", "Treatment 1", "Treatment 2")) # x軸必須透過scale_x_discrete()來改變
# 注意！！！start, end代表顏色的range
p + scale_fill_grey(start = .5, end = 1, labels = c("Control", "Treatment 1", "Treatment 2"))

p + scale_fill_discrete(limits = c("trt1", "trt2", "ctrl"),
                        labels = c("Treatment 1", "Treatment 2", "Control"))

p <- ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, colour = sex)) + 
        geom_point()
p
# 注意！！！
p + scale_shape_discrete(labels = c("Female", "Male")) # 一個變量映射兩個圖形屬性，要同時修改兩種屬性的標籤

p + scale_shape_discrete(labels = c("Female", "Male")) + 
        scale_colour_discrete(labels = c("Female", "Male"))



## 10.9 修改圖例標籤的外觀 theme(legend.text = element_text()); guides(fill = guide_legend(label.theme = element_text())) =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) + geom_boxplot()
p + theme(legend.text = element_text(
        face = "italic", family = "Times", colour = "red", size = 14))

p + guides(fill = guide_legend(label.theme = element_text(face = "italic", family = "Times", colour = "red", size = 14)))
# somehow, it doesn't work!



## 10.10 使用含多行文本的標籤 \n =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) + 
        geom_boxplot()
# 注意！！！
p + scale_fill_discrete(labels = c("Control", "Type 1\ntreatment", "Type 2\ntreatment"))
# 注意！！！legend.key.height
library(grid) # 使用unit()增加圖例說明的高度
p + scale_fill_discrete(labels = c("Control", "Type 1\ntreatment", "Type 2\ntreatment")) + 
        theme(legend.text = element_text(lineheight = .8), 
              legend.key.height = unit(1, "cm"))



# 11. 分面 =====

## 11.1 使用分面將數據分割繪製到子圖中 facet_grid(); facet_wrap() =====

p <- ggplot(mpg, aes(x = displ, y = hwy)) + geom_point()
p + facet_grid(drv ~ .) # 縱向
p + facet_grid(. ~ cyl) # 橫向
p + facet_grid(drv ~ cyl) # 同時

p + facet_wrap(~ class) # 依次橫向排列並換行
p + facet_wrap(~ class, nrow = 2)
p + facet_wrap(~ class, ncol = 4)



## 11.2 在不同座標軸下使用分面 facet_grid(scales) =====
# 注意！！！
p <- ggplot(mpg, aes(x = displ, y = hwy)) + geom_point()
p + facet_grid(drv ~ cyl, scales = "free_y") # y軸刻度範圍各自不同
p + facet_grid(drv ~ cyl, scales = "free") # 兩軸都自由



## 11.3 修改分面的文本標籤 facet_grid(labeller = ) =====

mpg2 <- mpg
mpg2$drv <- factor(mpg2$drv)
# 注意！！！
levels(mpg2$drv)[levels(mpg2$drv) == "4"] <- "4wd"
levels(mpg2$drv)[levels(mpg2$drv) == "f"] <- "Front"
levels(mpg2$drv)[levels(mpg2$drv) == "r"] <- "Rear"
ggplot(mpg2, aes(x = displ, y = hwy)) + geom_point() + facet_grid(drv ~ .)

ggplot(mpg2, aes(x = displ, y = hwy)) + 
        geom_point() + 
        facet_grid(drv ~ ., labeller = label_both) # facet_wrap()無法

mpg3 <- mpg
mpg3$drv <- factor(mpg3$drv)
levels(mpg3$drv)[levels(mpg3$drv) == "4"] <- "4^{wd}"
levels(mpg3$drv)[levels(mpg3$drv) == "f"] <- "- Front %.% e^{pi * i}"
levels(mpg3$drv)[levels(mpg3$drv) == "r"] <- "4^{wd} - Front"
ggplot(mpg3, aes(x = displ, y = hwy)) + geom_point() + 
        facet_grid(drv ~ ., labeller = label_parsed) # 輸入字符串當貼標



## 11.4 修改分面標籤和標題的外觀 theme(strip.) =====

library(gcookbook)
# 注意！！！element_rect()
ggplot(cabbage_exp, aes(x = Cultivar, y = Weight)) + geom_bar(stat = "identity") +
        facet_grid(. ~ Date) +
        theme(strip.text = element_text(face = "bold", size = rel(1.5)), 
              strip.background = element_rect(fill = "lightblue", colour = "black", size = 1))



# 12. 配色 =====

## 12.1 設置對象的顏色 fill, colour =====

ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point(colour = "red")
library(MASS)
ggplot(birthwt, aes(x = bwt)) + geom_histogram(fill = "red", colour = "black")



## 12.2 將變量映射到顏色上 fill, colour =====

library(gcookbook)
ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(colour = "black", position = "dodge", stat = "identity")

ggplot(cabbage_exp, aes(x = Date, y = Weight)) + 
        geom_bar(aes(fill = Cultivar), colour = "black",  # 等價做法
                 position = "dodge", stat = "identity")

ggplot(mtcars, aes(x = wt, y = mpg, colour = cyl)) + geom_point()
ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point(aes(colour = cyl)) # 等價
str(cabbage_exp)

ggplot(mtcars, aes(x = wt, y = mpg, colour = factor(cyl))) + geom_point()

m <- mtcars
m$cyl <- factor(m$cyl)
ggplot(m, aes(x = wt, y = mpg, colour = cyl)) + geom_point()



## 12.3 對離散型變量使用不同的調色板 scale_fill_ =====

library(gcookbook)
p <- ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup)) + 
        geom_area()
p
p + scale_fill_discrete()
p + scale_fill_hue()
p + scale_fill_brewer()

h <- ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + 
        geom_point()
h
h + scale_colour_hue(l = 45) # 亮度default為65(0 ~ 100)，調整l可改亮度

library(RColorBrewer)
display.brewer.all()

p + scale_fill_brewer(palette = "Oranges")
p + scale_fill_grey()
p + scale_fill_grey(start = .7, end = 0)



## 12.4 對離散型變量使用自定義調色板 scale_colour_manual(values = ) =====

library(gcookbook)
h <- ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + geom_point()
h + scale_colour_manual(values = c("red", "blue"))
h + scale_colour_manual(values = c("#CC6666", "#7777DD"))

levels(heightweight$sex)
h + scale_colour_manual(values = c(m = "blue", f = "red"))



## 12.5 使用色盲友好式的調色板 =====

library(gcookbook)
p <- ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup)) + geom_area()
cb_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")
p + scale_fill_manual(values = cb_palette)



## 12.6 對連續型變量使用自定義調色板 scale_colour_gradient(); scale_colour_gradient2(); scale_colour_gradientn() =====

library(gcookbook)
p <- ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = weightLb)) + 
        geom_point(size = 3)
p

p + scale_colour_gradient(low = "black", high = "white")

library(scales)
p + scale_colour_gradient2(low = muted("red"), mid = "white", high = muted("blue"), midpoint = 110)

p + scale_colour_gradientn(colours = c("darkred", "orange", "yellow", "white"))



## 12.7 根據數值設定陰影顏色 scale_fill_manual(values = c()) =====

library(gcookbook)
cb <- subset(climate, Source == "Berkeley")
cb$valence[cb$Anomaly10y >= 0] <- "pos"
cb$valence[cb$Anomaly10y < 0] <- "neg"
cb

ggplot(cb, aes(x = Year, y = Anomaly10y)) + 
        geom_area(aes(fill = valence)) +
        geom_line() + 
        geom_hline(yintercept = 0)
# 注意！！！
interp <- approx(cb$Year, cb$Anomaly10y, n = 1000) # 改善0附近的凌亂陰影

cbi <- data.frame(Year = interp$x, Anomaly10y = interp$y)
cbi$valence[cbi$Anomaly10y >= 0] <- "pos"
cbi$valence[cbi$Anomaly10y < 0] <- "neg"
ggplot(cbi, aes(x = Year, y = Anomaly10y)) + 
        geom_area(aes(fill = valence), alpha = .4) + 
        geom_line() + 
        geom_hline(yintercept = 0) + 
        scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"), guide = FALSE) + 
        scale_x_continuous(expand = c(0, 0))



# 13. 其他圖形 =====

## 13.1 相關係數矩陣圖 corrplot(, method) =====

mcor <- cor(mtcars)
round(mcor, digits = 2)
# install.packages("corrplot")
# 注意！！！
library(corrplot)
corrplot(mcor)
corrplot(mcor, method = "shade", shade.col = NA, tl.col = "black", tl.srt = 45)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# corrplot(mcor, method = "shade", shade.col = NA, tl.col = "black", 
#         tl.srt = 45, col = col(200), addCoef.col = "black", 
#         cl.pos = "no", order = "AOE") # Error



## 13.2 繪製函數曲線 stat_function(fun = ) =====

p <- ggplot(data.frame(x = c(-3, 3)), aes(x = x))
p + stat_function(fun = dnorm)

p + stat_function(fun = dt, args = list(df = 2))
myfun <- function(xvar) {
        1 / (1 + exp(-xvar + 10))
}
ggplot(data.frame(x = c(0, 20)), aes(x = x)) + stat_function(fun = myfun) # default會畫出x範圍內的101個點



## 13.3 在函數曲線下添加陰影 =====

dnorm_limit <- function(x) {
        y <- dnorm(x)
        y[x < 0 | x > 2] <- NA # 超出範圍變NA
        return(y)
}
p <- ggplot(data.frame(x = c(-3, 3)), aes(x = x))
p + stat_function(fun = dnorm_limit, geom = "area", fill = "blue", alpha = .2) +  # 畫陰影
        stat_function(fun = dnorm) # 畫線
# 注意！！！這是一個「製造函數的函數」！
limitRange <- function(fun, min, max) {
        function(x) { 
                y <- fun(x)
                y[x < min | x > max] <- NA
                return(y)
                }
}
dlimit <- limitRange(dnorm, 0, 2)
dlimit(-2:4)

p + stat_function(fun = dnorm) + 
        stat_function(fun = limitRange(dnorm, 0, 2), geom = "area", fill = "#2E4A71", alpha = .2, n = 200)



## 13.4 繪製網絡圖（igraph） graph(c()) =====

# install.packages("igraph")
library(igraph)
gd <- graph(c(1,2, 2,3, 2,4, 1,4, 5,5, 3,6))
plot(gd)
gu <- graph(c(1,2, 2,3, 2,4, 1,4, 5,5, 3,6), directed = FALSE)
plot(gu, vertex.label = NA)

str(gd) # 圖的結構
str(gu)

set.seed(229)
plot(gu)

library(gcookbook)
madmen2

g <- graph.data.frame(madmen2, directed = TRUE)
par(mar = c(0, 0, 0, 0))
plot(g, layout = layout.fruchterman.reingold, vertex.size = 8, edge.arrow.size = .5, vertex.label = NA)
# Fruchterman-Reingold佈局算法

g <- graph.data.frame(madmen, directed = FALSE)
par(mar = c(0, 0, 0, 0))
plot(g, layout = layout.circle, vertex.size = 8, vertex.label = NA)
# 無向圖



## 13.5 在網絡圖中使用文本標籤 plat(layout = ) =====

library(igraph)
library(gcookbook)
m <- madmen[1:nrow(madmen) %% 2 == 1, ] # 刪除偶數行
g <- graph.data.frame(m, directed = FALSE)
V(g)$name # 節點名稱
plot(g, layout = layout.fruchterman.reingold, 
     vertex.size = 4, # 縮小節點
     vertex.label = V(g)$name, # 設置標籤
     vertex.label.cex = .8, # 小號字體
     vertex.label.dist = .4, # 節點和標籤位置錯開
     vertex.label.color = "black")

# 另一種參數設定方式
V(g)$size <- 4
V(g)$label <- V(g)$name
V(g)$label.cex <- .8
V(g)$label.dist <- .4
V(g)$label.color <- "black"

g$layout <- layout.fruchterman.reingold
plot(g)


E(g)
E(g)[c(2, 11, 19)]$label <- "M" # 將線賦值
E(g)$color <- "grey70"
E(g)[c(2, 11, 19)]$color <- "red"
plot(g)

# ?igraph.plotting



## 13.6 如何繪製熱圖 geom_tile(); geom_raster() =====

str(presidents)
pres_rating <- data.frame(
        rating = as.numeric(presidents),
        year = as.numeric(floor(time(presidents))),
        quarter = as.numeric(cycle(presidents))
)
pres_rating
# 注意！！！
p <- ggplot(pres_rating, aes(x = year, y = quarter, fill = rating))
p + geom_tile() # 遺失值會用灰色代替
p + geom_raster() # 效率較高

p + geom_tile() + 
        scale_x_continuous(breaks = seq(1940, 1976, by = 4)) + 
        scale_y_reverse() + 
        scale_fill_gradient2(midpoint = 50, mid = "grey70", limits = c(0, 100)) # 調色



## 13.7 繪製三維散佈圖 plot3d(); segments3d(); axes3d(); mtext3d() =====

# install.packages("rgl")
library(rgl)
plot3d(mtcars$wt, mtcars$disp, mtcars$mpg, type = "s", size = .75, lit = FALSE)

interleave <- function(v1, v2) as.vector(rbind(v1, v2))
plot3d(mtcars$wt, mtcars$disp, mtcars$mpg, 
       xlab = "Weight", ylab = "Displacement", zlab = "MPG",
       size = .75, type = "s", lit = FALSE)

segments3d(interleave(mtcars$wt, mtcars$wt),
           interleave(mtcars$disp, mtcars$disp),
           interleave(mtcars$mpg, min(mtcars$mpg)),
           alpha = .4, col = "blue")

plot3d(mtcars$wt, mtcars$disp, mtcars$mpg,
       xlab = "", ylab = "", zlab = "",
       axes = FALSE,
       size = .75, type = "s", lit = FALSE)

segments3d(interleave(mtcars$wt, mtcars$wt),
           interleave(mtcars$disp, mtcars$disp),
           interleave(mtcars$mpg, min(mtcars$mpg)),
           alpha = .4, col = "blue")

rgl.bbox(color = "grey50",
         emission = "grey50",
         xlen = 0, ylen = 0, zlen = 0)
rgl.material(color = "black")
axes3d(edges = c("x--", "y+-", "z--"),
       ntick = 6, cex = .75)
mtext3d("Weight", edg = "x--", line = 2)
mtext3d("Displacement", edge = "y+-", line = 3)
mtext3d("MPG", edge = "z--", line = 3)



## 13.8 在三維圖上添加預測曲面 surface3d() =====

# 這一段都不work，重做！！！
predictgrid <- function(model, xvar, yvar, zvar, res = 16, type = NULL){
        xrnage <- range(model$model[[xvar]])
        yrange <- range(model$model[[yvar]])
        newdata <- expand.grid(x = seq(xrange[1], xrange[2], length.out = res),
                               y = seq(yrange[1], yrange[2], length.out = res))
        names(newdata) <- c(xvar, yvar)
        newdata[[zvar]] <- predict(model, newdata = newdata, type = type)
        newdata
}
df2mat <- function(p, xvar = NULL, yvar = NULL, zvar = NULL){
        if(is.null(xvar)) xvar <- names(p)[1]
        if(is.null(yvar)) yvar <- names(p)[2]
        if(is.null(zvar)) zvar <- names(p)[3]
        x <- unique(p[[xvar]])
        y <- unique(p[[yvar]])
        z <- matrix(p[[zvar]], nrow = length(y), ncol = length(x))
        m <- list(x, y, z)
        names(m) <- c(xvar, yvar, zvar)
        m
}
interleave <- function(v1, v2) as.vector(rbind(v1, v2))
library(rgl)
m <- mtcars
mod <- lm(mpg ~ wt + disp + wt:disp, data = m)
m$pred_mpg <- predict(mod)
mpgrid_df <- predictgrid(mod, "wt", "disp", "mpg")
mpgrid_list <- df2mat(mpgrid_df)
plot3d(m$wt, m$disp, m$mpg, type = "s", size = .5, lit = FALSE)
spheres3d(m$wt, m$disp, m$pred_mpg, alpha = .4, type = "s", size = .5, lit = FALSE)
segments3d(interleave(m$wt, m$wt),
           interleave(m$disp, m$disp),
           interleave(m$mpg, m$pred_mpg),
           alpha = .4, col = "red")
surface3d(mpgrid_list$wt, mpgrid_list$disp, mpgrid_list$mpg,
          alpha = .4, front = "lines", back = "lines")

# 修改圖形外觀
plot3d(mtcars$wt, mtcars$disp, mtcars$mpg,
       xlab = "", ylab = "", zlab = "",
       axes = FALSE, size = .5, type = "s", lit = FALSE)
spheres3d(m$wt, m$disp, m$pred_mpg, alpha = .4, type = "s", size = .5, lit = FALSE)
segments3d(interleave(m$wt, m$wt),
           interleave(m$disp, m$disp),
           interleave(m$mpg, m$pred_mpg),
           alpha = .4, col = "red")
surface3d(mpgrid_list$wt, mpgrid_list$disp, mpgrid_list$mpg, 
          alpha = .4, front = "lines", back = "lines")
rgl.bbox(color = "grey50",
         emission = "grey50",
         xlen = 0, ylen = 0, zlen = 0)
rgl.material(color = "black")
axes3d(edges = c("x--", "y+-", "z--"),
       ntick = 6,
       cex = .75)
mtext3d("Weight", edge = "x--", line = 2)
mtext3d("Displacement", edge = "y+-", line = 3)
mtext3d("MPG", edge = "z--", line = 3)



## 13.9 



## 13.10



## 13.11 繪製譜系圖(Dendrogram) plot(hclust()) =====

library(gcookbook)
c2 <- subset(countries, Year == 2009)
c2 <- c2[complete.cases(c2), ]
set.seed(201)
c2 <- c2[sample(1:nrow(c2), 25), ]
c2
rownames(c2) <- c2$Name
c2 <- c2[, 4:7]
c2
c3 <- scale(c2) # 做標準化，因為GDP資料值比infmortality大非常多，不標準化會使某些變數佔主導地位
c3
hc <- hclust(dist(c3))
plot(hc)
plot(hc, hang = -1)

# ?hclust


## 13.12 繪製向量場 geom_segment() =====

library(gcookbook)
head(isabel)
islice <- subset(isabel, z == min(z))
ggplot(islice, aes(x = x, y = y)) + 
        geom_segment(aes(xend = x + vx / 50, yend = y + vy / 50),
                     size = 0.25)

every_n <- function(x, by = 2){
        x <- sort(x)
        x[seq(1, length(x), by = by)]
}
keepx <- every_n(unique(isabel$x), by = 4)
keepy <- every_n(unique(isabel$y), by = 4)
islicesub <- subset(islice, x %in% keepx & y %in% keepy)
library(grid)
ggplot(islicesub, aes(x = x, y = y)) + 
        geom_segment(aes(xend = x + vx / 50, yend = y + vy / 50),
                     arrow = arrow(length = unit(0.1, "cm")), size = .25)

# 未完


## 13.13 繪製QQ圖 qqnorm(); qqline() =====

library(gcookbook)
qqnorm(heightweight$heightIn)
qqline(heightweight$heightIn)

qqnorm(heightweight$ageYear)
qqline(heightweight$ageYear)


## 13.14 繪製經驗累積分佈函數圖 ggplot() + stat_ecdf() =====

library(gcookbook)
ggplot(heightweight, aes(x = heightIn)) + stat_ecdf()
ggplot(heightweight, aes(x = ageYear)) + stat_ecdf()


## 13.15 創建馬賽克圖 mosaic() =====

UCBAdmissions
ftable(UCBAdmissions) # 平鋪後的列聯表

dimnames(UCBAdmissions)

# install.packages("vcd")
library(vcd)

mosaic(~ Admit + Gender + Dept, data = UCBAdmissions) # 依序分割數據
# 不同的分割順序會揭露不一樣的訊息

mosaic(~ Dept + Gender + Admit, data = UCBAdmissions, 
       highlighting = "Admit", highlighting_fill = c("lightblue", "pink"),
       direction = c("v", "h", "v"))

mosaic(~ Dept + Gender + Admit, data = UCBAdmissions, 
       highlighting = "Admit", highlighting_fill = c("lightblue", "pink"), 
       direction = c("v", "v", "h")) # 使用不同分割方向

mosaic(~ Dept + Gender + Admit, data = UCBAdmissions, 
       highlighting = "Admit", highlighting_fill = c("lightblue", "pink"), 
       direction = c("v", "h", "h"))

# 此為辛普森悖論的案例



## 13.16 繪製餅圖 pie() =====

library(MASS)
fold <- table(survey$Fold)
fold
pie(c(99, 18, 120), labels = c("L on R", "Neither", "R on L"))

# 圓餅圖非常容易被批判，不如考慮使用長條圖，但圓餅圖的優勢是人人都能判讀



## 13.17 創建地圖 =====

world_map <- map_data("world")
world_map
sort(unique(world_map$region)) # 看有沒有單獨地圖資料

euro <- map_data("world", region = c("UK", "France", "Netherlands", "Belgium"))
ggplot(euro, aes(x = long, y = lat, group = group, fill = region)) + 
        geom_polygon(colour = "black") + 
        scale_fill_brewer(palette = "Set2") + 
        scale_y_continuous(limits = c(40, 60)) + 
        scale_x_continuous(limits = c(-25, 25))

nz1 <- map_data("world", region = "New Zealand")
nz1 <- subset(nz1, long > 0, lat > -48)
ggplot(nz1, aes(x = long, y = lat, group = group)) + geom_path()


map()
?mappproject



## 13.18 繪製等值區域圖(Choropleth Map) ggplot() + geom_polyon() + coord_map() ===== 

crimes <- data.frame(state = tolower(rownames(USArrests)), USArrests)
crimes
library(maps)
states_map <- map_data("state")
crime_map <- merge(states_map, crimes, by.x = "region", by.y = "state")
head(crime_map)
library(plyr)
crime_map <- arrange(crime_map, group, order)
head(crime_map)
ggplot(crime_map, aes(x = long, y = lat, group = group, fill = Assault)) + 
        geom_polygon(colour = "black") + 
        coord_map("polyconic")

ggplot(crimes, aes(map_id = state, fill = Assault)) + 
        geom_map(map = states_map, colour = "black") + 
        scale_fill_gradient2(low = "#559999", mid = "grey90", high = "#BB650B", 
                             midpoint = median(crimes$Assault)) + 
        expand_limits(x = states_map$long, y = states_map$lat) + 
        coord_map("polyconic")
qa <- quantile(crimes$Assault, c(0, .2, .4, .6, .8, 1.0))
qa

crimes$Assault_q <- cut(crimes$Assault, qa, labels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%"),
                        include.lowest = TRUE)
crimes
pal <- colorRampPalette(c("#559999", "grey80", "#BB650B"))(5)
pal
ggplot(crimes, aes(map_id = state, fill = Assault_q)) + 
        geom_map(map = states_map, colour = "black") + 
        expand_limits(x = states_map$long, y = states_map$lat) + 
        coord_map("polyconic") + 
        labs(fill = "Assault Rate\nPercentile")

ggplot(crimes, aes(map_id = state, fill = Assault)) + 
        geom_map(map = states_map) + 
        expand_limits(x = states_map$long, y = states_map$lat) + 
        coord_map("polyconic")



## 13.19 創建空白背景的地圖 =====

theme_clean <- function(base_size = 12) {
        require(grid)
        theme_grey(base_size) %+replace%
                theme(
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        panel.background = element_blank(),
                        panel.grid = element_blank(),
                        axis.ticks.length = unit(0, "cm"),
                        axis.ticks.margin = unit(0, "cm"),
                        panel.margin = unit(0, "lines"),
                        plot.margin = unit(c(0, 0, 0, 0), "lines"),
                        complete = TRUE
                )
}

ggplot(crimes, aes(map_id = state, fill = Assault_q)) + 
        geom_map(map = states_map, colour = "black") + 
        scale_fill_manual(values = pal) + 
        expand_limits(x = states_map$long, y = states_map$lat) + 
        coord_map("polyconic") + 
        labs(fill = "Assault Rate\nPercentile") + 
        theme_clean()



## 13.20 基於空間數據格式（shapefile）創建地圖 ggplot() + geom_path() =====

# 這一段都不work！

# install.packages("maptools")
library(maptools)
uk_shp <- readShapePoly("GBR_adm/GBR_adm2.shp")
uk_map <- fortify(uk_shp)
ggplot(uk_ap, aes(x = long, y = lat, group = group)) + geom_path()

str(uk_shp)
uk_map
ggplot(uk_shp, aes(x = long, y = lat, group = group)) + geom_path()



# 14. 輸出圖形 =====

## 14.1 輸出為PDF向量文件 =====

pdf("myplot.pdf", width = 4, height = 4)
plot(mtcars$wt, mtcars$mpg)
print(ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point())
dev.off()

pdf("myplot.pdf", width = 8/2.54, height = 8/2.54)


ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()
ggsave("myplot.pdf", width = 8, height = 8, units = "cm") # ggsave()簡單保存ggplot最後畫的一張圖
# ggsave()無法用於創建多頁圖形



## 14.2 輸出為SVG向量文件 =====

svg("myplot.svg", width = 4, height = 4)
plot(mtcars$wt, mtcars$mpg)
dev.off()



## 14.3 輸出為WMF向量文件 =====

win.metafile("myplot.wmf", width = 4, height = 4)
plot(mtcars$wt, mtcars$mpg)
dev.off()

ggsave("myplot.wmf", width = 8, height = 8, units = "cm")



## 14.4 編輯向量格式的輸出文件 =====

pdf("myplot.pdf", width = 4, height = 4, useDingbats = FALSE)
ggsave("myplot.pdf", width = 4, height = 4, useDingbats = FALSE)
# 避免圖形被其他軟體辨識為字符而無法顯示的問題



## 14.5 輸出為點陣(PNG/TIFF)文件 =====

png("myplot.png", width = 400, height = 400)
plot(mtcars$wt, mtcars$mpg)
dev.off()

png("myplot-%d.png", width = 400, height = 400) # 輸出多幅圖形，加入%d
# width, height此處為像素
plot(mtcars$wt, mtcars$mpg)
print(ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point())
dev.off()

ppi <- 300
png("myplot.png", width = 4 * ppi, height = 4 * ppi, res = ppi)
plot(mtcars$wt, mtcars$mpg)
dev.off()

ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()
ggsave("myplot.png", width = 8, height = 8, unit = "cm", dpi = 300)

# install.packages("Cairo")
# CairoPNG("myplot.png")
# plot(mtcars$wt, mtcars$mpg)
# dev.off()



## 14.6 在PDF文件中使用字體 =====

# install.packages("extrafont") # 可以用於創建包含其它字體的pdf文件
library(extrafont)
font_import()
fonts()

library(extrafont)
loadfonts()

# Sys.seteny(R_GSCMD = "C:/Program Files/gs/gs9.05/bin/gswin32c.exe")

ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point() + 
        ggtitle("Title text goes here") + 
        theme(text = element_text(size = 16, family = "Impact"))
ggsave("myplot.pdf", width = 4, height = 4)
embed_fonts("myplot.pdf")



## 14.7 在Windows的點陣或螢幕輸出中使用字體 =====

# install.packages("extrafont")
library(extrafont)
font_import()
fonts()

library(extrafont)
loadfonts("win")

ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point() + 
        ggtitle("Title text goes here") + 
        theme(text = element_text(size = 16, family = "Georgia", face = "italic"))
ggsave("myplot.png", width = 4, height = 4, dpi = 300)


# 15. reshape =====

library(gcookbook)
heightweight
str(heightweight)



## 15.1 創建 data.frame =====

g <- c("A", "B", "C")
x <- 1:3
dat <- data.frame(g, x)
dat

lst <- list(group = g, value = x)
dat <- as.data.frame(lst)



## 15.2 從數據框中提取訊息 =====

str(ToothGrowth)
tg <- ToothGrowth
tg$supp <- as.character(tg$supp)
str(tg)
ToothGrowth$supp
tg$supp



## 15.3 向數據框添加欄位 =====

data$newcol <- NA
data$newcol <- vector()



## 15.4 從數據框中刪除一欄 =====

data$badcol <- NULL
data <- subset(data, select = -badcol)
data <- subset(data, select = c(-badcol, -othercol))



## 15.5 重命名數據框欄位名稱 =====

names(dat) <- c("name1", "name2", "name3")
library(gcookbook)
names(anthoming)
names(anthoming)[names(anthoming) == "ctrl"] <- c("Control")
names(anthoming)[names(anthoming) == "expt"] <- c("Experimental")
names(anthoming)

names(anthoming)[1] <- "Angle"
names(anthoming)



## 15.6 重排序數據框的欄位 =====

dat <- dat[c(1, 3, 2)]
dat <- dat[c("col1", "col3", "col2")]
library(gcookbook)
anthoming
anthoming[c(1, 3, 2)]
anthoming[, c(1, 3, 2)]
anthoming[3]
anthoming[, 3]
anthoming[, 3, drop = FALSE]
# 注意！！！ df[, , drop = FALSE]
class(anthoming[, 3]) # vector
class(anthoming[, 3, drop = FALSE]) # data.frame



## 15.7 從數據框提取子集 subset(, Source, select = c()) =====

library(gcookbook)
climate
subset(climate, Source == "Berkeley", select = c(Year, Anomaly10y))
subset(climate, Source == "Berkeley" & Year >= 1900 & Year <= 2000, 
       select = c(Year, Anomaly10y))
climate[climate$Source == "Berkeley" & climate$Year >= 1900 & climate$Year <= 2000, 
        c("Year", "Anomaly10y")]
climate[climate$Source == "Berkeley" & climate$Year >= 1900 & climate$Year <= 2000,
        c("Year", "Anomaly10y"), drop = FALSE]
climate[1:100, c(2, 5)]

# 注意！！！盡可能使用名稱做index



## 15.8 改變因子水平的順序 factor(, levels = c()) =====

sizes <- factor(c("small", "large", "large", "small", "medium"))
sizes # 注意！！！default情況因子是照字母排序

# 注意！！！
sizes <- factor(sizes, levels = c("small", "medium", "large"))
sizes

factor(sizes, levels = rev(levels(sizes))) # 顛倒因子順序



## 15.9 根據數據的值改變因子水平的順序 reorder(, , FUN = ) =====

# 注意！！！這超神！！！
iss <- InsectSprays
iss$spray
iss$spray <- reorder(iss$spray, iss$count, FUN = mean)
iss$spray



## 15.10 改變因子水平的名稱 revalue(, c()); mapvalues(, c(), c()) =====

library(plyr)
sizes <- factor(c("small", "large", "large", "small", "medium"))
sizes
levels(sizes)
sizes1 <- revalue(sizes, c(small = "S", medium = "M", large = "L"))
sizes1
revalue(sizes, c("small" = "S", "medium" = "M", "large" = "L")) # 若原因子名稱中有空格等特殊字符會很有用！！！
# 注意！！！
mapvalues(sizes, c("small", "medium", "large"), c("S", "M", "L"))


sizes <- factor(c("small", "large", "large", "small", "medium"))
levels(sizes)[levels(sizes) == "large"] <- "L" # 這種方法可以單改一個factor的名稱
levels(sizes)[levels(sizes) == "medium"] <- "M"
levels(sizes)[levels(sizes) == "small"] <- "S"
sizes

sizes <- factor(c("small", "large", "large", "small", "medium"))
levels(sizes) <- list(S = "small", M = "medium", L = "large")
sizes

sizes <- factor(c("small", "large", "large", "small", "medium"))
levels(sizes)[1] <- "L" # 盡量不要使用位置的方式來改，用名稱來改比較不會出錯
sizes

levels(sizes) <- c("L", "M", "S")
sizes



## 15.11 去掉因子中不再使用的水平 droplevels() =====

sizes <- factor(c("small", "large", "large", "small", "medium"))
sizes <- sizes[1:3] # medium沒用到了
sizes
# 注意！！！droplevels()
sizes <- droplevels(sizes)
sizes



## 15.12 在字符向量中改變元素的名稱 revalue(); mapvalues() =====

sizes <- c("small", "large", "large", "small", "medium")
sizes

sizes1 <- revalue(sizes, c(small = "S", medium = "M", large = "L"))
sizes1

revalue(sizes, c("small" = "S", "medium" = "M", "large" = "L"))
mapvalues(sizes, c("small", "medium", "large"), c("S", "M", "L"))

sizes <- c("small", "large", "large", "small", "medium")
sizes
sizes[sizes == "small"] <- "S"
sizes[sizes == "medium"] <- "M"
sizes[sizes == "large"] <- "L"
sizes



## 15.13 把一個分類變量轉化成另一個分類變量 interaction() =====

pg <- PlantGrowth[c(1, 2, 11, 21, 22), ]
pg
# 注意！！！有點難
oldvals <- c("ctrl", "trt1", "trt2")
newvals <- factor(c("No", "Yes", "Yes"))
pg$treatment <- newvals[match(pg$group, oldvals)]

pg <- PlantGrowth[c(1, 2, 11, 21, 22), ]
pg
pg$treatment[pg$group == "ctrl"] <- "no"
pg$treatment[pg$group == "trt1"] <- "yes"
pg$treatment[pg$group == "trt2"] <- "yes"
class(pg$treatment)
pg$treatment <- factor(pg$treatment)
pg

pg$weightcat[pg$group == "ctrl" & pg$weight < 5] <- "no_small"
pg$weightcat[pg$group == "ctrl" & pg$weight >= 5] <- "no_large"
pg$weightcat[pg$group == "trt1"] <- "yes"
pg$weightcat[pg$group == "trt2"] <- "yes"
pg$weightcat <- factor(pg$newcol)
pg

# 注意！！！
pg$weighttrt <- interaction(pg$weightcat, pg$treatment)
pg



## 15.14 連續變量轉變為分類變量 cut() =====

pg <- PlantGrowth[c(1, 2, 11, 21, 22), ]
pg
pg$wtclass<- cut(pg$weight, breaks = c(0, 5, 6, Inf)) # default為左開右閉，意即不包含最小值
pg
# 注意！！！
cut(pg$weight, breaks = c(0, 5, 6, Inf), right = FALSE)



## 15.15 變量轉換 transform() =====

library(gcookbook)
hw <- heightweight
hw
hw$heightCm <- hw$heightIn * 2.54
hw

hw <- transform(hw, heightCm = heightIn * 2.54, weightKg = weightLb / 2.204)
library(plyr)
# 注意！！！
hw <- mutate(hw, heightCm = heightIn * 2.54, weightKg = weightLb / 2.204)
hw

hw <- transform(hw, bmi = weightKg / (heightCm / 100) ^ 2)
hw <- mutate(hw, bmi = weightKg / (heightCm / 100) ^ 2)
hw$bmi <- hw$weightKg / (hw$heightCm / 100) ^ 2
hw



## 15.16 按組轉換數據 ddply() =====

library(MASS)
library(plyr)
# 注意！！！
cb <- ddply(cabbages, "Cult", transform, DevWt = HeadWt - mean(HeadWt))
cabbages
transform(cabbages, DevWt = HeadWt - mean(HeadWt))

library(plyr)
cb <- ddply(cabbages, "Cult", transform, DevWt = HeadWt - mean(HeadWt)) # 對各組做標準化
cb


ggplot(cb, aes(x = Cult, y = HeadWt)) + geom_boxplot()
ggplot(cb, aes(x = Cult, y = DevWt)) + geom_boxplot()

# 注意！！！可以一次多組
ddply(cabbages, c("Cult", "Date"), transform, 
      DevWt = HeadWt - mean(HeadWt), 
      DevVitC = VitC - mean(VitC))



## 15.17 分組匯總數據 ddply() =====

library(MASS)
library(plyr)
ddply(cabbages, c("Cult", "Date"), summarise, Weight = mean(HeadWt), VitC = mean(VitC))
cabbages

library(plyr)
summarise(cabbages, Weight = mean(HeadWt))

ddply(cabbages, "Cult", summarise, Weight = mean(HeadWt))
ddply(cabbages, c("Cult", "Date"), summarise, Weight = mean(HeadWt),
      sd = sd(HeadWt),
      n = length(HeadWt))

# 處理遺失值
c1 <- cabbages
c1$HeadWt[c(1, 20, 45)] <- NA
ddply(c1, c("Cult", "Date"), summarise, Weight = mean(HeadWt),
      sd = sd(HeadWt), 
      n = length(HeadWt))
ddply(c1, c("Cult", "Date"), summarise, Weight = mean(HeadWt, na.rm = TRUE),
            sd = sd(HeadWt, na.rm = TRUE),
            n = sum(!is.na(HeadWt)))

c2 <- subset(c1, !(Cult == "c52" & Date == "d21"))
c2a <- ddply(c2, c("Cult", "Date"), summarise, 
             Weight = mean(HeadWt, na.rm = TRUE), 
             sd = sd(HeadWt, na.rm = TRUE), 
             n = sum(!is.na(HeadWt)))
c2a

ggplot(c2a, aes(x = Date, fill = Cult, y = Weight)) + 
        geom_bar(position = "dodge", stat = "identity")

# 注意！！！.drop讓bar不會延展至無資料組別
c2b <- ddply(c2, c("Cult", "Date"), .drop = FALSE, summarise,
             Weight = mean(HeadWt, na.rm = TRUE),
             sd = sd(HeadWt, na.rm = TRUE),
             n = sum(!is.na(HeadWt)))
c2b
ggplot(c2b, aes(x = Date, fill = Cult, y = Weight)) + geom_bar(position = "dodge", stat = "identity")



## 15.18 使用標準誤差和信賴區間來匯總數據 =====

library(MASS)
library(plyr)
ca <- ddply(cabbages, c("Cult", "Date"), summarise, 
            Weight = mean(HeadWt, na.rm = TRUE),
            sd = sd(HeadWt, na.rm = TRUE),
            n = sum(!is.na(HeadWt)),
            se = sd / sqrt(n))
ca

ddply(cabbages, c("Cult", "Date"), summarise,
      Weight = mean(HeadWt, na.rm = TRUE),
      sd = sd(HeadWt, na.rm = TRUE),
      n = sum(!is.na(HeadWt)),
      se = sd / sqrt(n))

ciMult <- qt(.975, ca$n - 1)
ciMult
ca$ci <- ca$se * ciMult

ca$ci95 <- ca$se * qt(.975, ca$n)

summarySE <- function(data = NULL, measurevar, groupvars = NULL,
                      conf.interval = .95, na.rm = FALSE, .drop = TRUE) {
        require(plyr)
        length2 <- function(x, na.rm = FALSE) {
                if (na.rm) sum(!is.na(x))
                else length(x)
        }
        
        # 匯總
        datac <- ddply(data, groupvars, .drop = .drop,
                       .fun = function(xx, col, na.rm) {
                               c(n = length2(xx[, col], na.rm = na.rm),
                                 mean = mean(xx[, col], na.rm = na.rm),
                                 sd = sd(xx[, col], na.rm = na.rm))
                       }, 
                       measurevar,
                       na.rm)
        # 注意！！！rename()
        datac <- rename(datac, c("mean" = measurevar))
        datac$se <- datac$sd / sqrt(datac$n)
        ciMult <- qt(conf.interval / 2 + .5, datac$n-1)
        datac$ci <- datac$se * ciMult
        return(datac)
}

c2 <- subset(cabbages, !(Cult == "c52" & Date == "d21"))
c2$HeadWt[c(1, 20, 45)] <- NA
summarySE(c2, "HeadWt", c("Cult", "Date"), conf.interval = .99,
          na.rm = TRUE, .drop = FALSE)



## 15.19 把數據框從寬變長 melt(, id.vars, measure.vars, variable.name, value.name) =====

library(gcookbook)
anthoming
library(reshape2)
melt(anthoming, id.vars = "angle", variable.name = "condition", value.name = "count")

drunk
melt(drunk, id.vars = "sex", measure.vars = c("0-29", "30-39"), variable.name = "name", value.name = "count")

plum_wide
# 注意！！！measure.vars的default設定為id.vars「以外的所有」變量
melt(plum_wide, id.vars = c("length", "time"), variable.name = "survival", value.name = "count")

co <- corneas
co # 注意！！！沒有id.vars的資料集可以自己加入id.vars

co$id <- 1:nrow(co)
melt(co, id.vars = "id", variable.name = "eye", value.name = "thickness")
# 注意！！！用數值作為id.vars可能會給後續分析造成問題，建議轉成character或factor



## 15.20 把數據框從長變寬 dcast() =====

library(gcookbook)
plum
library(reshape2)
dcast(plum, length + time ~ survival, value.var = "count")



## 15.21 把時間序列數據對象拆分成時間和數據 =====

nhtemp
class(nhtemp) # ts
as.numeric(time(nhtemp))
as.numeric(nhtemp)
nht <- data.frame(year = as.numeric(time(nhtemp)), temp = as.numeric(nhtemp))
nht

presidents
pres_rating <- data.frame(
        year = as.numeric(time(presidents)), 
        rating = as.numeric(presidents)
)
# 注意！！！time(), cycle()
pres_rating2 <- data.frame(
        year = as.numeric(floor(time(presidents))),
        quarter = as.numeric(cycle(presidents)),
        rating = as.numeric(presidents)
)
time()



# R資料採礦與數據分析 =====

# install.packages("RGtk2")
# install.packages("rattle")
# install.packages("iClick")
library(iClick)
library(rattle)



# 實用R程式設計 =====
# R繪圖 =====
# 探索資料圖形 =====

data(iris)
str(iris)
summary(iris)
head(iris, 5)
tail(iris, 5)
iris[iris$Species == "setosa", 1]
subset(iris, Species == "setosa", select = Sepal.Length)
iris[iris$Species == "setosa", 1:2]
subset(iris, Species == "setosa", select = c(Sepal.Length, Sepal.Width))

x <- iris[, 1]
c(min(x), max(x))
range(x)
summary(x)
fivenum(x)
mean(x)
median(x)
IQR(x)
mean(x, trim = .1)

sd(x)
cv <- sd(x) / mean(x)

old.par <- par(mfrow = c(1, 2), mex = .5, mar = c(4, 4, 3, 2) + .1)
boxplot(x)
rug(x, side = 4)
boxplot(x, horizontal = TRUE)
rug(x, side = 1)
par(old.par)

old.par <- par(mfrow = c(1, 2), mex = .8, mar = c(5, 5, 3, 1) + .1)
hist(x, freq = TRUE, breaks = "Sturges")
rug(x, side = 1)
hist(x, prob = TRUE, breaks = "Sturges", col = "lightblue", border = "magenta")
rug(x, side = 1)
par(old.par)

hist(x, prob = TRUE, breaks = "Sturges")
hist(x, prob = TRUE, breaks = "Scott")
hist(x, prob = TRUE, breaks = "Freedman-Diaconis")

hist(x, prob = TRUE, breaks = 10)
hist(x, prbo = TRUE, breaks = seq(from = 4, to = 8, by = .25))

library(MASS)
old.par <- par(mfrow = c(1, 2), mex = .8, mar = c(5, 5, 3, 1) + .1)
truehist(x, prob = FALSE, ylab = "Frequency", main = "Histogram")
truehist(x, prob = TRUE, ylab = "Density", main = "Histogram")
par(old.par)

stem(x, scale = .5)
sum(x == 4.4)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
stripchart(x, method = "overplot", at = .7)
text(6, .65, "overplot")
stripchart(x, method = "stack", add = TRUE, at = .85)
text(6, .8, "stack")
stripchart(x, method = "jitter", add = TRUE, at = 1.2)
text(6, 1.05, "jitter")
title(main = "strip chart")
par(old.par)

y <- cut(x, breaks = 6)
z <- table(y)
old.par <- par(mfrow = c(2, 2), mex = .2, mar = c(3, 3, 3, 2) + .1)
pie(z)
pie(z, clockwise = TRUE)
pie(z, col = terrain.colors(6))
pie(z, col = gray(seq(from = .4, to = 1.0, length = 6)))
par(old.par)

old.par <- par(mfrow = c(1, 2), mex = .8, mar = c(5, 4, 3, 1) + .1)
plot(density(x), col = "red", main = "Kernel density estimate")
rug(x, side = 1)
hist(x, prob = TRUE, breaks = "Sturges", main = "Histogram and KDE")
lines(density(x), col = "red")
rug(x, side = 1)
par(old.par)

f <- density(x)
class(f)
names(f)
print(f)
plot(f, type = "n")
polygon(f, col = "wheat")

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
plot.ecdf(x) # 累積機率密度圖(emperical cumulative density function)
par(old.par)

F <- ecdf(x)
class(F)
names(F)
print(F)
summary(F)
plot(F)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
qqnorm(x)
qqline(x, col = "red", lwd = 2)
par(old.par)

data(VADeaths)
str(VADeaths)
class(VADeaths)
VADeaths
names <- c("RM", "RF", "UM", "UF")
colors <- c("lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk")
old.par <- par(mfrow = c(2, 2), mex = .8, mar = c(3, 3, 3, 2) + .1)
barplot(VADeaths, names.arg = names)
barplot(VADeaths, names.arg = names, horiz = TRUE)
barplot(VADeaths, names.arg = names, col = colors, border = "blue")
barplot(VADeaths, names.arg = names, col = colors, border = "blue", space = 1.5)
par(old.par)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
barplot(VADeaths, beside = TRUE, col = c("lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk"), legend.text = rownames(VADeaths), ylim = c(0, 100))
title(main = "Death Rates in Virginia", font.main = 4)
par(old.par)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
barplot(t(VADeaths), beside = TRUE, col = c("lightblue", "mistyrose", "lightcyan", "lavender"), legend.text = rownames(t(VADeaths)), ylim = c(0, 80), args.legend = list(x = "topleft"))
title(main = "Death Rates in Virginia", font.main = 4)
par(old.par)

colnames(VADeaths) <- c("RM", "RF", "UM", "UF")
old.par <- par(mfrow = c(1, 2), mex = .8, mar = c(5, 4, 3, 1) + .1)
dotchart(VADeaths, xlim = c(0, 100), xlab = "Deaths per 1000", main = "Death rates")
dotchart(t(VADeaths), xlim = c(0, 100), xlab = "Deaths per 1000", main = "Death Rates")
par(old.par)

data(warpbreaks)
str(warpbreaks)
with(warpbreaks, tapply(breaks, INDEX = wool, FUN = sum))
with(warpbreaks, tapply(breaks, INDEX = tension, FUN = sum))
with(warpbreaks, tapply(breaks, INDEX = list(wool, tension), FUN = sum))
xtabs(breaks ~ wool, data = warpbreaks)
xtabs(breaks ~ tension, data = warpbreaks)
xtabs(breaks ~ wool + tension, data = warpbreaks)
ftable(xtabs(breaks ~ wool + tension, data = warpbreaks))

t1 <- with(warpbreaks, tapply(breaks, INDEX = list(wool, tension), FUN = sum))
t2 <- with(warpbreaks, tapply(breaks, INDEX = list(wool, tension), FUN = mean))
old.par <- par(mfrow = c(1, 2), mex = .8, mar = c(5, 4, 3, 1) + .1)
barplot(t1, beside = TRUE, col = c("lightblue", "mistyrose"), main = "counts", legend.text = rownames(t1), ylim = c(0, max(t1)))
barplot(t2, beside = TRUE, col = c("lightblue", "mistyrose"), main = "means", legend.text = rownames(t2), ylim = c(0, max(t2)))
par(old.par)

brks <- as.integer(xtabs(breaks ~ wool + tension, data = warpbreaks))
label <- c("AL", "BL", "AM", "BM", "AH", "BH")
old.par <- par(mfrow = c(1, 2), mex = .2, mar = c(3, 3, 3, 2) + .1)
pie(brks, label = label)
pie(brks, label = label, col = gray(seq(from = 0.4, to = 1.0, length = 8)))
par(old.par)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
pairs(iris[, 1:4], panel = panel.smooth)
par(old.par)

panel.hist <- function(x, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5))
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$counts
        y <- y / max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
        lines(density(x, na.rm = TRUE), col = "red")
}
old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
pairs(iris[, 1:4], panel = panel.smooth, pch = 1, bg = "lightcyan", diag.panel = panel.hist, font.labels = 2, cex.labels = 1.2)
par(old.par)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
pairs(iris[, 1:4], pch = c(1, 2, 4)[iris$Species], col = c("red", "green", "blue")[iris$Species])
par(old.par)

library(lattice)
old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
xyplot(Sepal.Length + Sepal.Width ~ Petal.Length + Petal.Width | Species, data = iris, layout = c(2, 2), scales = list(cex = .5, cex.lab = .5, relation = "free"), auto.key = TRUE)

setosa <- iris[iris$Species == "setosa", 1:4]
old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
boxplot(setosa, names = c("sep.len", "sep.wid", "pet.len", "pet.wid"), main = "Iris setosa")
par(old.par)

old.par <- par(mfrow = c(1, 2), mex = .5, mar = c(5, 4, 4, 2) + .1)
with(iris, boxplot(Sepal.Length ~ Species, main = "Sepal length"))
with(iris, boxplot(Sepal.Length ~ Species, notch = TRUE, main = "Sepal length"))
par(old.par)

sx <- with(iris, split(Sepal.Length, Species))
old.par <- par(mfrow = c(1, 2), mex = .5, mar = c(5, 4, 4, 2), + .1)
boxplot(sx, main = "Sepal length")
boxplot(sx, notch = TRUE, main = "Sepal length")
par(old.par)

sx <- with(iris, split(Sepal.Length, Species))
sy <- with(iris, split(Sepal.Width, Species))

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
plot(0, xlim = range(sx), ylim = range(sy), type = "n", xlab = "x", ylab = "y")
points(sx[[1]], sy[[1]], pch = 1, col = 1)
points(sx[[2]], sy[[2]], pch = 2, col = 2)
points(sx[[3]], sy[[3]], pch = 3, col = 3)
for(i in 1:3) abline(lm(sy[[i]] ~ sx[[i]]), col = i)
legend("topright", legend = c("setosa", "versicolor", "virginica"), lty = 1, pch = 1:3, col = 1:3)
par(old.par)

x <- iris[[1]]
y <- iris[[2]]
species <- iris[[5]]
library(lattice)
old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
xyplot(y ~ x, groups = species, type = c("g", "p", "r"), auto.key = TRUE)
par(old.par)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
xyplot(y ~ x | species, type = c("g", "p", "r"), auto.key = TRUE)
par(old.par)

x <- iris[[1]]
y <- iris[[2]]

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
qqplot(x, y)
abline(a = 0, b = 1, col = "red")
par(old.par)

f <- function(x, y) {
        r <- sqrt(x^2 + y^2)
        10 * sin(r) / r
}
x <- seq(from = -10, to = 10, length = 50)
y <- x
z <- outer(x, y, f)
z[is.na(z)] <- 1

old.par <- par(mfrow = c(1, 2), mex = .4, mar = c(5, 4, 3, 1) + .1, bg = "white")
persp(x, y, z, theta = 30, phi = 30, expand = .5, col = "lightgreen")
persp(x, y, z, theta = 30, phi = 30, expand = .5, col = "lightblue", ltheta = 120, shade = .75, ticktype = "detailed", xlab = "x", ylab = "y", zlab = "f(x, y)")
par(old.par)

persp(x, y, z, theta = 30, phi = 30, expand = .5, col = "lightgreen", box = FALSE)

library(lattice)
f <- function(x, y) {
        r <- sqrt(x^2 + y^2);
        10 * sin(r) / r
}
y <- x <- seq(from = -10, to = 10, length = 50)
h <- expand.grid(x = x, y = y)
z <- f(h$x, h$y)
z[is.na(z)] <- 1
old.par <- par(mex = .8, mar = c(5, 5, 3, 1) + .1)
wireframe(z ~ h$x * h$y, xlab = "x", ylab = "y", zlab = "f(x, y)")
par(old.par)

old.par <- par(mex = .8, mar = c(5, 5, 3, 1) + .1, cex.axis = .5)
wireframe(z ~ h$x * h$y, xlab = "x", ylab = "y", zlab = "f(x, y)", scales = list(arrows = FALSE), light.source = c(10, 0, 10), col = "lightblue")
par(old.par)

wireframe(z ~ h$x * h$y, xlab = "x", ylab = "y", zlab = "f(x, y)", scales = list(arrows = TRUE), shade = TRUE, light.source = c(10, 0, 10), col = "lightblue")
wireframe(z ~ h$x * h$y, xlab = "x", ylab = "y", zlab = "f(x, y)", scales = list(arrows = FALSE), shade = TRUE, light.source = c(10, 0, 10), col = "lightblue", colorkey = TRUE)
wireframe(z ~ h$x * h$y, xlab = "x", ylab = "y", zlab = "f(x, y)", scales = list(arrows = FALSE), shade = TRUE, light.source = c(10, 0, 10), col = "lightblue", colorkey = TRUE, drape = TRUE, screen = list(z = 30, x = -60))

library(rgl)

f <- function(x, y) {
        r <- sqrt(x^2 + y^2);
        10 * sin(r) / r
}
y <- x <- seq(from = -10, to = 10, length = 50)
z <- outer(x, y, f)
z[is.na(z)] <- 1

persp3d(x, y, z, xlab = "x", ylab = "y", zlab = "f(x, y)", col = "lightblue")

rgl.snapshot("/Users/thomas/Documents/R Files/persp1.png")
snapshot3d("/Users/thomas/Documents/R Files/persp2.png")

persp3d(x, y, z, xlab = "x", ylab = "y", zlab = "f(x, y)", col = "lightblue")
play3d(spin3d(axis = c(0, 0, 1), rpm = 8), duration = 10)
play3d(spin3d(axis = c(0, 1, 0), rpm = 6), duration = 10)
play3d(spin3d(axis = c(0, 1, 1), rpm = 8), duration = 5)

f <- function(x, y) {
        r <- sqrt(x^2 + y^2);
        10 * sin(r) / r
}
x <- seq(from = -10, to = 10, length = 100)
y <- x
z <- outer(x, y, f)
z[is.na(z)] <- 1

old.par <- par(mfrow = c(1, 2), mex = .5, mar = c(5, 4, 3, 1) + .1)
image(x, y, z, col = cm.colors(20))
contour(x, y, z, add = TRUE)
col <- gray(seq(from = .4, to = 1.0, length = 6))
image(x, y, z, col = col)
contour(x, y, z, add = TRUE)
par(old.par)

old.par <- par(mfrow = c(1, 2), mex = .5, mar = c(6, 4, 3, 1) + .1)
image(x, y, z, col = cm.colors(20))
contour(x, y, z, levels = c(-5, -2, 0, 2, 5), lwd = c(1, 1, 2, 1, 1), add = TRUE)
image(x, y, z, col = cm.colors(20))
contour(x, y, z, levels = 0, lwd = 2, add = TRUE)
par(old.par)

old.par <- par(mex = .8, mar = c(3, 3, 2, 1) + .1)
filled.contour(x, y, z)
par(old.par)

old.par <- par(mex = .8, mar = c(3, 3, 2, 1) + .1)
filled.contour(x, y, z, levels = c(-5, -2, 0, 2, 5), lwd = c(1, 1, 2, 1, 1))
par(old.par)

data(iris)
x <- iris[, 1]
y <- iris[, 2]
z <- iris[, 3]

library(lattice)
old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
cloud(z ~ x * y, groups = iris$Species, pch = 1:3, col = 1:3, scales = list(arrows = FALSE), light.source = c(10, 0, 10))
par(old.par)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
cloud(z ~ x * y, groups = iris$Species, pch = 1:3, col = 1:3, scales = list(arrows = FALSE), light.source = c(10, 0, 10), shade = TRUE, colorkey = TRUE)
par(old.par)

library(scatterplot3d)
old.par <- par(mfrow = c(1, 2), mex = .5, mar = c(5, 5, 3, 1) + .1)
scatterplot3d(x, y, z, xlab = "x", ylab = "y", zlab = "z", angle = 30, y.margin.add = .1, scale.y = .7, pch =c(1, 2, 3)[iris$Species], color = c("red", "green", "blue")[iris$Species])
scatterplot3d(x, y, z, xlab = "x", ylab = "y", zlab = "z", type = "h", angle = 30, y.margin.add = .1, scale.y = .7, pch = c(1, 2, 3)[iris$Species], color = c("red", "green", "blue")[iris$Species])
par(old.par)

library(rgl)
plot3d(x, y, z, pch = c(1, 2, 3)[iris$Species], col = c("red", "green", "blue")[iris$Species])
plot3d(x, y, z, pch = c(1, 2, 3)[iris$Species], col = c("red", "green", "blue")[iris$Species])
play3d(spin3d(axis = c(0, 0, 1), rpm = 8), duration = 10)
play3d(spin3d(axis = c(0, 1, 0), rpm = 8), duration = 10)
play3d(spin3d(axis = c(0, 1, 1), rpm = 8), duration = 10)

# install.packages("SemiPar")
library(SemiPar)
data(pig.weights, package = "SemiPar")
str(pig.weights)
xtabs(weight ~ id.num + num.weeks, data = pig.weights)
ftable(xtabs(weight ~ id.num + num.weeks, data = pig.weights))

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
plot(weight ~ num.weeks, data = pig.weights)
for(i in 1:48) lines(weight ~ num.weeks, subset(pig.weights, id.num == i), col = i)
par(old.par)

library(lattice)
xyplot(weight ~ num.weeks, data = pig.weights, groups = id.num, type = "b")

library(nlme)
pig.growth <- groupedData(weight ~ num.weeks | id.num, data = pig.weights)
class(pig.growth)
plot(pig.growth, outer = ~1, key = FALSE)
plot(pig.growth, outer = ~1, key = TRUE)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
xyplot(weight ~ num.weeks | id.num, data = pig.weights, type = c("b", "g", "p", "r"), scales = list(cex = .5, cex.axis = .5, cex.lab = .5), auto.key = TRUE)
par(old.par)

library(nlme)
pig.growth <- groupedData(weight ~ num.weeks | id.num, data = pig.weights)
class(pig.growth)
plot(pig.growth, pch = 16)



# R bloggers =====

# ggExtra =====
# An awesome RStudio addin for selecting colours, and another for adding marginal density plots to ggplot2
# https://www.r-bloggers.com/an-awesome-rstudio-addin-for-selecting-colours-and-another-for-adding-marginal-density-plots-to-ggplot2/
# install.packages("ggExtra")
library(ggExtra)
myplot <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
mynewplot <- ggMarginalGadget(myplot)


# ggrepel =====
# Avoid overlapping labels in ggplot2 charts
# https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
# install.packages("ggrepel")
library(ggrepel)
set.seed(42)
ggplot(mtcars) +
        geom_point(aes(wt, mpg), color = 'red') +
        geom_text_repel(aes(wt, mpg, label = rownames(mtcars))) +
        theme_classic(base_size = 16)

set.seed(42)
ggplot(mtcars) +
        geom_point(aes(wt, mpg), size = 5, color = 'grey') +
        geom_label_repel(
                aes(wt, mpg, fill = factor(cyl), label = rownames(mtcars)),
                fontface = 'bold', color = 'white',
                box.padding = unit(0.35, "lines"),
                point.padding = unit(0.5, "lines"),
                segment.color = 'grey50'
        ) +
        theme_classic(base_size = 16)


# Fishbone Diagram 魚骨圖; Cause-and-Effect Diagram 因果圖; Ishikawa Diagram 石川圖 =====

# install.packages("qcc")
library(qcc)
cause.and.effect(
        cause = list(
                INGREDIENTS = c("Coffee", "Sugar", "Water"), 
                PEOPLE = c("You"),
                ENERGY = c("Gas"), 
                TOOLS = c("Spoon", "Cooker", "Cup")
                ), 
        effect = "A cup of coffee", 
        title = "Fishbone chart for a cup of coffee"
)



# Venn Diagram 文氏圖 =====

# https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html
# install.packages("VennDiagram")
library(VennDiagram)
grid.newpage()
draw.pairwise.venn(area1 = 22, area2 = 9, cross.area = 3, category = c("Dog People", 
                                                                         "Cat People"))
?draw.pairwise.venn



# DataCamp ==================================================

## Data Visualization with ggplot2 (Part 3) =====

## Refresher (1) -----
# Create movies_small
library(ggplot2movies)
set.seed(123)
movies_small <- movies[sample(nrow(movies), 1000), ]
movies_small$rating <- factor(round(movies_small$rating))
# Explore movies_small with str()
str(movies_small)
# Build a scatter plot with mean and 95% CI
ggplot(movies_small, aes(x = rating, y = votes)) +
        geom_point() +
        stat_summary(fun.data = "mean_cl_normal",
                     geom = "crossbar",
                     width = 0.2,
                     col = "red") + 
        scale_y_log10()

## Refresher (2) -----
# Reproduce the plot
ggplot(diamonds, aes(x = carat, y = price, col = color)) +
        geom_point(alpha = 0.5, size = 0.5, shape = 16) +
        scale_x_log10(expression(log[10](Carat)), limits = c(0.1, 10)) +
        scale_y_log10(expression(log[10](Price)), limits = c(100, 100000)) +
        scale_color_brewer(palette = "YlOrRd") +
        coord_equal() +
        theme_classic()


## Transformations! -----
# movies_small is available
# Add a boxplot geom
d <- ggplot(movies_small, aes(x = rating, y = votes)) +
        geom_point() +
        geom_boxplot() +
        stat_summary(fun.data = "mean_cl_normal",
                     geom = "crossbar",
                     width = 0.2,
                     col = "red")
# Untransformed plot
d
# Transform the scale
d + scale_y_log10()
# Transform the coordinates
d + coord_trans(y = "log10")





# 輕鬆學習R語言 =============================================================

# 探索資料分析 =====

nrow(iris); ncol(iris)
dim(iris)
head(iris); tail(iris)
summary(iris)
str(iris) # structure

hist(rnorm(1000))
boxplot(Sepal.Length ~ Species, data = iris)


x <- seq(from = as.Date("2017-01-01"), to = as.Date("2017-01-31"), by = 1)
set.seed(123)
y <- sample(1:100, size = 31, replace = TRUE)
plot(x, y, type = "l", xaxt = "n")
axis.Date(1, at = x, format = "%Y-%m-%d")


class(AirPassengers)
class(LakeHuron)
plot(AirPassengers)
plot(LakeHuron)


plot(cars$speed, cars$dist)

plot(iris)


ice_cream_flavor <- rep(NA, times = 100)
for (i in 1:100) {
        ice_cream_flavor[i] <- sample(c("vanilla", "chocolate", "matcha", "other"), size = 1)
}
ice_cream_flavor
table(ice_cream_flavor)
barplot(table(ice_cream_flavor))


curve(sin, from = -pi, to = pi)

my_sqr <- function(x) {
        return(x ^ 2)
}
curve(my_sqr, from = -3, to = 3)


plot(cars, main = "Car speed vs. breaking distance", xlab = "Car speed(mph)", ylab = "Breaking distance(ft)")
grid()


barplot(table(ice_cream_flavor), horiz = TRUE, las = 1) # 0 ~ 3
barplot(table(ice_cream_flavor), horiz = TRUE, las = 1, cex.name = 0.8, cex.axis = 1.2)
# cex: character expansion factor

norm_dist <- rnorm(1000)
hist(norm_dist, freq = FALSE)
lines(density(norm_dist))


plot(cars, pch = 2, col = "red") # pch: plotting character


iris_pch <- c(1, 2, 3)[as.numeric(iris$Species)]
plot(iris$Sepal.Length, iris$Sepal.Width, col = iris$Species, pch = iris_pch)


par(mfrow = c(2, 2)) # matrix of figures entered row-wise
boxplot(iris$Sepal.Length ~ iris$Species, main = "Sepal length by species")
boxplot(iris$Sepal.Width ~ iris$Species, main = "Sepal width by species")
boxplot(iris$Petal.length ~ iris$Species, main = "Petal length by species")
boxplot(iris$Petal.Width ~ iris$Species, main = "Petal width by species")



# 探索資料分析2 =====

library(ggplot2) # gg: grammer of graphics
ggplot(cars, aes(x = speed, y = dist))

ggplot(cars, aes(x = speed, y = dist)) + geom_point()

set.seed(123)
norm_nums <- rnorm(1000)
hist_df <- data.frame(norm_nums = norm_nums)
ggplot(hist_df, aes(x = norm_nums)) + geom_histogram()

ggplot(hist_df, aes(x = norm_nums)) + geom_histogram(binwidth = 0.1)
ggplot(hist_df, aes(x = norm_nums)) + geom_histogram(binwidth = 0.5)

ggplot(iris, aes(x = Species, y = Sepal.Length)) + geom_boxplot()


x <- seq(from = as.Date("2017-01-01"), to = as.Date("2017-01-31"), by = 1)
set.seed(123)
y <- sample(1:100, size = 31, replace = TRUE)
line_df <- data.frame(x = x, y = y)
ggplot(line_df, aes(x = x, y = y)) + geom_line() + scale_x_date(date_labels = "%m.%d")


ggplot(cars, aes(x = speed, y = dist)) + geom_point()


ice_cream_flavor <- rep(NA, times = 100)
for(i in 1:100) {
        ice_cream_flavor[i] <- sample(c("vanilla", "chocolate", "matcha", "other"), size = 1)
}
ice_cream_df <- data.frame(ice_cream_flavor = ice_cream_flavor)
ggplot(ice_cream_df, aes(x = ice_cream_flavor)) + geom_bar() # 如果資料未經彙整


flavor <- names(table(ice_cream_flavor))
votes <- as.vector(unname(table(ice_cream_flavor))) # unname()
ice_cream_df <- data.frame(flavor = flavor, votes = votes)
ice_cream_df
ggplot(ice_cream_df, aes(x = flavor, y = votes)) + geom_bar(stat = "identity")
# 因為ggplot2會協助計算類別個數，所以資料框本身已經算好，就要加stat = "identity"


sin_df <- data.frame(x = c(-pi, pi))
sin_df
ggplot(sin_df, aes(x = x)) + stat_function(fun = sin, geom = "line")

my_sqr <- function(x) {
        return(x ^ 2)
}
my_sqr_df <- data.frame(x = c(-3, 3))
ggplot(my_sqr_df, aes(x = x)) + stat_function(fun = my_sqr, geom = "line")


ggplot(cars, aes(x = speed, y = dist)) + 
        geom_point() + 
        ggtitle("Car speed vs. braking distance") + 
        xlab("Car speed(mph)") + 
        ylab("Breaking distance(ft)")


ggplot(cars, aes(x = speed, y = dist)) + 
        geom_point() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())


ice_cream_flavor <- rep(NA, times = 100)
for(i in 1:100) {
        ice_cream_flavor[i] <- sample(c("vanilla", "chocolate", "matcha", "other"), size = 1)
}
ice_cream_df <- data.frame(ice_cream_flavor = ice_cream_flavor)
ggplot(ice_cream_df, aes(x = ice_cream_flavor)) + 
        geom_bar() + 
        coord_flip()


set.seed(123)
norm_nums <- rnorm(100)
hist_df <- data.frame(norm_nums = norm_nums)
ggplot(hist_df, aes(x = norm_nums)) + 
        geom_histogram(binwidth = 0.5, aes(y = ..density..), alpha = 0.5) + 
        geom_density()

ggplot(hist_df, aes(x = norm_nums, y = ..density..)) + 
        geom_histogram(binwidth = 0.5, alpha = 0.5) + 
        geom_density()


ggplot(cars, aes(x = speed, y = dist)) + 
        geom_point(shape = 2, colour = "red")


ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width)) + 
        geom_point(aes(shape = Species, colour = Species))


library(ggplot2)
library(gridExtra)
g1 <- ggplot(iris, aes(x = Species, y = Sepal.Length)) + 
        geom_boxplot()
g2 <- ggplot(iris, aes(x = Species, y = Sepal.Width)) + 
        geom_boxplot()
g3 <- ggplot(iris, aes(x = Species, y = Petal.Length)) + 
        geom_boxplot()
g4 <- ggplot(iris, aes(x = Species, y = Petal.Width)) + 
        geom_boxplot()
grid.arrange(g1, g2, g3, g4, nrow = 2, ncol = 2)





### R軟體資料分析基礎與應用 =====

# 07 統計繪圖 =====

require(ggplot2)
data(diamonds)
head(diamonds)
hist(diamonds$carat, main = "Carat Histogram", xlab = "Carat")
plot(price ~ carat, data = diamonds)
plot(diamonds$carat, diamonds$price)
identical(plot(price ~ carat, data = diamonds), plot(diamonds$carat, diamonds$price))
boxplot(diamonds$carat)

ggplot(data = diamonds) + geom_histogram(aes(x = carat))
ggplot(data = diamonds) + geom_point(aes(x = carat, y = price))
ggplot(data = diamonds, aes(x = carat, y = price)) + geom_point()
g <- ggplot(diamonds, aes(x = carat, y = price))
g + geom_point(aes(color = color))
g + geom_point(aes(color = color)) + facet_wrap(~ color)
g + geom_point(aes(color = color)) + facet_grid(cut ~ clarity)
ggplot(diamonds, aes(x = carat)) + geom_histogram() + facet_wrap(~color)
ggplot(diamonds, aes(y = carat, x = 1)) + geom_boxplot()
ggplot(diamonds) + geom_boxplot(aes(y = carat, x = clarity))
ggplot(diamonds, aes(x = cut, y = carat)) + geom_boxplot()
ggplot(diamonds, aes(x = cut, y = carat)) + geom_violin()
ggplot(diamonds, aes(x = cut, y = carat)) + geom_violin() + geom_point()

ggplot(economics, aes(x = date, y = pop)) + geom_line()
head(economics)

require(lubridate)
economics$year <- year(economics$date)
economics$month <- month(economics$date)
head(economics)

econ2000 <- economics[which(economics$year >= 2000), ]

require(scales)
g <- ggplot(econ2000, aes(x = month, y = pop))
g <- g + geom_line(aes(color = factor(year), group = year))
g <- g + scale_color_discrete(name = "Year")
g <- g + scale_y_continuous(labels = comma)
g <- g + labs(title = "Population Growth", x = "Month", y = "Population")
g

g2 <- ggplot(diamonds, aes(x = carat, y = price)) + geom_point(aes(color = color))
g2
# install.packages("ggthemes")
require(ggthemes)
g2 + theme_economist()
g2 + theme_economist() + scale_colour_economist()
g2 + theme_excel()
g2 + theme_excel() + scale_colour_excel()
g2 + theme_tufte()
g2 + theme_wsj()


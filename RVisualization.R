# 1. R 基礎 =====


# 2. 快速探索數據 =====


## 2.1 散佈圖 =====

library(ggplot2)
qplot(mtcars$wt, mtcars$mpg)
qplot(wt, mpg, data = mtcars)
ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()


## 2.2 折線圖 =====

plot(pressure$temperature, pressure$pressure, type = "l")

plot(pressure$temperature, pressure$pressure, type = "l")
points(pressure$temperature, pressure$pressure) # 記得是point"s"
lines(pressure$temperature, pressure$pressure/2, col = "red")
points(pressure$temperature, pressure$pressure/2, col = "red")

library(ggplot2)
qplot(pressure$temperature, pressure$pressure, geom = "line") # geom = 
qplot(temperature, pressure, data = pressure, geom = "line") # 方法2
ggplot(pressure, aes(x = temperature, y = pressure)) + geom_line() #方法3

qplot(temperature, pressure, data = pressure, geom = c("line", "point"))


# 2.3 長條圖 =====

barplot(BOD$demand, names.arg = BOD$Time)

barplot(table(mtcars$cyl))

library(ggplot2)
# qplot(BOD$Time, BOD$demand, geom = "bar", stat = "identity") # Error
# qplot(factor(BOD$Time), BOD$demand, geom = "bar", stat = "identity") # Error

# qplot(Time, demand, data = BOD, geom = "bar", stat = "identity") # Error
ggplot(BOD, aes(x = Time, y = demand)) + geom_bar(stat = "identity")

qplot(factor(cyl), data = mtcars)
ggplot(mtcars, aes(x = factor(cyl))) + geom_bar()


## 2.4 直方圖 =====

hist(mtcars$mpg)
qplot(mtcars$mpg)

qplot(mpg, data = mtcars, binwidth = 4) # binwidth =


## 2.5 箱型圖 =====

plot(ToothGrowth$supp, ToothGrowth$len) # X為factor，所以自動畫出boxplot
qplot(ToothGrowth$supp, ToothGrowth$len, geom = "boxplot") # another way
ggplot(ToothGrowth, aes(x = supp, y = len)) + geom_boxplot() # the same

# 使用 interaction() 創造更多分組
qplot(interaction(ToothGrowth$supp, ToothGrowth$dose), ToothGrowth$len, geom = "boxplot")
qplot(interaction(supp, dose), len, data = ToothGrowth, geom = "boxplot") # another way
ggplot(ToothGrowth, aes(x = interaction(supp, dose), y = len)) + geom_boxplot() # the same


## 2.6 函數圖形 =====

curve(x ^ 3 - 5 * x, from = -4, to = 4)

myfun <- function(xvar) {
        1 / (1 + exp(-xvar + 10))
}
curve(myfun(x), from = 0, to = 20)
curve(1 - myfun(x), add = TRUE, col = "red")

# qplot(c(0, 20), fun = myfun, stat = "function", geom = "line") # Error
ggplot(data.frame(x = c(0, 20)), aes(x = x)) + stat_function(fun = myfun, geom = "line")


# 3. 長條圖 =====

## 3.1 簡單長條圖 =====

library(gcookbook)
ggplot(pg_mean, aes(x = group, y = weight)) + geom_bar(stat = "identity")

ggplot(BOD, aes(x = Time, y = demand)) + geom_bar(stat = "identity")
ggplot(BOD, aes(x = factor(Time), y = demand)) + geom_bar(stat = "identity")

ggplot(pg_mean, aes(x = group, y = weight)) + geom_bar(stat = "identity", fill = "lightblue", colour = "black")


## 3.2 簇狀條形圖 =====

cabbage_exp
class(cabbage_exp$Cultivar); class(cabbage_exp$Date)
ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) +  # fill要用factor變數
        geom_bar(position = "dodge",  # 水平排列
                 stat = "identity")

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(position = "dodge", stat = "identity", colour = "black") +
        scale_fill_brewer(palette = "Pastel1") # scale_fill_brewer(palette = "")


ce <- cabbage_exp[1:5, ]
ce

ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(position = "dodge", stat = "identity", colour = "black") +
        scale_fill_brewer(palette = "Pastel1")


## 3.3 頻度長條圖 =====

ggplot(diamonds, aes(x = cut)) + geom_bar() # 離散型會得到bar chart

ggplot(diamonds, aes(x = carat)) + geom_bar() # 連續型會得到histogram


## 3.4 長條圖著色 =====

upc <- subset(uspopchange, rank(Change) > 40)
upc

ggplot(upc, aes(x = Abb, y = Change, fill = Region)) + geom_bar(stat = "identity")

ggplot(upc, aes(x = reorder(Abb, Change), y = Change, fill = Region)) + 
        geom_bar(stat = "identity", colour = "black") +
        scale_fill_manual(values = c("#669933", "#FFCC66")) +  # scale_fill_manual()
        xlab("State")


## 3.5 正負長條圖 =====

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


## 3.6 調整長條寬度和間距 =====

ggplot(pg_mean, aes(x = group, y = weight)) + geom_bar(stat = "identity")


ggplot(pg_mean, aes(x = group, y = weight)) + geom_bar(stat = "identity", width = .5)
ggplot(pg_mean, aes(x = group, y = weight)) + geom_bar(stat = "identity", width = 1)

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", width = .5, position = "dodge")

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", width = .5, position = position_dodge(.7))
# position = position_dodge()


## 3.7 堆積長條圖 =====

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity")
cabbage_exp

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity") + 
        guides(fill = guide_legend(reverse = TRUE))
# guides() 調整圖例，調整圖例的填充色順序用 fill = guide_legend(reverse = TRUE)

library(plyr) # 為了使用 order = desc()
ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar, order = desc(Cultivar))) + 
        geom_bar(stat = "identity")
# 調整堆疊順序用 order = desc()

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", colour = "black") + 
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_fill_brewer(palette = "Pastel1")


## 3.8 百分比堆積長條圖 =====

library(plyr)
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


## 3.9 加上資料標籤 =====

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
        geom_text(aes(y = Weight + .1, label = Weight)) # 設定標籤的y軸位置

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label = Weight), vjust = 1.5, colour = "white", 
                  position = position_dodge(.9), size = 3) # 分類間距的default為0.9


# 直向堆積長條圖要加入標籤之前，要先對各組資料求出累積和！
library(plyr)
ce <- arrange(cabbage_exp, Date, Cultivar)
ce <- ddply(ce, "Date", transform, label_y = cumsum(Weight))
ce

ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity") +
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


## 3.10 Cleveland點圖 =====

# Check the name!!!
library(gcookbook)

tophit <- tophitters2001[1:25, ]
ggplot(tophit, aes(x = avg, y = name)) + geom_point()


# 4. 折線圖 =====

## 4.1 簡單折線圖 =====

ggplot(BOD, aes(x = Time, y = demand)) + geom_line()
BOD

BOD1 <- BOD
BOD1$Time <- factor(BOD1$Time) # 連續型變量轉factor
ggplot(BOD1, aes(x = Time, y = demand, group = 1)) + geom_line() # 對於factor變量，要加入 group = 1 確認為同一組資料


# 拓寬y軸範圍
ggplot(BOD, aes(x = Time, y = demand)) + geom_line() + ylim(0, max(BOD$demand))
ggplot(BOD, aes(x = Time, y = demand)) + geom_line() + expand_limits(y = 0)


# 4.2 加上資料標記 =====

ggplot(BOD, aes(x = Time, y = demand)) + geom_line() + geom_point()


library(gcookbook)
ggplot(worldpop, aes(x = Year, y = Population)) + geom_line() + geom_point()
ggplot(worldpop, aes(x = Year, y = Population)) + geom_line() + geom_point() + scale_y_log10()


# 4.3 多重折線圖 =====

library(plyr)
tg <- ddply(ToothGrowth, c("supp", "dose"), summarise, length = mean(len))
ggplot(tg, aes(x = dose, y = length, colour = supp)) + geom_line()
ggplot(tg, aes(x = dose, y = length, linetype = supp)) + geom_line() # 改變線的形狀


# 一定一定要加 group = supp 讓電腦知道資料是一組的
ggplot(tg, aes(x = factor(dose), y = length, colour = supp, group = supp)) + geom_line()

# 錯誤示範
ggplot(tg, aes(x = dose, y = length)) + geom_line() # 沒有正確分組造成一個x對應不只一個點


ggplot(tg, aes(x = dose, y = length, shape = supp)) + geom_line() + 
        geom_point(size = 4)

ggplot(tg, aes(x = dose, y = length, fill = supp)) + geom_line() + 
        geom_point(size = 4, shape = 21)

# 標記可能互相重疊，可以適當地左右移動 position = position_dodge()
ggplot(tg, aes(x = dose, y = length, shape = supp)) + 
        geom_line(position = position_dodge(.2)) + 
        geom_point(position = position_dodge(.2), size = 4)


# 4.4 修改線條樣式 =====

ggplot(BOD, aes(x = Time, y = demand)) + 
        geom_line(linetype = "dashed", size = 1, colour = "blue")

library(plyr)
tg <- ddply(ToothGrowth, c("supp", "dose"), summarise, length = mean(len))
ggplot(tg, aes(x = dose, y = length, colou = supp)) + 
        geom_line() + scale_colour_brewer(palette = "Set1")

ggplot(tg, aes(x = dose, y = length, group = supp)) + 
        geom_line(linetype = "dashed") + 
        geom_point(shape = 22, size = 3, fill = "white")


# 4.5 修改資料標記 =====

ggplot(BOD, aes(x = Time, y = demand)) + geom_line() + geom_point(size = 4, shape = 22, colour = "darkred", fill = "pink")
ggplot(BOD, aes(x = Time, y = demand)) + geom_line() + geom_point(size = 4, shape = 21, fill = "white")

library(plyr)
tg <- ddply(ToothGrowth, c("supp", "dose"), summarise, length = mean(len))
pd <- position_dodge(.2)
ggplot(tg, aes(x = dose, y = length, fill = supp)) + 
        geom_line(position = pd) + 
        geom_point(shape = 21, size = 3, position = pd) +
        scale_fill_manual(values = c("black", "white"))


# 4.6 面積圖 =====

sunspotyear <- data.frame(
        Year = as.numeric(time(sunspot.year)),
        Sunspots = as.numeric(sunspot.year))
ggplot(sunspotyear, aes(x = Year, y = Sunspots)) + geom_area()

ggplot(sunspotyear, aes(x = Year, y = Sunspots)) + 
        geom_area(colour = "black", fill = "blue", alpha = .2) # 調整填充色

ggplot(sunspotyear, aes(x = Year, y = Sunspots)) + 
        geom_area(fill = "blue", alpha = .2) +  # 不設定colour就沒有底部框線
        geom_line()


# 4.7 堆積面積圖 =====

ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup)) + geom_area()

# default 的堆積順序與圖標有時是相反的，可以用 breaks 參數調整
ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup)) + geom_area(colour = "black", size = .2, alpha = .4) +
        scale_fill_brewer(palette = "Blues", breaks = rev(levels(uspopage$AgeGroup)))

library(plyr)
# 可以用 order = desc() 反轉堆積順序
ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup, order = desc(AgeGroup))) +
        geom_area(colour = "black", size = .2, alpha = .4) + 
        scale_fill_brewer(palette = "Blues")

ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup, order = desc(AgeGroup))) +
        geom_area(colour = NA, alpha = .4) +
        scale_fill_brewer(palette = "Blues") +
        geom_line(position = "stack", size = .2) # 堆積面積圖加線


# 4.8 百分比堆積面積圖 =====

library(gcookbook)
library(plyr)
uspopage_prop <- ddply(uspopage, "Year", transform, Percent = Thousands / sum(Thousands) * 100)
# ddply 將 uspopage 資料集，按照 Year 拆成多個獨立的 data frames

ggplot(uspopage_prop, aes(x = Year, y = Percent, fill = AgeGroup)) + 
        geom_area(colour = "black", size = .2, alpha = .4) + 
        scale_fill_brewer(palette = "Blues", breaks = rev(levels(uspopage$AgeGroup)))


# 4.9 增加信賴區間 =====

clim <- subset(climate, Source == "Berkeley", select = c("Year", "Anomaly10y", "Unc10y"))
clim
# geom_ribbon() 要先畫，才不會讓 geom_line() 糊掉！
ggplot(clim, aes(x = Year, y = Anomaly10y)) + 
        geom_ribbon(aes(ymin = Anomaly10y - Unc10y, ymax = Anomaly10y + Unc10y), alpha = .2) +
        geom_line()

ggplot(clim, aes(x = Year, y = Anomaly10y)) + 
        geom_line(aes(y = Anomaly10y - Unc10y), colour = "grey50", linetype = "dotted") +
        geom_line(aes(y = Anomaly10y + Unc10y), colour = "grey50", linetype = "dotted") +
        geom_line()


# 5. 散佈圖 =====


# 5.1 散佈圖 =====

library(gcookbook)
heightweight[, c("ageYear", "heightIn")]
ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point()
ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point(shape = 21)
ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point(size = 1.5)


# 5.2 修改點的樣式 =====

heightweight[, c("sex", "ageYear", "heightIn")]

ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + geom_point()
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex)) + geom_point()

ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, colour = sex)) + geom_point()

ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, colour = sex)) + geom_point() + 
        scale_shape_manual(values = c(1, 2)) +  # 人工選擇點的形狀
        scale_colour_brewer(palette = "Set1") # 人工選擇點的顏色


# 5.3 使用非內建的點形 =====

ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point(shape = 3)
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex)) + geom_point(size = 3) + scale_shape_manual(values = c(1, 4))


hw <- heightweight
hw$weightGroup <- cut(hw$weightLb, breaks = c(-Inf, 100, Inf), labels = c("< 100", ">= 100"))
ggplot(hw, aes(x = ageYear, y = heightIn, shape = sex, fill = weightGroup)) + 
        geom_point(size = 2.5) +
        scale_shape_manual(values = c(21, 24)) +
        scale_fill_manual(values = c(NA, "black"), 
                          guide = guide_legend(override.aes = list(shape = 21)))


# 5.4 將連續型變數映射到點的顏色或大小 =====

heightweight[, c("sex", "ageYear", "heightIn", "weightLb")]

ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = weightLb)) + geom_point()
ggplot(heightweight, aes(x = ageYear, y = heightIn, size = weightLb)) + geom_point()

# 人天生對顏色和大小的變化不太敏銳，因此當一個變數不需要太精密的解釋時才適合用

ggplot(heightweight, aes(x = ageYear, y = heightIn, fill = weightLb)) + 
        geom_point(shape = 21, size = 2.5) + 
        scale_fill_gradient(low = "black", high = "white") # 黑白漸層

ggplot(heightweight, aes(x = ageYear, y = heightIn, fill = weightLb)) + 
        geom_point(shape = 21, size = 2.5) +
        scale_fill_gradient(low = "black", high = "white", breaks = seq(70, 170, by = 20), guide = guide_legend())
# 黑白漸層搭配離散圖標

ggplot(heightweight, aes(x = ageYear, y = heightIn, size = weightLb, colour = sex)) + 
        geom_point(alpha = .5) +
        scale_size_area() +  # 使面積與資料值成正比
        scale_colour_brewer(palette = "Set1")

# 不同形狀的點很難比較面積大小，因此不要同時操作


# 5.5 處理資料點重疊問題 =====

# 圖形重疊（overplotting）的解決方案：
## 半透明的點
## 矩形資料分箱
## 六邊形資料分箱
## 箱型圖

sp <- ggplot(diamonds, aes(x = carat, y = price))
sp + geom_point()

sp + geom_point(alpha = .1) # 90%透明度
sp + geom_point(alpha = .01) # 99%透明度


sp + stat_bin2d() # 分別在x軸y軸分割30組，共900個箱子
sp + stat_bin2d(bins = 50) + scale_fill_gradient(low = "lightblue", high = "red", limits = c(0, 6000)) # 2500箱


library(hexbin) # 使用六邊形箱子
sp + stat_binhex() + scale_fill_gradient(low = "lightblue", high = "red", limits = c(0, 8000))
sp + stat_binhex() + scale_fill_gradient(low = "lightblue", high = "red", breaks = c(0, 250, 500, 1000, 2000, 4000, 6000), limits = c(0, 6000))
# 範圍外會變成灰色箱子


# 當其中一軸或兩軸為離散型變數時，也會出現overplotting，可用 position_jitter() 增加隨機擾動

sp1 <- ggplot(ChickWeight, aes(x = Time, y = weight))
sp1 + geom_point()
sp1 + geom_point(position = "jitter") # 跟 position_jitter() 意思一樣

sp1 + geom_point(position = position_jitter(width = .5, height = 0)) # width 和 height 調整擾動值的精度

sp1 + geom_boxplot(aes(group = Time))


# 5.6 模型擬合線 =====

library(gcookbook)
sp <- ggplot(heightweight, aes(x = ageYear, y = heightIn))
sp + geom_point() + stat_smooth(method = lm) # default 為95%信賴區間

sp + geom_point() + stat_smooth(method = lm, level = .99) # 調整 level 可調整信賴區間

sp + geom_point(colour = "grey60") + stat_smooth(method = lm, se = FALSE, colour = "black")

# stat_smooth 的 default 是 loess曲線
sp + geom_point(colour = "grey60") + stat_smooth()
sp + geom_point(colour = "grey60") + stat_smooth(method = loess) # the same


library(MASS)
b <- biopsy
b$classn[b$class == "benign"] <- 0
b$classn[b$class == "malignant"] <- 1
b

ggplot(b, aes(x = V1, y = classn)) + 
        geom_point(position = position_jitter(width = .3, height = .06), alpha = .4, shape = 21, size = 1.5) + 
        stat_smooth(method = glm, method.args = list(family = "binomial"))


sps <- ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + 
        geom_point() + 
        scale_colour_brewer(palette = "Set1")

sps + stat_smooth() # stat_smooth() 的範圍限定在預測資料對應的範圍內

sps + stat_smooth(method = lm, se = FALSE, fullrange = TRUE) # 可外推的模型要加入參數 fullrange = TRUE


# 5.7 既有模型散佈圖加入擬合線 =====

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
        if (is.null(xrange)) {
                if (any(class(model) %in% c("lm", "glm")))
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


# 5.8 多模型擬合線 =====

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


# 5.9 散佈圖加入模型係數 =====

library(gcookbook)
model <- lm(heightIn ~ ageYear, heightweight)
summary(model)

pred <- predictvals(model, "ageYear", "heightIn")
sp <- ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point() + geom_line(data = pred)
sp + annotate("text", label = "r^2 = 0.42", x = 16.5, y = 52) # 加入係數標籤
sp + annotate("text", label = "r^2 == 0.42", parse = TRUE, x = 16.5, y = 52) # 用R的方式表示數學符號


eqn <- as.character(as.expression(
        substitute(italic(y) == a + b * italic(x) * "," ~~ italic(r)^2 ~ "=" ~ r2, 
                   list(a = format(coef(model)[1], digits = 3),
                        b = format(coef(model)[2], digits = 3),
                        r2 = format(summary(model)$r.squared, digits = 2)
                        ))))
eqn

parse(text = eqn)


sp + annotate("text", label = eqn, parse = TRUE, x = Inf, y = -Inf, hjust = 1.1, vjust = -.5)
# x = Inf, y = -Inf 使公式置於右下角


# 5.10 散佈圖加入邊際地毯 =====

ggplot(faithful, aes(x = eruptions, y = waiting)) + geom_point() + geom_rug()
# marginal rugs 本質上是一個一維的散佈圖
ggplot(faithful, aes(x = eruptions, y = waiting)) + geom_point() + 
        geom_rug(position = "jitter", size = .2)
# position = "jitter", size = .2 調整線寬及減輕重疊程度


# 5.11 散佈圖加標籤 =====

library(gcookbook)
subset(countries, Year == 2009 & healthexp > 2000)
sp <- ggplot(subset(countries, Year == 2009 & healthexp > 2000), 
             aes(x = healthexp, y = infmortality)) + geom_point()
sp + annotate("text", x = 4350, y = 5.4, label = "Canada") +
        annotate("text", x = 7400, y = 6.8, label = "USA")


sp + geom_text(aes(label = Name), size = 4) # geom_text() 可用factor或char類型的向量製作標籤


sp + geom_text(aes(label = Name), size = 4, vjust = 0)
sp + geom_text(aes(y = infmortality + .1, label = Name), size = 4, vjust = 0)

# 左對齊 hjust = 0; 右對齊 hjust = 1
# 用這種方法，較長的標籤會有較大的移動，此時最好用 x 增減一個值來調整

sp + geom_text(aes(label = Name), size = 4, hjust = 0)
sp + geom_text(aes(x = healthexp + 100, label = Name), size = 4, hjust = 0)


# 如果不想要全部加上標籤，可以複製一個新的標籤
cdat <- subset(countries, Year == 2009 & healthexp > 2000)
cdat$Name1 <- cdat$Name

idx <- cdat$Name1 %in% c("Canada", "Ireland", "United Kingdom", "United States", "New Zealand", "Iceland",
                         "Japan", "Luxembourg", "Netherlands", "Switzerland")
idx
cdat$Name1[!idx] <- NA
cdat

ggplot(cdat, aes(x = healthexp, y = infmortality)) + geom_point() + 
        geom_text(aes(x = healthexp + 100, label = Name1), size = 4, hjust = 0) + xlim(2000, 10000)


# 5.12 氣泡圖 =====

library(gcookbook)
cdat <- subset(countries, Year == 2009 & Name %in% c("Canada", "Ireland", "United Kingdom", 
                                                     "United States", "New Zealand", "Iceland","Japan", 
                                                     "Luxembourg", "Netherlands", "Switzerland"))
cdat
p <- ggplot(cdat, aes(x = healthexp, y = infmortality, size = GDP)) + geom_point(shape = 21, colour = "black", fill = "cornsilk")
p
p + scale_size_area(max_size = 15) # 以GDP決定面積


# 當x軸y軸都是類別變數，氣泡圖可以用來表示變量值
hec <- HairEyeColor[,, "Male"] + HairEyeColor[,, "Female"] # HairEyeColor是個list
library(reshape2)
hec <- melt(hec, value.name = "count")
ggplot(hec, aes(x = Eye, y = Hair)) + geom_point(aes(size = count), shape = 21, colour = "black", fill = "cornsilk") +
        scale_size_area(max_size = 20, guide = FALSE) + 
        geom_text(aes(y = as.numeric(Hair) - sqrt(count) / 22, label = count), vjust = 1, colour = "grey60", size = 4)
# 此處的y座標是計算得出


# 5.13 散佈圖矩陣 =====

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


# 6. 敘述統計 =====

# 6.1 簡單直方圖 =====

ggplot(faithful, aes(x = waiting)) + geom_histogram()

w <- faithful$waiting
ggplot(NULL, aes(x = w)) + geom_histogram()

ggplot(faithful, aes(x = waiting)) + 
        geom_histogram(binwidth = 5, 
                       fill = "white",  # 填滿
                       colour = "black") # 邊框顏色

binsize <- diff(range(faithful$waiting)) / 15 # 計算binwidth
ggplot(faithful, aes(x = waiting)) + geom_histogram(binwidth = binsize, fill = "white", colour = "black")


h <- ggplot(faithful, aes(x = waiting)) # 儲存成變量以便重複利用

h + geom_histogram(binwidth = 8, fill = "white", colour = "black", origin = 31) # 設定組邊界
h + geom_histogram(binwidth = 8, fill = "white", colour = "black", origin = 35)

# geom_bar(stat = "bin") 可得相同結果


## 6.2 多組資料直方圖 =====

library(MASS) # 取用數據
ggplot(birthwt, aes(x = bwt)) + geom_histogram(fill = "white", colour = "black") + 
        facet_grid(smoke ~ .)

birthwt1 <- birthwt
birthwt1$smoke <- factor(birthwt1$smoke)
levels(birthwt1$smoke) # 標籤是0, 1不方便

library(plyr) # 為了使用 revalue()
birthwt1$smoke <- revalue(birthwt1$smoke, c("0" = "No Smoke", "1" = "Smoke"))

ggplot(birthwt1, aes(x = bwt)) + geom_histogram(fill = "white", colour = "black") +
        facet_grid(smoke ~ .)


ggplot(birthwt, aes(x = bwt)) + geom_histogram(fill = "white", colour = "black") + facet_grid(race ~ .)
ggplot(birthwt, aes(x = bwt)) + geom_histogram(fill = "white", colour = "black") + 
        facet_grid(race ~ ., scales = "free") # scales = "free" 可以單獨設定各y軸尺度

birthwt1$smoke <- factor(birthwt1$smoke)
ggplot(birthwt1, aes(x = bwt, fill = smoke)) + 
        geom_histogram(position = "identity",  # 使直方圖重疊
                       alpha = .4) # 設定透明度


## 6.3 密度曲線 =====

ggplot(faithful, aes(x = waiting)) + geom_density()

ggplot(faithful, aes(x = waiting)) + 
        geom_line(stat = "density") +  # 不想要邊線和底線可以用 geom_line()
        expand_limits(y = 0) # 擴大y軸範圍，包含0


w <- faithful$waiting
ggplot(NULL, aes(x = w)) + geom_density() # 可以傳遞向量作為參數


ggplot(faithful, aes(x = waiting)) +  # 帶寬不同，曲線的光滑程度不同
        geom_line(stat = "density", adjust = .25, colour = "red") +
        geom_line(stat = "density") + 
        geom_line(stat = "density", adjust = 2, colour = "blue")


ggplot(faithful, aes(x = waiting)) + geom_density(fill = "blue", alpha = .2) + 
        xlim(35, 105) # 限制x軸範圍

ggplot(faithful, aes(x = waiting)) + geom_density(fill = "blue", colour = NA, alpha = .2) +
        geom_line(stat = "density") + xlim(35, 105)


ggplot(faithful, aes(x = waiting, y = ..density..)) +  # y = ..density.. 可以縮小直方圖的尺度以配合密度曲線
        geom_histogram(fill = "cornsilk", colour = "grey60", size = .2) +
        geom_density() + 
        xlim(35, 105)


## 6.4 多組資料密度曲線 =====

library(MASS)
birthwt1 <- birthwt
birthwt1$smoke <- factor(birthwt1$smoke) # 必須先轉化為factor
ggplot(birthwt1, aes(x = bwt, colour = smoke)) + geom_density()

ggplot(birthwt1, aes(x = bwt, fill = smoke)) + geom_density(alpha = .3)


ggplot(birthwt1, aes(x = bwt)) + geom_density() + facet_grid(smoke ~ .)
levels(birthwt1$smoke)
library(plyr)
birthwt1$smoke <- revalue(birthwt1$smoke, c("0" = "No Smoke", "1" = "Smoke"))
ggplot(birthwt1, aes(x = bwt)) + geom_density() + facet_grid(smoke ~ .)

ggplot(birthwt1, aes(x = bwt, y = ..density..)) + 
        geom_histogram(binwidth = 200, fill = "cornsilk", colour = "grey60", size = .2) +
        geom_density() + 
        facet_grid(smoke ~ .)


## 6.5 頻次多邊形 =====

# freqpoly 跟 histogram 非常類似，而核密度曲線只是一個估計

ggplot(faithful, aes(x = waiting)) + geom_freqpoly()

ggplot(faithful, aes(x = waiting)) + geom_freqpoly(binwidth = 4) # 控制組距

binsize <- diff(range(faithful$waiting)) / 15
ggplot(faithful, aes(x = waiting)) + geom_freqpoly(binwidth = binsize)


## 6.6 基本箱型圖（盒鬚圖） =====

library(MASS)
ggplot(birthwt, aes(x = factor(race), y = bwt)) + geom_boxplot()

ggplot(birthwt, aes(x = factor(race), y = bwt)) + geom_boxplot(width = .5) # 調整箱型圖寬度

ggplot(birthwt, aes(x = factor(race), y = bwt)) + 
        geom_boxplot(outlier.size = 1.5, outlier.shape = 21) # 修改outlier的外觀


ggplot(birthwt, aes(x = 1, y = bwt)) + geom_boxplot() + 
        scale_x_continuous(breaks = NULL) +  #只有1組，x軸不標記
        theme(axis.title.x = element_blank())


## 6.7 缺口箱型圖 =====

ggplot(birthwt, aes(x = factor(race), y = bwt)) + geom_boxplot(notch = TRUE)

# 若各組缺口不重疊，表示各組中位數有差異


## 6.8 箱型圖加平均值 =====

ggplot(birthwt, aes(x = factor(race), y = bwt)) + geom_boxplot() + 
        stat_summary(fun.y = "mean", geom = "point", shape = 23, size = 3, fill = "white")


## 6.9 小提琴圖 =====

# install.packages("gcookbook")
library(gcookbook)
p <- ggplot(heightweight, aes(x = sex, y = heightIn))
p + geom_violin()

p + geom_violin() +  # 同時畫小提琴圖跟箱型圖
        geom_boxplot(width = .1, fill = "black", outlier.colour = NA) +
        stat_summary(fun.y = median, geom = "point", fill = "white", shape = 21, size = 2.5)


p + geom_violin(trim = FALSE)

# default情況是不同組的小提琴圖面積會相同，可以用scale = "count"調整

p + geom_violin(scale = "count")

p + geom_violin(adjust = 2) # adjust修改平滑程度
p + geom_violin(adjust = .5)


## 6.10 Wilkinson點圖 =====

library(gcookbook)
countries2009 <- subset(countries, Year == 2009 & healthexp > 2000)

p <- ggplot(countries2009, aes(x = infmortality))
p + geom_dotplot()

p + geom_dotplot(binwidth = .25) + geom_rug() +  # geom_rug() 標示資料點的具體位置
        scale_y_continuous(breaks = NULL) +  # 移除y軸尺度
        theme(axis.title.y = element_blank()) # 移除座標軸標籤
# 數據在水平方向上並非均勻分佈，而是堆在他表示的資料點的中心位置

p + geom_dotplot(method = "histodot", binwidth = .25) + geom_rug() +
        scale_y_continuous(breaks = NULL) + theme(axis.title.y = element_blank())
# 用像直方圖的固定間距分組算法


# 奇數數量與偶數數量保持一致的中心堆疊方式，設定stackdir = "center"或stackdir = "centerwhole"
p + geom_dotplot(binwidth = .25, stackdir = "center") +
        scale_y_continuous(breaks = NULL) + theme(axis.title.y = element_blank())
p + geom_dotplot(binwidth = .25, stackdir = "centerwhole") +
        scale_y_continuous(breaks = NULL) + theme(axis.title.y = element_blank())


## 6.11 多組資料點圖 =====

library(gcookbook)
ggplot(heightweight, aes(x = sex, y = heightIn)) + geom_dotplot(binaxis = "y", binwidth = .5, stackdir = "center")

# 同時畫箱型圖跟點圖
ggplot(heightweight, aes(x = sex, y = heightIn)) + geom_boxplot(outlier.colour = NA, width = .4) + 
        geom_dotplot(binaxis = "y", binwidth = .5, stackdir = "center", fill = NA)


# 把箱型圖畫在點圖旁邊
ggplot(heightweight, aes(x = sex, y = heightIn)) + 
        geom_boxplot(aes(x = as.numeric(sex) + .2, group = sex), width = .25) +
        geom_dotplot(aes(x = as.numeric(sex) - .2, group = sex), binaxis = "y", binwidth = .5, stackdir = "center") +
        scale_x_continuous(breaks = 1:nlevels(heightweight$sex), labels = levels(heightweight$sex))


## 6.12 二維資料密度圖 =====

p <- ggplot(faithful, aes(x = eruptions, y = waiting))
p + geom_point() + stat_density2d()


p + stat_density2d(aes(colour = ..level..)) # 把密度曲面的高度映射給等高線的顏色


# 將密度估計映射給填充色或瓦片圖的透明度
# 柵格: raster
# 瓦片: tile
p + stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE)
p + geom_point() + stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE)

p + stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE, h = c(.5, 5))


# 7. 註解 =====

## 7.1 添加文本註解 =====

p <- ggplot(faithful, aes(x= eruptions, y = waiting)) + geom_point()
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
# Inf, -Inf 在繪圖區邊緣放置文本，hjust, vjust 調整文本相對於邊緣的位置


## 7.2 在註解中用數學方程式 =====

p <- ggplot(data.frame(x = c(-3, 3)), aes(x = x)) + stat_function(fun = dnorm)
p + annotate("text", x = 2, y = .3, parse = TRUE, label = "frac(1, sqrt(2 * pi)) * e ^ {-x ^ 2 / 2}")
# 加 parse = TRUE
# 內部引號閉合的每一部分都會被當作數學方程式的一個變量，要讓乘號 * 變可見要用 %*%

p + annotate("text", x = 0, y = .05, parse = TRUE, size = 4, label = "'Function: ' * y == frac(1, sqrt(2 * pi)) * e ^ {-x ^ 2 / 2}")
# ?plotmath 看更多數學方程
# ?demo(plotmath) 看數學方程的圖示


## 7.3 添加直線 =====

library(gcookbook)
p <- ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + geom_point()
p + geom_hline(yintercept = 60) + geom_vline(xintercept = 14)
p + geom_abline(intercept = 37.4, slope = 1.75)

library(plyr)
hw_means <- ddply(heightweight, "sex", summarise, heightIn = mean(heightIn))
hw_means

p + geom_hline(aes(yintercept = heightIn, colour = sex), data = hw_means, linetype = "dashed", size = 1)


# 若 x 軸為離散型變數
pg <- ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_point()
pg + geom_vline(xintercept = 1)
pg + geom_vline(xintercept = which(levels(PlantGrowth$group) == "ctrl"))


## 7.4 添加線段和箭頭 =====

library(gcookbook)
p <- ggplot(subset(climate, Source == "Berkeley"), aes(x = Year, y = Anomaly10y)) + geom_line()
p + annotate("segment", x = 1950, xend = 1980, y = -.25, yend = -.25)

library(grid)
p + annotate("segment", x = 1850, xend = 1820, y = -.8, yend = -.95, colour = "blue", size = 2, arrow = arrow()) + 
        annotate("segment", x = 1950, xend = 1980, y = -.25, yend = -.25, arrow = arrow(ends = "both", angle = 90, length = unit(.2, "cm")))
# 前者為30度箭頭，後者為90度雙箭頭


## 7.5 添加矩形陰影 =====

library(gcookbook)
p <- ggplot(subset(climate, Source == "Berkeley"), aes(x = Year, y = Anomaly10y)) + geom_line()
p + annotate("rect", xmin = 1950, xmax = 1980, ymin = -1, ymax = 1, alpha = .1, fill = "blue")


## 7.6 高亮某一元素 =====

pg <- PlantGrowth
pg$h1 <- "no"
pg$h1[pg$group == "trt2"] <- "yes"
ggplot(pg, aes(x = group, y = weight, fill = h1)) + geom_boxplot() + 
        scale_fill_manual(values = c("grey85", "#FFDDCC"), guide = FALSE)

ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) + geom_boxplot() + 
        scale_fill_manual(values = c("grey85", "grey85", "#FFDDCC"), guide = FALSE)


## 7.7 添加誤差線 =====

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
        geom_errorbar(aes(ymin = Weight - se, ymax = Weight + se), position = "dodge", width = .2)

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(position = "dodge", stat = "identity") + 
        geom_errorbar(aes(ymin = Weight - se, ymax = Weight + se), position = position_dodge(.9), width = .2)

# position = "dodge"即為position = position_dodge()的縮寫版本

pd <- position_dodge(.3)
ggplot(cabbage_exp, aes(x = Date, y = Weight, colour = Cultivar, group = Cultivar)) + 
        geom_errorbar(aes(ymin = Weight - se, ymax = Weight + se), width = .2, size = .25, colour = "black", position = pd) + 
        geom_line(position = pd) + 
        geom_point(position = pd, size = 2.5)


## 7.8 向獨立分面添加註解 =====

p <- ggplot(mpg, aes(x = displ, y = hwy)) + geom_point() + facet_grid(. ~ drv)
f_labels <- data.frame(drv = c("4", "f", "r"), label = c("4wd", "Front", "Rear"))
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
        geom_text(x = 3, y = 34, aes(label = r2), data = labels, parse = TRUE, hjust = 0)

labels <- ddply(mpg, "drv", summarise, r2 = cor(displ, hwy) ^ 2)
labels$r2 <- sprintf("italic(R ^ 2) == %.2f", labels$r2)


# 8. 座標軸 =====

## 8.1 交換x軸和y軸 =====

ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot()
ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot() + coord_flip()

ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot() + coord_flip() + 
        scale_x_discrete(limits = rev(levels(PlantGrowth$group))) # 反轉因子的排列順序


## 8.2 設置連續型座標軸的值域 =====

p <- ggplot(PlantGrowth, aes(x = group, y= weight)) + geom_boxplot()
p

p + ylim(0, max(PlantGrowth$weight))


ylim(0, 10)
scale_y_continuous(limits = c(0, 10)) # 兩式等價


p + ylim(0, 10) + scale_y_continuous(breaks = c(0, 5, 10))
p + scale_y_continuous(breaks = c(0, 5, 10)) + ylim(0, 10) # 同時使用會產生錯誤
p + scale_y_continuous(limits = c(0, 10), breaks = c(0, 5, 10)) # 想同時生效可加在參數

p + scale_y_continuous(limits = c(5, 6.5)) # 限制y的值域會使得部分資料被剪除
p + coord_cartesian(ylim = c(5, 6.5)) # 此法純粹縮放


p + expand_limits(y = 0) # 擴展值域


# 8.3 反轉一條連續型座標軸 =====

ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot() + 
        scale_y_reverse() # y軸刻度反轉

ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot() + ylim(6.5, 3.5)


ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot() + scale_y_reverse(limits = c(8, 0))
# ylim()與scale_y_reverse()同樣不能配合，要設定在參數


## 8,4 修改類別型座標軸上項目的順序 =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot()
p + scale_x_discrete(limits = c("trt1", "ctrl", "trt2"))


p + scale_x_discrete(limits = c("ctrl", "trt1")) # 只展示部分類別
p + scale_x_discrete(limits = rev(levels(PlantGrowth$group))) # 反轉類別順序


## 8.5 設置x軸和y軸的縮放比例 =====

library(gcookbook)
sp <- ggplot(marathon, aes(x = Half, y = Full)) + geom_point()
sp + coord_fixed() # 等長縮放座標軸

sp + coord_fixed() + 
        scale_y_continuous(breaks = seq(0, 420, 30)) + 
        scale_x_continuous(breaks = seq(0, 420, 30)) # 指定位置放刻度線

sp + coord_fixed(ratio = 1/2) +  # 讓座標軸延伸
        scale_y_continuous(breaks = seq(0, 420, 30)) + 
        scale_x_continuous(breaks = seq(0, 420, 15))

## 8.6 設置刻度線的位置 =====

ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot()
ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot() + 
        scale_y_continuous(breaks = c(4, 4.25, 4.5, 5, 6, 8))

seq(4, 7, by = .5)

ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot() + 
        scale_x_discrete(limits = c("trt2", "ctrl"), breaks = "ctrl") # 用limits篩選組別


## 8.7 移除刻度線和標籤 =====

p <- ggplot(PlantGrowth, aes(x = group, y = weight)) + geom_boxplot()
p + theme(axis.text.y = element_blank()) # 移除刻度標籤

p + theme(axis.ticks = element_blank(), axis.text.y = element_blank()) # 移除兩軸刻度線

p + scale_y_continuous(breaks = NULL) # 移除刻度線、刻度標籤和網格線


## 8.8 修改刻度標籤的文本 =====

library(gcookbook)
hwp <- ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point()
hwp + scale_y_continuous(breaks = c(50, 56, 60, 66, 72), labels = c("Tiny", "Really\nshort", "Short", "Medium", "Tallish"))
# \n 指的是換行！


footinch_formatter <- function(x) {
        foot <- floor(x / 12)
        inch <- x %% 12
        return(paste(foot, "'", inch, "\"", sep = "")) # \"加反斜線是為了區分引號字符與參數引號
}


# 9. 整體外觀 =====

# 10. 圖例 =====

# 11. 分面 =====

# 12. 配色 =====

# 13. 其他圖形 =====

## 13.1 相關係數矩陣圖 =====

mcor <- cor(mtcars)
round(mcor, digits = 2)
library(corrplot)
corrplot(mcor)

corrplot(mcor, method = "shade", shade.col = NA, tl.col = "black", tl.srt = 45)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# corrplot(mcor, method = "shade", shade.col = NA, tl.col = "black", 
#         tl.srt = 45, col = col(200), addCoef.col = "black", 
#         cl.pos = "no", order = "AOE") # Error


# 14. 輸出圖形 =====

# 15. reshape =====

library(gcookbook)
heightweight
str(heightweight)


# 15.1 創建 data.frame =====

g <- c("A", "B", "C")
x <- 1:3
dat <- data.frame(g, x)
dat

lst <- list(group = g, value = x)
dat <- as.data.frame(lst)


# 15.2 從數據框中提取訊息 =====


# 15.2 提取資料
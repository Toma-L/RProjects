# FB爬蟲 =====

library(Rfacebook)
#從此取得token: https://developers.facebook.com/tools/explorer 
token <- "EAACEdEose0cBALmokzrbBMQGYzekvU4PZBMgcy9qtTu5ZA9PZAfbXI8tFE1VKIMdOnrhq1vTZAjZB7POFjvZCeeVEX61Dev9gNlZBNljNuoEqMAmGzTxr9ZBvssc74xeaFoVjuzAknTfmLsRcqZB5OMdwsBbvLm3P73iP9VCEIzF5UwZDZD"

fb_oauth <- fbOAuth(app_id = "1758018904467958", app_secret = "6d6bd9221d7562cc5a23e00266fe582b")
save(fb_oauth, file = "fb_oauth")
load("fb_oauth")
me <- getUsers("me", token=fb_oauth)

# getGroup =====

load("fb_oauth")
ids <- searchGroup(name = "字嗨", token = token)
ids[1,] # id = 18533493739
group <- getGroup(group_id = 149874075167476, n = 30000, token = token, since = '2011/01/01', until='2016/07/24')
dim(group)
head(group)
tail(group)

names(group)[8:10] <- c("like", "comment", "share")
names(group)[1:4] <- c("id", "name", "message", "time")

group[group$name == "黃治豪", ]

write.csv(group, "enjoyfonts.csv", row.names = F, fileEncoding = "Big-5")


# getPage =====

page <- getPage(page = "bwnet.fans", token = token, n = 80, since = NULL, until = NULL, feed = FALSE)
names(page)
head(page)
head(page[order(page$shares_count, decreasing = TRUE), c("message", "shares_count")], 10)

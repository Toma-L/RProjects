#Getting and Cleaning Data==================================================

#Reading local flat files

if(!file.exists("data")) {dir.create("data")}
fileUrl <- "http://data.baltimorecity.gov/api/views/dz54-2aru/rows.csv?accessType=DOWNLOAD"
download.file(fileUrl, destfile = "cameras.csv")
dateDownloaded <- date()

cameraData <- read.table("./data/cameras.csv")
head(cameraData)

cameraData <- read.table("./data/cameras.csv", sep = ",", header = TRUE)
head(cameraData)

cameraData <- read.csv("./data/cameras.csv")
head(cameraData)

#Reading Excel files

fileUrl <- "http://data.baltimorecity.gov/api/views/dz54-2aru/rows.xlsx?accessType=DOWNLOAD"
download.file(fileUrl, destfile = "./data/cameras.xlsx")
dateDownloaded <- date()

library(xlsx)
cameraData <- read.xlsx("./data/cameras.xlsx", sheetIndex = 1, header = TRUE)
head(cameraData)

colIndex <- 2:3
rowIndex <- 1:4
cameraDataSubset <- read.xlsx("./data/cameras.xlsx", sheetIndex = 1, colIndex = colIndex, rowIndex = rowIndex)
cameraDataSubset

#Reading XML

library(XML)
fileUrl <- "http://www.w3schools.com/xml/simple.xml"
doc <- xmlTreeParse(fileUrl, useInternal = TRUE)
rootNode <- xmlRoot(doc)
xmlName(rootNode)

rootNode[[1]]
rootNode[[1]][[1]]

xmlSApply(rootNode, xmlValue)
xmlSApply(rootNode, xmlValue)


xpathSApply(rootNode, "//name", xmlValue)

fileUrl <- "http://espn.go.com/nfl/team/_/name/bal/baltimore-ravens"
doc <- htmlTreeParse(fileUrl, useINternal = TURE)
scores <- xpathSApply(doc, "//li[@class='score']", xmlValue)
teams <- xpathSApply(doc, "//li[@class='team'name']", xmlValue)
scores
teams

#Reading JSON

library(jsonlite)
jsonData <- fromJSON("http://api.github.com/users/jtleek/repos")
names(jsonData)
names(jsonData$owner)
jsonData$owner$login

myjson <- toJSON(iris, pretty = TRUE) #把資料轉成JSON，讓某些api用，pretty會編排精美
cat(myjson) #cat就是print()

iris2 <- fromJSON(myjson) #從JSON變data.frame
head(iris2)

#Reading mySQL

ucscDb <-dbConnect(MySQL(), user = "genome", host = "genome-mysql.cse.ucsc.edu")
result <- dbGetQuery(ucscDb, "show databases;")
dbDisconnect(ucscDb) #每次擷取完資料要記得輸入這個！！！
result

hg19 <- dbConnect(MySQL(), user = "genome", db = "hg19", host = "genome-mysql.cse.ucsc.edu")
allTables <- dbListTables(hg19) #dbListTables()
length(allTables)
allTables[1:5]

dbListFields(hg19, "affyU133Plus2") #從hg19資料庫，連上affyU133Plus2表格
dbGetQuery(hg19, "select count(*) from affyU133Plus2") #dbGetQuery()

affyData <- dbReadTable(hg19, "affyU133Plus2") #dbReadTable()
head(affyData)

query <- dbSendQuery(hg19, "select * from affyU133Plus2 where misMatches between 1 and 3") #dbSendQuery
affyMis <- fetch(query) #還沒結束query
quantile(aff)
affyMisSmall <- fetch(query, n = 10)
dbClearResult(query) #停止query
dim(affyMisSmall)

dbDisconnect(hg19)

#Reading HDF5

source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
library(rhdf5)
created <- h5createFile("example.h5")
created

created <- h5createGroup("example.h5", "foo")
created <- h5createGroup("example.h5", "baa")
created <- h5createGroup("example.h5", "foo/foobaa")
h5ls("example.h5")

A <- matrix(1:10, nr = 5, nc = 2)
h5write(A, "example.h5", "foo/A")
B <- array(seq(0.1, 2.0, by = .1), dim = c(5, 2, 2))
attr(B, "scale") <- "liter"
h5write(B, "example.h5", "foo/foobaa/B")
h5ls("example.h5")

df <- data.frame(1L:5L, seq(0, 1, length.out = 5), 
                 c("ab", "cde", "fghi", "a", "s"), stringsAsFactors = FALSE)
h5write(df, "example.h5", "df")
h5ls("example.h5")
readA <- h5read("example.h5", "foo/A")
readB <- h5read("example.h5", "foo/foobaa/B")
readdf <- h5read("example.h5", "df")
readA

h5write(c(12, 13, 14), "example.h5", "foo/A", index = list(1:3, 1))
h5read("example.h5", "foo/A")

#Reading data from the web

con <- url("http://scholar.google.com/citations?user=HI-I6C0AAAAJ&hl=en")
htmlCode <- readLines(con)
close(con)
htmlCode

library(XML)
url <- "http://scholar.google.com/citations?user=HI-I6C0AAAAJ&hl=en"
html <- htmlTreeParse(url, useInternalNodes = TRUE)
xpathSApply(html, "//title", xmlValue)
xpathSApply(html, "//td[@id='col-citedby']", xmlValue)

library(httr)
html2 <- GET(url)
content2 <- content(html2, as = "text")
parsedHtml <- htmlParse(content2, asText = TRUE)
xpathSApply(parsedHtml, "//title", xmlValue)

pg1 <- GET("http://httpbin.org/basic-auth/user/passwd")
pg1

pg2 <- GET("http://httpbin.org/basic-auth/user/passwd", authenticate("user", "passwd"))
pg2

google <- handle("http://google.com")
pg1 <- GET(handle = google, path = "/")
pg2 <- GET(handle = google, path = "search")

#Reading data from APIs

myapp <- oauth_app("twitter", key = "yourConsumerKeyHere", secret = "yourConsumerSecretHere")
sig <- sign_oauth1.0(myapp, token = "yourTokenHere", token_secret = "yourTokenSecretHere")
homeTL <- GET("http://api.twitter.com/1.1/statuses/home_timeline.json", sig)

json1 <- content(homeTL)
json2 <- jsonlite::fromJSON(toJSON(json1))
json2
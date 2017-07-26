# install.packages("mailR")
library(mailR)

# 到下列網址設定：「登入和安全性」 >>> 「允許安全性較低的應用程式」
# https://myaccount.google.com/security?pli=1#connectedapps

send.mail(from = "@gmail.com",  # 寄
          to = "@gmail.com",  # 收
          subject = "標題測試",
          body = "
內文測試123
換行測試
abc
", # 換行必須靠齊左側，換行才不會產生空格
          encoding = "utf-8",  # 編碼
          #夾帶檔案
          # attach.files = c("FBCrawler.R"),
          smtp = list(host.name = "smtp.gmail.com", # gmail通用
                      port = 465, # gmail通用
                      user.name = "ACCOUNT@gmail.com", 
                      passwd = "PASSWORD", 
                      ssl = TRUE), # gmail通用
          authenticate = TRUE,
          send = TRUE)

# try gmailr
# try sendmailR

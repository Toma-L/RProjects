# R High Performance Programming =====

system.time(runif(1e8))

# install.packages("rbenchmark")
library(rbenchmark)

bench1 <- benchmark(runif(1e8), replications = 10)
bench1

within(bench1, {
        elapsed.mean <- elapsed/replications
        user.self.mean <- user.self/replications
        sys.self.mean <- sys.self/replications
})

benchmark(runif(1e8), replications = rep.int(1, 10))

# install.packages("microbenchmark")
library(microbenchmark)

microbenchmark(runif(1e8), times = 10)



# 6. 減少內存使用的簡單方法 =====

x <- runif(1e6)
print(object.size(x), units = "auto")
y <- list(x, x)
print(object.size(y), units = "auto")

# install.packages("pryr")
library(pryr)
object_size(x)
object_size(y) # 另一種衡量內存方式，結果大不相同


library(kappalab)

## a capacity
mu <- capacity(c(0,0.1,0.6,rep(0.9,4),1))
## the cdf of the Choquet integral at 0.7
cdf.Choquet(mu,0.7)

## the same but empirically
m <- 10000
ch <- numeric(m)
for (i in 1:m) {
     f <- runif(3) 
     ch[i] <- Choquet.integral(mu,f)
}
sum(ifelse(ch<=0.7,1,0))/m

require("MASS")


hp <- list('pow.ex.w'=log(10),'pow.ex.v'=log(5),
           'linear.a'=log(10), 'linear.i'=log(10), 'vv'=log(1))
n <- 40
input <- seq(0, 1, length.out=n)
idx <- sort(sample(1:40, 21))
X <- as.matrix(c[idx])
Sigma <- cov.linear(hyper=hp, Data=input) + 
  cov.pow.ex(hyper=hp, Data=input, gamma=2) + 
  diag(exp(hp$vv), n, n)
Y <- matrix(mvrnorm(n=1, mu=rep(0,n), Sigma=Sigma), ncol=1)*0.1+sin(input*6)
Y <- as.matrix(Y[idx])
x <- as.matrix(seq(0,1,by=0.03))
a <- gpr(Data=X, response=Y, Cov=c('linear','pow.ex'), gamma=2, trace=4)
plot(a)
b1 <- gppredict(a,x, noiseFreePred=F)
plot(b1)
b2 <- gppredict(a,x, noiseFreePred=T)
plot(b2)

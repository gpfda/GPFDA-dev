

# rm(list = ls())
library("GPFDA")
require("MASS")

library("useful")
library("rbenchmark")
library("Rcpp")
library("inline")


n <- 100 # training set size + test set size
nSamples <- 1
# hp <- list('linear.a'=log(10), 'pow.ex.v'=log(5),'pow.ex.w'=log(10),'vv'=log(1))
hp <- list('pow.ex.v'=log(2),'pow.ex.w'=log(10),'vv'=log(0.3))
h <- seq(0,1,len=n)
idx <- sort(sample(1:n,floor(n/2)))
X <- as.matrix(h[idx])
# Y <- (mvrnorm(n=n,mu=h-h,Sigma=(cov.linear(hp,h)+cov.pow.ex(hp,h)))[,1]) # *0.1+sin(h*6)
gamma <- 0.5 # 1 for SqExp
Sigma <- cov.pow.ex(hp,h, gamma = gamma)

# source("gp.functions5NEW.R")

source("CovMaternCpp_sq.R")
source("CovMaternCpp.R")

cc <- exp(hp$pow.ex.v)
A <- as.matrix(exp(hp$pow.ex.w))
nu <- 0.5
SigmaRcpp <- CovMaternCpp_sq(X=as.matrix(h), cc = cc, A = A, nu = nu)
SigmaRcpp_nm <- CovMaternCpp(X=as.matrix(h), Xnew=as.matrix(h[1:15]), cc = cc, A = A, nu = nu)

corner(Sigma)
corner(SigmaRcpp)
corner(SigmaRcpp_nm)

ts.plot(Sigma[1,])
ts.plot(SigmaRcpp[1,])
ts.plot(SigmaRcpp[1,1:15])
ts.plot(SigmaRcpp_nm[1,])

dim(SigmaRcpp_nm)

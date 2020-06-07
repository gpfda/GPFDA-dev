# rm(list = ls()) 

library("GPFDA")
require("MASS")

library("useful")
library("rbenchmark")
library("Rcpp")
library("inline")

source("DistMat_sq.R")
source("DistMat.R")
# source("DistMatLinear_sq.R")
# source("DistMatLinear.R")
source("gp.functions5NEW_1.R")


n <- 2000 # training set size + test set size
nSamples <- 1
hp <- list('pow.ex.v'=log(2),'pow.ex.w'=log(10^2),'vv'=log(0.3)) # for gam=0.5
gam <- 2
h <- seq(0,1,len=n)
idx <- sort(sample(1:n,floor(n/2)))
X <- as.matrix(h[idx])

Sigma <- cov.pow.ex(hyper = hp, Data = h, gamma = gam)
diag(Sigma) <- diag(Sigma) + exp(hp$vv)
ts.plot(Sigma[1,])


if(nSamples==1){
  Y <- matrix(mvrnorm(n=nSamples,mu=h-h,Sigma=Sigma), ncol=1)
}else{
  Y <- t(mvrnorm(n=nSamples,mu=h-h,Sigma=Sigma))
}
str(Y)
str(X)


# plot(h, Y)

Y <- as.matrix(Y[idx,])

str(Y)
str(X)

#############################################################################
##### Estimation

Data=X
response=Y
Cov=c('pow.ex')
trace=2
gamma=gam
nInitCandidates = 10
hyper=NULL
m=NULL
mean=0
itermax=100;reltol=8e-10;trace=0
useGradient=T
nu=NULL

unlist(hp)

t0 <- proc.time()
aNew <- gprNEW(Data=X, response=Y, Cov=c('pow.ex'),trace=2, useGradient = T,
               gamma=gam, nInitCandidates = 30)
newGPFDA_time <- proc.time() - t0

exp(unlist(hp)) # true
(newGPFDA_estimates <- round(exp(unlist(aNew$hyper)), 4))

SigmaHat <- cov.pow.ex(hyper = aNew$hyper, Data = h, gamma=gam)
diag(SigmaHat) <- diag(SigmaHat) + exp(aNew$hyper$vv)
ts.plot(SigmaHat[1,])
# ts.plot(Sigma[1,])





###############################################################
##### Prediction

plot(X, Y, type="p", pch=20)

x <- as.matrix(seq(0,1,by=0.01))

train <- aNew
Data.new <- x
# aNew$gamma <- gam/2 # different parametrisations in older version of GPFDA package
aNewGPFDA <- aNew
aNewGPFDA$gamma <- aNew$gamma/2 # different parametrisation
aNewGPFDA$hyper$pow.ex.w <- aNew$hyper$pow.ex.w*1.155  # different parametrisation
b <- GPFDA:::gppredict(train=aNewGPFDA, Data.new=x)
plot(x, b$pred.mean, ylim=range(Y), type="l")

aNew$gamma <- gam   # put back the new parametrisation
bNew <- gppredictNEW(train=aNew, Data.new=x, noiseFreePred = F)
lines(x, bNew$pred.mean)

GPFDA:::plot.gpr(b)
plot(bNew, ylim=range(Y))
# text(0.8, 3, "new results", cex=2)

max(abs(b$pred.mean - bNew$pred.mean))
max(abs(b$pred.sd - bNew$pred.sd))


plot(x, b$pred.sd, ylim=range(c(b$pred.sd, bNew$pred.sd)), cex=3)
lines(x, bNew$pred.sd, col="blue", type="p", pch=20)


covPowExp_OLD <- GPFDA::cov.pow.ex(hyper = aNewGPFDA$hyper, Data = x, gamma=gam/2)
plot(covPowExp_OLD[1,1:20], type="l", col="red")
covPowExp_NEW <- cov.pow.ex(hyper = aNew$hyper, Data = x, gamma=gam)
lines(covPowExp_NEW[1,1:20], col="blue")

log(covPowExp_OLD[1,1:20])/log(covPowExp_NEW[1,1:20])


bNew <- gppredictNEW(train=aNew, Data.new=x, noiseFreePred = T)
plot(bNew, ylim=range(Y))


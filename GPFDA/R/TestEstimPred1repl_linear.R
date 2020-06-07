# rm(list = ls())

library("GPFDA")
require("MASS")

library("useful")
library("rbenchmark")
library("Rcpp")
library("inline")

# source("DistMat_sq.R")
# source("DistMat.R")
source("DistMatLinear_sq.R")
source("DistMatLinear.R")
source("gp.functions5NEW_1.R")


n <- 3000 # training set size + test set size
nSamples <- 1
hp <- list('linear.a'=log(50), 'vv'=log(0.01))
h <- seq(0,1,len=n)
idx <- sort(sample(1:n,floor(n/2)))
X <- as.matrix(h[idx])

Sigma <- cov.linear(hyper = hp, Data = h)
diag(Sigma) <- diag(Sigma) + exp(hp$vv)
corner(Sigma)
# ts.plot(Sigma[1,])
# ts.plot(Sigma[1,1:10])
# 
# Sigma[100,100:110]
# ts.plot(Sigma[100,100:110])
ts.plot(diag(Sigma))

if(nSamples==1){
  Y <- matrix(mvrnorm(n=nSamples,mu=h-h,Sigma=Sigma), ncol=1)
}else{
  Y <- t(mvrnorm(n=nSamples,mu=h-h,Sigma=Sigma))
}
str(Y)
str(X)
# ts.plot(Y)
# plot(h, Y)

Y <- as.matrix(Y[idx,])

str(Y)
str(X)

#############################################################################
##### Estimation

Data=X
response=Y
Cov=c('linear')
trace=2
nInitCandidates = 10
hyper=NULL
m=NULL
mean=0
itermax=100;reltol=8e-10;trace=0


unlist(hp)
t0 <- proc.time()
aNew <- gprNEW(Data=X, response=Y, Cov=c('linear'),trace=2, useGradient=T, itermax=50,
               nInitCandidates = 30)
newGPFDA_time <- proc.time() - t0

exp(unlist(hp)) # true
(newGPFDA_estimates <- round(exp(unlist(aNew$hyper)), 4))



SigmaHat <- cov.linear(hyper = aNew$hyper, Data = h)
diag(SigmaHat) <- diag(SigmaHat) + exp(aNew$hyper$vv)
ts.plot(diag(SigmaHat))
# ts.plot(diag(Sigma))



#############################################################################
##### Prediction

plot(X, Y, type="p", pch=20)

x <- as.matrix(seq(0,1,by=0.01))

train <- aNew
Data.new <- x
bNew <- gppredictNEW(train=aNew, Data.new=x, noiseFreePred = F)
plot(bNew, ylim=range(Y))
# text(0.8, 3, "new results", cex=2)
bNew <- gppredictNEW(train=aNew, Data.new=x, noiseFreePred = T)
plot(bNew, ylim=range(Y))


plot(x, bNew$pred.mean, ylim=range(Y), type="l")
plot(x, bNew$pred.sd)



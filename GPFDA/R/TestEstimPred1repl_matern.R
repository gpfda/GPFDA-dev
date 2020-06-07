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

source("CovMaternCpp_sq.R")
source("CovMaternCpp.R")

n <- 2000 # training set size + test set size
nSamples <- 1
hp <- list('matern.v'=log(2),'matern.w'=log(10),'vv'=log(0.05))
nu <- 1.5
h <- seq(0,1,len=n)


# # ####################
# # checking if Mat with nu=0.5 is identical to Exp cov fun (pow.exp gam=1)
# hp <- list('matern.v'=log(2),'matern.w'=log(10^2))
# SigmaMatern <- cov.matern(hyper = hp, Data = h, nu = 0.5)
# corner(SigmaMatern)
# ts.plot(SigmaMatern[1,])
# SigmaExp <- cov.pow.ex(hyper = list('pow.ex.v'=log(2),
#                                     'pow.ex.w'=log(10^2),
#                                     'vv'=log(0.3)), Data = h, gamma = 1)
# max(abs(SigmaMatern - SigmaExp))
# 
# ####################
# # checking if Mat with large nu is close to  Squared-Exp cov fun (pow.exp gam=2)
# hp <- list('matern.v'=log(2),'matern.w'=log(10))
# SigmaMatern <- cov.matern(hyper = hp, Data = h, nu = 25)
# corner(SigmaMatern)
# ts.plot(SigmaMatern[1,])
# SigmaSqExp <- cov.pow.ex(hyper = list('pow.ex.v'=log(2),
#                                     'pow.ex.w'=log(100),
#                                     'vv'=log(0.3)), Data = h, gamma = 2)
# max(abs(SigmaMatern - SigmaSqExp))
# corner(SigmaSqExp)
# ts.plot(SigmaSqExp[1,])

# ####################
# # checking if Mat with nu=5/2 is close to Matern with nu=5/2+1e-8
# hp <- list('matern.v'=log(2),'matern.w'=log(10))
# SigmaMatern <- cov.matern(hyper = hp, Data = h, nu = 2.5)
# corner(SigmaMatern)
# ts.plot(SigmaMatern[1,])
# SigmaMaternClose <- cov.matern(hyper = hp, Data = h, nu = 2.5+1e-8)
# max(abs(SigmaMatern - SigmaMaternClose))
# corner(SigmaMaternClose)
# ts.plot(SigmaMaternClose[1,])



####################

idx <- sort(sample(1:n,floor(n/2)))
X <- as.matrix(h[idx])

Sigma <- cov.matern(hyper = hp, Data = h, nu = nu)
diag(Sigma) <- diag(Sigma) + exp(hp$vv)
ts.plot(Sigma[1,])
corner(Sigma[1,])


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
Cov=c('matern')
trace=2
nu=nu
nInitCandidates = 10
# hyper=NULL
m=NULL
mean=0
itermax=50;reltol=8e-10;trace=0
useGradient=T
gamma=NULL

unlist(hp)

t0 <- proc.time()
aNew <- gprNEW(Data=X, response=Y, Cov=c('matern'),useGradient=T, itermax=50,
               trace=2, nu=nu, nInitCandidates = 30)
newGPFDA_time <- proc.time() - t0

exp(unlist(hp)) # true
(newGPFDA_estimates <- round(exp(unlist(aNew$hyper)), 4))



SigmaHat <- cov.matern(hyper = aNew$hyper, Data = h, nu=nu)
diag(SigmaHat) <- diag(SigmaHat) + exp(aNew$hyper$vv)
ts.plot(SigmaHat[1,])
# ts.plot(Sigma[1,])



###############################################################
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



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


n1 <- 40 # training set size + test set size
n2 <- 40
(n <- n1*n2)
x1 <- seq(0,1,len=n1)
x2 <- seq(0,1,len=n2)

Cov <- c('matern')
# hp <- list('matern.v'=log(2),'matern.w'=c(log(25), log(20)),'vv'=log(0.1))
hp <- list('matern.v'=log(2),'matern.w'=c(log(10), log(15)),'vv'=log(0.5))
nu <- 1.5
# hp <- list('matern.v'=log(2),'matern.w'=c(log(25), log(20)),'vv'=log(0.5))
# nu <- 2.5


nSamples <- 1
X <- as.matrix(expand.grid(x1=x1, x2=x2))
X0 <- X

Sigma <- cov.matern(hyper = hp, Data = X0, nu = 0.5)
diag(Sigma) <- diag(Sigma) + exp(hp$vv)
ts.plot(Sigma[1,])
tail(eigen(Sigma)$values)


mu <- rep(0, n)

if(nSamples==1){
  Y <- matrix(mvrnorm(n=nSamples,mu=mu,Sigma=Sigma), ncol=1)
}else{
  Y <- t(mvrnorm(n=nSamples,mu=mu,Sigma=Sigma))
}
str(Y)
str(X)
Y0 <- Y
# plot(h, Y)

idx <- sort(sample(1:n,floor(n/2))) # take about half for training 
X <- X[idx,,drop=F]
Y <- as.matrix(Y[idx,])

str(Y)
str(X)

#############################################################################
##### Estimation

Data=X
response=Y

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
aNew <- gprNEW(Data=X, response=Y, Cov=c('matern'),trace=2, useGradient = F,
               itermax=itermax, nu=nu, nInitCandidates = 30)
newGPFDA_time <- proc.time() - t0

exp(unlist(hp)) # true
(newGPFDA_estimates <- round(exp(unlist(aNew$hyper)), 4))

SigmaHat <- cov.matern(hyper = aNew$hyper, Data = X, nu=nu)
diag(SigmaHat) <- diag(SigmaHat) + exp(aNew$hyper$vv)
ts.plot(SigmaHat[1,])
ts.plot(Sigma[1,])





###############################################################
##### Prediction


train <- aNew
Data.new <- X0
bNew <- gppredictNEW(train=aNew, Data.new=X0, noiseFreePred = T)

predMean <- bNew$pred.mean
predMeanMat <- matrix(bNew$pred.mean, nrow=n1, ncol=n2, byrow=F)

ObsMat <- matrix(Y0, nrow=n1, ncol=n2, byrow=F)


# packages useful for plotting eigensurfaces
library("akima")
library("fields")


set.panel()
par(mar=c(2.1,2.1,4.1,4.1))
par(oma=c( 0,0,0,0))
set.panel( 1,2)
myRange <- range(c(ObsMat, predMeanMat))
myRange <- myRange + 0.2*c(-1,1)

for(met in 1:2){
    if(met==1){
      SurfToPlot <- ObsMat
      CovFunName <- "observed"
    }
    if(met==2){
      SurfToPlot <- predMeanMat
      CovFunName <- "prediction"
    }
  
  akima.li <- interp(x=X0[,1], y=X0[,2], z=as.numeric(SurfToPlot),
                     nx=200, ny=200, linear = F)
  
  image(x = akima.li$x, y=akima.li$y, z=akima.li$z, col=tim.colors(),
        xlab=NA, ylab=NA, cex.lab=2.5, zlim=myRange)

  title(main = paste(CovFunName), font.main = 2, cex.main=2.5, line=1)

  fields::image.plot( legend.only=TRUE, zlim=myRange, legend.width=2.5, 
                      axis.args=list(cex.axis=2))

}
set.panel() # reset plotting device



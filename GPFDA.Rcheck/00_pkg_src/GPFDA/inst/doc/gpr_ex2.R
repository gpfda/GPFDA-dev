## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
        echo = TRUE, results = 'hold', warning=F, cache=F, 
  #dev = 'pdf', 
  message=F, 
  fig.width=5, fig.height=5,
  tidy.opts=list(width.cutoff=75), tidy=FALSE
)
options(scipen = 1, digits = 4)

## ----setup---------------------------------------------------------------
library(GPFDA)
require(MASS)
# packages required for visualisation:
require(interp)
require(fields)

## ------------------------------------------------------------------------
set.seed(123)
nrep <- 10
n1 <- 30
n2 <- 30
n <- n1*n2
input1 <- seq(0,1,len=n1)
input2 <- seq(0,1,len=n2)
input <- as.matrix(expand.grid(input1=input1, input2=input2))
hp <- list('matern.v'=log(2),'matern.w'=c(log(20), log(25)),'vv'=log(0.2))
nu <- 1.5
Sigma <- cov.matern(hyper = hp, input = input, nu = nu) + diag(exp(hp$vv), n)
Y <- t(mvrnorm(n=nrep, mu=rep(0,n), Sigma=Sigma))

## ------------------------------------------------------------------------
idx <- expand.grid(1:n1, 1:n2)
n1test <- floor(n1*0.8)
n2test <- floor(n2*0.8)
idx1 <- sort(sample(1:n1, n1test))
idx2 <- sort(sample(1:n2, n2test))
whichTest <- idx[,1]%in%idx1 & idx[,2]%in%idx2

inputTest <- input[whichTest, ]
Ytest <- Y[whichTest, ]
inputTrain <- input[!whichTest, ]
Ytrain <- Y[!whichTest, ]

## ------------------------------------------------------------------------
fit <- gpr(input=inputTrain, response=Ytrain, Cov='matern', trace=4, useGradient=T,
            iter.max=50, nu=nu, nInitCandidates=50)

## ------------------------------------------------------------------------
sapply(fit$hyper, exp)

## ------------------------------------------------------------------------
pred <- gprPredict(train=fit, input.new=inputTest, noiseFreePred=T)

## ---- fig.width=6, fig.height=5------------------------------------------
zlim <- range(c(pred$pred.mean, Ytest))

plotImage(response = Ytest, input = inputTest, realisation = 1, 
            n1 = n1test, n2 = n2test,
            zlim = zlim, main = "observed")

plotImage(response = pred$pred.mean, input = inputTest, realisation = 1, 
            n1 = n1test, n2 = n2test,
            zlim = zlim, main = "prediction")


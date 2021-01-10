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

## ------------------------------------------------------------------------
set.seed(123)
nrep <- 30
n <- 15
input <- seq(0, 1, length.out=n)
hp <- list('linear.a'=log(40), 'linear.i'=log(10),
           'pow.ex.v'=log(5), 'pow.ex.w'=log(15),
           'vv'=log(0.3))
Sigma <- cov.linear(hyper=hp, input=input) + 
  cov.pow.ex(hyper=hp, input=input, gamma=2) + 
  diag(exp(hp$vv), n, n)
Y <- t(mvrnorm(n=nrep, mu=rep(0,n), Sigma=Sigma))

## ------------------------------------------------------------------------
set.seed(111)
fitNoGrad <- gpr(input=input, response=Y, Cov=c('linear','pow.ex'), gamma=2, 
               trace=4, nInitCandidates = 1, useGradient = F)

## ------------------------------------------------------------------------
set.seed(111)
fit <- gpr(input=input, response=Y, Cov=c('linear','pow.ex'), gamma=2, 
         trace=4, nInitCandidates = 1, useGradient = T)

## ------------------------------------------------------------------------
sapply(fit$hyper, exp)

## ------------------------------------------------------------------------
plot(fit, realisation=10)

## ------------------------------------------------------------------------
input.new <- seq(0, 1, length.out = 1000)
pred1 <- gprPredict(train=fit, input.new=input.new, noiseFreePred=T)
plot(pred1, realisation=10)

## ------------------------------------------------------------------------
pred2 <- gprPredict(train=fit, input.new=input.new, noiseFreePred=F)
plot(pred2, realisation=10)


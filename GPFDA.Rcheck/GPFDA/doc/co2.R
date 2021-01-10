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

## ------------------------------------------------------------------------
library(GPFDA)
data(co2, package = "GPFDA")

## ------------------------------------------------------------------------
y <- c(t(data.matrix(co2[1:5, 1:12])))
n <- length(y)
input <- 5*(1:n)/n
whichNeg <- (y<0) # removing negative values
y <- y[!whichNeg]
input <- input[!whichNeg]

## ------------------------------------------------------------------------
input.new <- 5*seq(1/n, 1, length.out=1000) # input test for prediction

## ---- warnings=F, results=F, fig.width=8, fig.height=5-------------------
set.seed(789)
fit1 <- gpr(input=input, response=y, Cov='pow.ex', gamma=1, meanModel='t',
                      nInitCandidates=50, trace=2)

## ------------------------------------------------------------------------
sapply(fit1$hyper, exp)

## ---- warnings=F, results=F, fig.width=8, fig.height=5-------------------
plot(fit1, main="exponential kernel - fitted model")

## ---- warnings=F, results=F, fig.width=8, fig.height=5-------------------
pred1 <- gprPredict(train = fit1, input.new=input.new)
plot(pred1, main="exponential kernel - predictions")

## ------------------------------------------------------------------------
cov.custom <- function(hyper, input, input.new=NULL){
  hyper <- lapply(hyper, exp)
  if(is.null(input.new)){
    input.new <- as.matrix(input)
  }else{
    input.new <- as.matrix(input.new)
  }
  input <- as.matrix(input)
  A1 <- distMat(input=input, inputNew=input.new, A=as.matrix(hyper$custom.w), 
                        power=2)
  sept <- outer(input[,1], input.new[,1], "-")
  A2 <- hyper$custom.u*(sin(pi*sept))^2
  customCov <- hyper$custom.v*exp(-A1-A2)
  return(customCov)
}

## ------------------------------------------------------------------------
Dloglik.custom.w <- function(hyper, input, AlphaQ){
  A1 <- distMatSq(input=input, A=as.matrix(exp(hyper$custom.w)), power=2)
  Dcov <- - cov.custom(hyper, input)*A1
  out <- 0.5*sum(diag(AlphaQ%*%Dcov))
  return(out)
}

Dloglik.custom.u <- function(hyper, input, AlphaQ){
  sept <- outer(input[,1], input[,1], "-")
  A2 <- exp(hyper$custom.u)*(sin(pi*sept))^2
  Dcov <- - cov.custom(hyper, input)*A2
  out <- 0.5*sum(diag(AlphaQ%*%Dcov))
  return(out)
}

Dloglik.custom.v <- function(hyper, input, AlphaQ){
  Dcov <- cov.custom(hyper, input)
  out <- 0.5*sum(diag(AlphaQ%*%Dcov))
  return(out)
}

## ------------------------------------------------------------------------
D2custom.w <- function(hyper, input, inv.Q, Alpha.Q){
  Cov <- cov.custom(hyper, input)
  A1 <- distMatSq(input=input, A=as.matrix(exp(hyper$custom.w)), power=2)
  D1cov <- -Cov*A1
  D2cov <- -(D1cov + Cov)*A1
  D2c.w <- D2(D1cov, D2cov, inv.Q, Alpha.Q)
  return(D2c.w)
}

D2custom.u <- function(hyper, input, inv.Q, Alpha.Q){
  Cov <- cov.custom(hyper, input)
  sept <- outer(input[,1], input[,1], "-")
  A2 <- exp(hyper$custom.u)*(sin(pi*sept))^2
  D1cov <- - Cov*A2
  D2cov <- -(D1cov + Cov)*A2
  D2c.u <- D2(D1cov, D2cov, inv.Q, Alpha.Q)
  return(D2c.u)
}

D2custom.v <- function(hyper, input, inv.Q, Alpha.Q){
  out <- cov.custom(hyper, input)
  return(out)
}

## ------------------------------------------------------------------------
diag.custom <- function(hyper, input){
  Qstar <- rep(exp(hyper$custom.v), nrow(input))
  return(Qstar)
}

## ---- message=F, results=F, fig.width=8, fig.height=5--------------------
fit2 <- gpr(input=input, response=y, Cov='custom',
            NewHyper=c('custom.w','custom.u','custom.v'), gamma=2, meanModel='t',
            nInitCandidates=50, trace=2, useGradient = F)

## ------------------------------------------------------------------------
sapply(fit2$hyper, exp)

## ---- message=F, fig.width=8, fig.height=5-------------------------------
plot(fit2, main="customised kernel - fitted model")

## ---- message=F, fig.width=8, fig.height=5-------------------------------
pred2 <- gprPredict(train = fit2, input.new=input.new)
plot(pred2, main="customised kernel - predictions")


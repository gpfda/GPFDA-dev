## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
        echo = TRUE, results = 'hold', warning=F, cache=F, #dev = 'pdf', 
  message=F, 
  fig.width=5, fig.height=5,
  tidy.opts=list(width.cutoff=75), tidy=FALSE
)
options(scipen = 1, digits = 4)

## ----setup---------------------------------------------------------------
library(GPFDA)
require(MASS)

## ------------------------------------------------------------------------
set.seed(100)

M <- 20
n <- 50
p <- 2  # number of scalar covariates

hp <- list('pow.ex.v'=log(10), 'pow.ex.w'=log(1),'vv'=log(1))

## Training data: M realisations -----------------
 
tt <- seq(-4,4,len=n)
b <- sin((0.5*tt)^3)

scalar_train <- matrix(NA, M, p)
t_train <- matrix(NA, M, n)
x_train <- matrix(NA, M, n)
response_train <- matrix(NA, M, n)
for(i in 1:M){
  u0 <- rnorm(1)
  u1 <- rnorm(1,10,5)
  x <- exp(tt) + rnorm(n, 0, 0.1)
  Sigma <- cov.pow.ex(hp, cbind(x))
  diag(Sigma) <- diag(Sigma) + exp(hp$vv)
  y <- u0+u1*b + mvrnorm(n=1, mu=rep(0,n), Sigma=Sigma)
  scalar_train[i,] <- c(u0,u1)
  t_train[i,] <- tt
  x_train[i,] <- x
  response_train[i,] <- y
}

## Test data (M+1)-th realisation ------------------
n_new <- 100
t_new <- seq(-4,4,len=n_new)
b_new <- sin((0.5*t_new)^3)
u0_new <- rnorm(1)
u1_new <- rnorm(1,10,5)
scalar_new <- cbind(u0_new, u1_new)
x_new <- exp(t_new) + rnorm(n_new, 0, 0.1)
Sigma_new <- cov.pow.ex(hp,cbind(x_new))
diag(Sigma_new) <- diag(Sigma_new) + exp(hp$vv)
response_new <- u0_new + u1_new*b_new + mvrnorm(n=1, mu=rep(0,n_new), 
                                                Sigma=Sigma_new)

## ---- include=F, eval=F--------------------------------------------------
#  dataExampleGPFR <- list(tt=tt,
#                          response_train=response_train,
#                          x_train=x_train,
#                          scalar_train=scalar_train,
#                          t_new=t_new,
#                          response_new=response_new,
#                          x_new=x_new,
#                          scalar_new=scalar_new)
#  save(dataExampleGPFR, file = "data/dataExampleGPFR.rda")

## ---- results=F----------------------------------------------------------
a1 <- gpfr(response = response_train, time = tt, lReg = scalar_train,
           fReg = NULL, gpReg = list(x_train),
           fyList = list(nbasis = 23, lambda = 0.0001),
           fbetaList_l = list(list(lambda = 0.0001, nbasi = 23)),
           Cov = 'pow.ex',  fitting = T)

## ------------------------------------------------------------------------
unlist(lapply(a1$hyper,exp))

## ------------------------------------------------------------------------
plot(a1, type='raw')

## ------------------------------------------------------------------------
plot(a1, type='raw', realisations = 1:3)

## ------------------------------------------------------------------------
plot(a1, type = 'meanFunction', realisations = 1:3)

## ------------------------------------------------------------------------
plot(a1, type = 'fitted', realisations = 1:3)

## ---- results=F----------------------------------------------------------
b1 <- gpfrPredict(a1, TestData = as.matrix(x_new), NewTime = t_new,
               lReg = scalar_new, fReg = NULL,
               gpReg = list('response' = response_new,
                            'input' = x_new, 'time' = t_new))

plot(b1, type = 'prediction')
lines(t_new, response_new, type = 'b', col = 4, pch = 19, cex = 0.6, lty = 3, lwd = 2)

## ---- results=F----------------------------------------------------------
b2 <- gpfrPredict(a1, TestData = as.matrix(x_new), NewTime = t_new,
               lReg = scalar_new, fReg = NULL,
               gpReg = list('response' = response_new[1:20],
                          'input' = x_new[1:20],
                          'time' = t_new[1:20]))

plot(b2, type = 'prediction')
lines(t_new, response_new, type = 'b', col = 4, pch = 19, cex = 0.6, lty = 3, lwd = 2)

## ---- results=F----------------------------------------------------------
b3 <- gpfrPredict(a1, TestData = as.matrix(x_new), NewTime = t_new,
               lReg = scalar_new, fReg = NULL, gpReg = NULL)

plot(b3, type = 'prediction')
lines(t_new, response_new, type='b', col = 4, pch = 19, cex = 0.6, lty = 3, lwd = 2)


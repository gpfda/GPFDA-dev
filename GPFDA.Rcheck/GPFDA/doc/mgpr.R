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
set.seed(123)
nrep <- 30
n1 <- 250
n2 <- 250
n3 <- 250
N <- 3
n <- n1+n2+n3
input1 <- sapply(1:n1, function(x) (x - min(1:n1))/max(1:n1 - min(1:n1)))
input2 <- input1
input3 <- input1

# storing input vectors in a list
Data <- list()
Data$input <- list(input1, input2, input3)

# true hyperparameter values
nu0s <- c(6, 4, 2)
nu1s <- c(0.1, 0.05, 0.01)
a0s <- c(500, 500, 500)
a1s <- c(100, 100, 100)
sigm <- 0.05
hp <- c(nu0s, log(nu1s), log(a0s), log(a1s), log(sigm))

# Calculate covariance matrix
Psi <- mgpCovMat(Data=Data, hp=hp)

## ------------------------------------------------------------------------
ns <- sapply(Data$input, length)
idx <- c(unlist(sapply(1:N, function(i) rep(i, ns[i])))) 

## ------------------------------------------------------------------------
# Plotting an auto-covariance function
plotmgpCovFun(type="Cov", output=1, outputp=1, Data=Data, hp=hp, idx=idx)
# Plotting a cross-covariance function
plotmgpCovFun(type="Cov", output=1, outputp=2, Data=Data, hp=hp, idx=idx)

## ------------------------------------------------------------------------
# Plotting an auto-correlation function
plotmgpCovFun(type="Cor", output=1, outputp=1, Data=Data, hp=hp, idx=idx)
# Plotting a cross-correlation function
plotmgpCovFun(type="Cor", output=1, outputp=2, Data=Data, hp=hp, idx=idx)

## ------------------------------------------------------------------------
mu <- c( 5*input1, 10*input2, -3*input3)
Y <- t(mvrnorm(n=nrep, mu=mu, Sigma=Psi))
response <- list()
for(j in 1:N){
  response[[j]] <- Y[idx==j,,drop=F]
}
# storing the response in the list
Data$response <- response

## ------------------------------------------------------------------------
res <- mgpr(Data=Data, m=200, meanModel = 't')

## ------------------------------------------------------------------------
n_star <- 60*3
input1star <- seq(min(input1), max(input1), length.out = n_star/N)
input2star <- seq(min(input2), max(input2), length.out = n_star/N)
input3star <- seq(min(input3), max(input3), length.out = n_star/N)
Data.new <- list()
Data.new$input <- list(input1star, input2star, input3star)

## ------------------------------------------------------------------------
obsSet <- list()
obsSet[[1]] <- c(5, 10, 23, 50, 80, 200)
obsSet[[2]] <- c(5, 10, 23, 200)
obsSet[[3]] <- c(5, 10, 23, 200)

Data.obs <- list()
Data.obs$input[[1]] <- Data$input[[1]][obsSet[[1]]]
Data.obs$input[[2]] <- Data$input[[2]][obsSet[[2]]]
Data.obs$input[[3]] <- Data$input[[3]][obsSet[[3]]]
Data.obs$response[[1]] <- Data$response[[1]][obsSet[[1]], ]
Data.obs$response[[2]] <- Data$response[[2]][obsSet[[2]], ]
Data.obs$response[[3]] <- Data$response[[3]][obsSet[[3]], ]

## ------------------------------------------------------------------------
# Calculate predictions for the test set given some observations
predCGP <- mgprPredict(train=res, Data.obs=Data.obs, Data.new=Data.new)
str(predCGP)

## ---- fig.width=9, fig.height=4------------------------------------------
plot(res, Data.obs=Data.obs, Data.new=Data.new, realisation=5)


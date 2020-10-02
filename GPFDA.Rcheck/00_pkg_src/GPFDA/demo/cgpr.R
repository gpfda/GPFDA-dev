require("MASS")

set.seed(123)

# settings
nSamples <- 30

Q <- 1
n1 <- 250
n2 <- 250
n3 <- 250
N <- 3
n <- n1+n2+n3
input1 <- sapply(1:n1, function(x) (x - min(1:n1))/max(1:n1 - min(1:n1)))
input2 <- input1
input3 <- input1

Data <- list()
Data$input <- list(input1, input2, input3)

ns <- sapply(Data$input, length)
idx <- c(unlist(sapply(1:N, function(i) rep(i, ns[i]))))

# true hyperparameter values
nu0s <- c(6, 4, 2)
nu1s <- c(0.1, 0.05, 0.01)
a0s <- c(500, 500, 500)
a1s <- c(100, 100, 100)
sigm <- 0.05
hp <- c(nu0s, log(nu1s), log(a0s), log(a1s), log(sigm))

# Calculate covariance matrix
Psi <- CGPCovMat(Data=Data, hp=hp)

# Plot auto-covariance functions
plot.CGPCovFun(type="Cov", output=1, outputp=1, Data=Data, hp=hp)
plot.CGPCovFun(type="Cov", output=2, outputp=2, Data=Data, hp=hp)
plot.CGPCovFun(type="Cov", output=3, outputp=3, Data=Data, hp=hp)
# Plot cross-covariance functions
plot.CGPCovFun(type="Cov", output=1, outputp=2, Data=Data, hp=hp)
plot.CGPCovFun(type="Cov", output=1, outputp=3, Data=Data, hp=hp)
plot.CGPCovFun(type="Cov", output=3, outputp=1, Data=Data, hp=hp)

# Plot auto-correlation functions
plot.CGPCovFun(type="Cor", output=1, outputp=1, Data=Data, hp=hp)
plot.CGPCovFun(type="Cor", output=2, outputp=2, Data=Data, hp=hp)
plot.CGPCovFun(type="Cor", output=3, outputp=3, Data=Data, hp=hp)
# Plot cross-correlation functions
plot.CGPCovFun(type="Cor", output=1, outputp=2, Data=Data, hp=hp)
plot.CGPCovFun(type="Cor", output=1, outputp=3, Data=Data, hp=hp)
plot.CGPCovFun(type="Cor", output=2, outputp=3, Data=Data, hp=hp)


# Simulate data
mu <- c( 5*input1, 10*input2, -3*input3)
Y <- t(mvrnorm(n=nSamples, mu=mu, Sigma=Psi))
response <- list()
for(j in 1:N){
  response[[j]] <- Y[idx==j,,drop=F]
}
Data$response <- response


# CGPR estimation 
system.time(res <- CGPR(Data=Data, m=200, meanModel = 't'))


# Test set
n_star <- 60*3
input1star <- seq(min(input1), max(input1), length.out = n_star/N)
input2star <- seq(min(input2), max(input2), length.out = n_star/N)
input3star <- seq(min(input3), max(input3), length.out = n_star/N)
Data.new <- list()
Data.new$input <- list(input1star, input2star, input3star)

# Subset of training set
trainSet <- list()
trainSet[[1]] <- c(5, 10, 23, 50, 80, 200)
trainSet[[2]] <- c(5, 10, 23, 200)
trainSet[[3]] <- c(5, 10, 23, 200)

Data.train <- list()
Data.train$input[[1]] <- Data$input[[1]][trainSet[[1]]]
Data.train$input[[2]] <- Data$input[[2]][trainSet[[2]]]
Data.train$input[[3]] <- Data$input[[3]][trainSet[[3]]]
Data.train$response[[1]] <- Data$response[[1]][trainSet[[1]], ]
Data.train$response[[2]] <- Data$response[[2]][trainSet[[2]], ]
Data.train$response[[3]] <- Data$response[[3]][trainSet[[3]], ]

# Calculate predictions for the test set given training set observations
# predCGP <- CGPprediction(train=res,
#                          Data.train=Data.train,
#                          Data.new=Data.new)
# str(predCGP)

# Plot prediction of the 5-th realisation
plot.CGPprediction(train=res,
                   Data.train=Data.train,
                   Data.new=Data.new, i=5,
                   cex=2, cex.lab=2, cex.axis=2)


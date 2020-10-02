require("MASS")
library("akima")
library("fields")


# Simulate data
n1 <- 30
n2 <- 30
n <- n1*n2
x1 <- seq(0,1,len=n1)
x2 <- seq(0,1,len=n2)

hp <- list('matern.v'=log(2),'matern.w'=c(log(20), log(25)),'vv'=log(0.2))
nu <- 1.5

nSamples <- 10
X <- as.matrix(expand.grid(x1=x1, x2=x2))
Sigma <- cov.matern(hyper = hp, Data = X, nu = nu) + diag(exp(hp$vv), n)
mu <- rep(0,n)
set.seed(1234)
Y <- t(mvrnorm(n=nSamples, mu=mu, Sigma=Sigma))

# Split the dataset into train and test sets
idxAll <- expand.grid(1:n1,1:n2)
n1test <- floor(n1*0.8)
n2test <- floor(n2*0.8)

idx1 <- sort(sample(1:n1, n1test))
idx2 <- sort(sample(1:n2, n2test))
whichTest <- idxAll[,1]%in%idx1 & idxAll[,2]%in%idx2

Xtest <- X[whichTest,,drop=F]
Ytest <- Y[whichTest,,drop=F]
Xtrain <- X[!whichTest,,drop=F]
Ytrain <- Y[!whichTest,,drop=F]


##### Estimation
aNew <- gpr(Data=Xtrain, response=Ytrain, Cov='matern', trace=2, useGradient=T,
            itermax=50, nu=nu, nInitCandidates=50)

##### Prediction
bNew <- gppredict(train=aNew, Data.new=Xtest, noiseFreePred=T)

#### Plot predictions for the i-th realisation
i <- 2
predMeanMat <- matrix(bNew$pred.mean[,i], nrow=n1test, ncol=n2test, byrow=F)
YtrainMat <- matrix(Ytest[,i], nrow=n1test, ncol=n2test, byrow=F)

# packages useful for plotting eigensurfaces
library("akima")
library("fields")


opar <- par(no.readonly = TRUE)
par(mfrow=c(1,2), mar=c(2.1,2.1,3.1,4.1), oma=c( 0,1,0,2))
myRange <- range(c(YtrainMat, predMeanMat)) + 0.2*c(-1,1)
for(met in 1:2){
  if(met==1){
    SurfToPlot <- YtrainMat
    CovFunName <- "observed"
  }
  if(met==2){
    SurfToPlot <- predMeanMat
    CovFunName <- "prediction"
  }
  
  akima.li <- interp(x=Xtest[,1], y=Xtest[,2], z=as.numeric(SurfToPlot),
                     nx=200, ny=200, linear = F)
  
  image(x = akima.li$x, y=akima.li$y, z=akima.li$z, col=tim.colors(),
        xlab=NA, ylab=NA, cex.lab=2.5, zlim=myRange)
  
  title(main = paste(CovFunName), font.main=2, cex.main=2, line=1)
  
  fields::image.plot( legend.only=TRUE, zlim=myRange, legend.width=4,
                      axis.args=list(cex.axis=1.5))
  
}
par(opar)

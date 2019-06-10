\name{gpfrpred}
\alias{gpfrpred}
\title{Prediction of the Gaussian Process using functional regression}
\description{
Predict the new points in Gaussian Process using the training results}

\usage{
gpfrpred(object,TestData,NewTime=NULL,lReg=NULL,fReg=NULL,gpReg=NULL,
        GP_predict=TRUE)
}
\arguments{
  \item{object}{
The result from training with class`gpfda'. If missing, function stops running.
}
  \item{TestData}{
The test data. Must be matrix or fd object. 
}
  \item{NewTime}{
New time for test data. If NULL, default setting will be applied.
}
\item{lReg}{
The test scale data for functional regression with scale covariates.
}
  \item{fReg}{
The test functional data for functional regression with functional covariates.
}
  \item{gpReg}{
List of three items. The names of the items must be `response', `input', `time'. For type I prediction, `response' is the observed response for a new batch, `input' is the observed functional covariates for a new batch, `time' is the observation time for the previous two. If NULL, type II prediction will be carried out.
}
  \item{GP_predict}{
Logical. If true, GP prediction is carried out, otherwise functional prediction is carried out. Default to be True.
}
}
\details{
Two types of prediction are supplied. Type one is the new batch has a few observations, type two is the new batch has no observations.}
\value{
\item{ypred}{matrix of predicted value with confidence interval. First column is the fitted value, second and third are the confidence interval.}
\item{ypred.mean}{The mean value of the prediction.}
\item{ypred.sd}{The standard deviation of the prediction.}
\item{time}{time of test data}
\item{object}{all items trained from gpfr if exists}
}
\references{
  Shi, J Q., and Choi, T. (2011), \emph{Gaussian Process Regression Analysis for Functional Data}, Springer, New York.
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}
\seealso{
  \code{\link{gpr}}
}
\examples{
library(GPFDA)

# code from: demo('gpfr') 

traindata <- vector('list',20)
for(i in 1:20) traindata[[i]]=i
n <- 50
traindata <- lapply(traindata,function(i) {
  x <- seq(-3,3,len=n)
  y <- sin(x^2)-x+0.2*rnorm(n,0,3)
  x1 <- 0.5*x^3+exp(x)+rnorm(n,0,3)
  x2 <- cos(x^3)+0.2*rnorm(n,0,3)
  mat <- cbind(x,x1,x2,y)
  colnames(mat) <- c('time','x1','x2','y')
  scale <- t(c(2*(mean(y)>0.25)-1,(var(y)>3.6)*2-1,(sd(y)-sd(x)>1.4)*2-1))
  i <- list(mat,scale)
})

n <- 800 #test input
x <- seq(-3,3,len=n)
y <- sin(x^2)-x+0.2*rnorm(n,0,3)
x1 <- 0.5*x^3+exp(x)+rnorm(n,0,3)
x2 <- cos(x^3)+0.2*rnorm(n,0,3)
mat <- cbind(x,x1,x2,y)
colnames(mat) <- c('time','x1','x2','y')
scale <- t(c(2*(mean(y)>0.25)-1,(var(y)>3.6)*2-1,(sd(y)-sd(x)>1.4)*2-1))
# testdata[[1]]=vector('list',3)
n <- 100 # test new points
xt <- seq(1,3,len=n)
yt <- sin(xt^2)-xt+0.2*rnorm(n,0,3)
xt1 <- 0.5*xt^3+exp(xt)+rnorm(n,0,3)
xt2 <- cos(xt^3)+0.2*rnorm(n,0,3)
mat_t <- cbind(xt,xt1,xt2)
colnames(mat_t) <- c('time','xt1','xt2')
td <- list(mat,scale,mat_t)



lx=do.call('rbind',lapply(traindata,function(i)i[[2]]))
fx1=do.call('rbind',lapply(traindata,function(i)i[[1]][,2]))
fx2=do.call('rbind',lapply(traindata,function(i)i[[1]][,3]))
fy1=do.call('rbind',lapply(traindata,function(i)i[[1]][,4]))
time_old=traindata[[1]][[1]][,1]

pfx=td[[1]][,c(2,3)]
pfy=td[[1]][,4]
ptime=td[[1]][,1]
time_new=td[[3]][,1]
tfx=td[[3]][,c(2,3)]
tx=td[[2]]

## comment out because running time is a bit long
# system.time(a1<-gpfr(response=(fy1),lReg=lx,fReg=NULL,gpReg=list(fx1,fx2),
# fyList=list(nbasis=23,lambda=0.1),fbetaList_l=list(list(lambda=100,
# nbasi=17)),hyper=NULL,Cov=c('pow.ex','linear'),fitting=TRUE,
# time=seq(-3,3,len=50),rPreIdx=TRUE,concurrent=TRUE))

# type I prediction
# system.time(b1<-gpfrpred(a1,TestData=(tfx),NewTime=time_new,lReg=tx,
# fReg=NULL,gpReg=list('response'=(pfy),'input'=(pfx),'time'=ptime)))

# type II prediction
# system.time(b2<-gpfrpred(a1,TestData=(tfx),NewTime=time_new,lReg=tx,
# fReg=NULL,gpReg=NULL))


}

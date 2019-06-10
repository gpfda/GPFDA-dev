\name{gppredict}
\alias{gppredict}
\title{Prediction of the Gaussian Process}
\description{
Predict the new points in Gaussian Process using the training results or manual input}

\usage{
gppredict(train=NULL,Data.new=NULL,hyper=NULL, Data=NULL, Y=NULL, 
      Cov=NULL,gamma=NULL,lrm=NULL,mean=0)
}
\arguments{
  \item{train}{
The result from training which is a 'gpr' object. Default to be NULL. If NULL, do training based on the other given arguments; if TURE, other arguments (except for Data.new) will replaced by NULL; if FALSE, only do prediction based on the other given arguments.}
  \item{Data.new}{
The test data. Must be a vector or a matrix.
}
  \item{hyper}{
Hyper-parameter estimated from training. Can use manual input. Default to be NULL.
}
  \item{Data}{
The data from training. Must be a vector or a matrix. Default to be NULL. 
}
  \item{Y}{
The response from training. Must be a vector or a matrix. Default to be NULL. 
}
  \item{Cov}{
Names of covariance functions used. Default to be NULL. 
}
  \item{gamma}{
Parameter used in power exponential covariance function. Default to be NULL. 
}
  \item{lrm}{
The linear trend from learning. Default to be lrm. If lrm exists from learning, NULL will be replaced.
}
  \item{mean}{
Is the mean taken out when analysis? Default to be 0, which assumes the mean is zero. if assume mean is a constant, mean=1; if assume mean is a linear trend, mean='t'.
}
}
\details{
Use the result from training to predict the value for new points. 
}
\value{
\item{CovFun}{Covariance function type}
\item{fitted}{Fitted value of training data}
\item{fitted.sd}{Standard deviation of the fitted value of training data}
\item{gamma}{Parameter used in powered exponential covariance function}
\item{hyper}{Hyper-parameter estimated from training data}
\item{I}{Variance of the estimated hyper-parameters}
\item{pred.mean}{Estimated prediction mean}
\item{pred.sd}{Estimated prediction variance}
\item{train.x}{Training covariates}
\item{train.y}{Training response, may be transformed, for prediction use only}
\item{train.yOri}{Original training response}
}
\references{
  Shi, J Q., and Choi, T. (2011), \emph{Gaussian Process Regression Analysis for Functional Data}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}

\seealso{
  \code{\link{gpr}}
}

\examples{
library(GPFDA)
library(MASS) ## used to generate data
hp <- list('pow.ex.w'=log(10),'linear.a'=log(10),'pow.ex.v'=log(5),
      'vv'=log(1))
c <- seq(0,1,len=40)
idx <- sort(sample(1:40,21))
X <- as.matrix(c[idx])
Y <- (mvrnorm(n=40,mu=c-c,Sigma=(cov.linear(hp,c)+cov.pow.ex(hp,c)))[,1])+
      sin(c*6)
Y <- as.matrix(Y[idx])
x <- as.matrix(seq(0,1,by=0.03))
a <- gpr(X,Y,c('linear','pow.ex'))
b <- gppredict(a,x)
}


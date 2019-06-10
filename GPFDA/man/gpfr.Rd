\name{gpfr}
\alias{gpfr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gaussian Process in functional data.
}
\description{
Use functional regression to be the mean structure and Gaussian Process to be the covariance structure. 
}

\usage{
gpfr(response,lReg=NULL,fReg=NULL,fyList=NULL,fbetaList_l=NULL,
      fxList=NULL,fbetaList=NULL,concurrent=TRUE,fbetaList_f=NULL,
      gpReg=NULL,hyper=NULL,Cov=c('pow.ex','linear'),gamma=1,
      time=NULL,NewHyper=NULL,accuracy=c('high','normal','low'),
      trace.iter=5,fitting=FALSE,rPreIdx=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{response}{
The training response. can be an fd object or a matrix with nrow samples, ncol time points
}
  \item{lReg}{
The input variable for functional linear regression with scale covariates. Expected to be a matrix with nrow samples.
  }
    \item{fReg}{
The input variable for functional linear regression with functional covariates. Expected to be a matrix with nrow samples, or an fd object, or a list of matrices or fd objects.
  }
    \item{fyList}{
The list to control the smoothing of response. See details for more info.
  }
    \item{fbetaList_l}{
The list to control the smoothing of beta for functional regression with scale covariates. See details for more info.
}
    \item{fxList}{
The list to control the smoothing of functional covariates for functional regression with functional covariates. See details for more info.
  }
    \item{fbetaList_f}{
The list to control the smoothing of beta for functional regression with functional covariates. See details for more info.
  }  
    \item{fbetaList}{
The list to control the smoothing of functional covariates for functional regression with functional covariates and scale response. Not available for now.
  }
    \item{concurrent}{
Logical. If True concurrent functional regression will be carried out, otherwise the full functional regression will be carried out.
}
    \item{gpReg}{
Data for training Gaussian Process. Expecting matrix, fd object, list of matrices or list of fd objects.
}  
  \item{Cov}{
Kernel function or covariance function type(s).
}
  \item{hyper}{
Hyper parameter initial value. Default to be NULL.
}
  \item{NewHyper}{
Vector of the names of the new hyper parameters from customized kernel function.
}
  \item{gamma}{
Power parameter used in powered exponential.
}
  \item{time}{
Time used in globle setting for functional objects.
  }
  \item{accuracy}{
Optimization accuracy. Default to be high. 
}
  \item{trace.iter}{
Print the processing of iterations of optimization. 
}
  \item{fitting}{
Is fitting required or not. Default to be F.
}
  \item{rPreIdx}{
  Logical. If True, do random selected index for pre-optimization, otherwise use fixed index.
  }
}
\details{
fyList is a list with items: `time': a sequence of time points default to be 100 points from 0 to 1; `nbasis': number of basis functions used in smoothing, default to be less or equal to 23; `norder': the order of the functional curves default to be 6, `bSpline': logical, if True, b-spline is used, otherwise use Fourier basis, default to be True; `Pen': default to be c(0,0), means that the penalty is on the second order derivative of the curve, since the weight for zero-th and first order derivatives of the curve are zero, `lambda':default to be 1e-4, the smoothing parameter for the penalty.

fxList is a list of lists which are similar to fyList. Because it may contain different information for more than one functional covariates.

fbetaList, fbetaList_l and fbetaList_f are similar to each other. They are also expected to be list of lists. The items in each sub-list are: `rtime': range of time, default to be 0 and 1; `nbasis': number of basis functions used in smoothing, default to be less or equal to 19; norder: the order of the functional curves default to be 6;`bSpline': logical, if True, b-spline is used, otherwise use Fourier basis, default to be True; `Pen': default to be c(0,0);`lambda':default to be 1e4;`bivar':logical, if True, the bivariate basis will be calculated, otherwise normal basis, default to be False; `lambdas':the smoothing parameter for the penalty of the additional basis, default to be 1e4.

Note that user only write the item they need to change in the list, all items have default settings. See example below.


}
\value{
A list of
\item{hyper}{Estimated hyper-parameters}
\item{I}{A vector of estimated standard deviation of hyper-parameters}
\item{modellist}{List of models fitted before Gaussian process}
\item{CovFun}{Covariance function}
\item{gamma}{gamma used in Gaussian process powered exponential kernel}
\item{init_resp}{Initial response value}
\item{resid_resp}{Residual after the fitted value from models has been taken out}
\item{fitted}{Fitted value}
\item{fitted.sd}{Standard deviation of the fitted value}
\item{ModelType}{The model applied in the function.}
\item{lTrain}{Training data for functional regression with scalper covariates}
\item{fTrain}{Training data for functional regression with functional covariates}
\item{mfTrainfd}{List of fd objects that from training data for functional regression with functional covariates}
\item{gpTrain}{Training data for Gaussian Process}
\item{time}{Time used in training in Gaussian Process}
\item{iuuL}{Inverse of covariance matrix for lReg}
\item{iuuF}{Inverse of covariance matrix for fReg}
\item{fittedFM}{Fitted value from functional regression}
\item{fyList}{fyList used in the function}
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

traindata=vector('list',20)
for(i in 1:20) traindata[[i]]=i
n=50
traindata=lapply(traindata,function(i) {
  x=seq(-3,3,len=n)
  y=sin(x^2)-x+0.2*rnorm(n,runif(1,-3,3),runif(1,0.5,3))
  x1=0.5*x^3+exp(x)+rnorm(n,runif(1,-3,3),runif(1,0.5,5))
  x2=cos(x^3)+0.2*rnorm(n,runif(1,-3,3),runif(1,0.5,5))
  mat=cbind(x,x1,x2,y)
  colnames(mat)=c('time','x1','x2','y')
  scale=t(c(2*(mean(y)>0.25)-1,(var(y)>3.6)*2-1,(sd(y)-sd(x)>1.4)*2-1))
  i=list(mat,scale)
})

lx=do.call('rbind',lapply(traindata,function(i)i[[2]]))
fx1=do.call('rbind',lapply(traindata,function(i)i[[1]][,2]))
fx2=do.call('rbind',lapply(traindata,function(i)i[[1]][,3]))
fy1=do.call('rbind',lapply(traindata,function(i)i[[1]][,4]))
time_old=traindata[[1]][[1]][,1]

## comment out because running time is a bit long
# system.time(a1<-gpfr(response=(fy1),lReg=lx,fReg=NULL,gpReg=list(fx1,fx2)
#                ,fyList=list(nbasis=23,lambda=0.1),fbetaList_l=
#                list(list(lambda=.01,nbasi=17)),hyper=NULL,
#                Cov=c('pow.ex','linear'),fitting=TRUE,
#                time=seq(-3,3,len=50),rPreIdx=TRUE,concurrent=TRUE))

}

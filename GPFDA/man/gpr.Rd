\name{gpr}
\alias{gpr}

\title{
Gaussian Process regression for single curve
}
\description{
Gaussian Process regression for single curve with train data.
}
\usage{
gpr(Data, response, Cov=c('linear','pow.ex'), hyper=NULL, NewHyper=NULL,
      mean=0, gamma=1,itermax=100,reltol=8e-10,trace=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
The input data from train data. Matrix or vectors are both acceptable. Some data.frames are not acceptable.
}
  \item{response}{
The response data from train data. Matrix or vectors are both acceptable. Some data.frames are not acceptable.
}
  \item{Cov}{
The kernel functions or covariance functions to use. Default is the sum of 'linear' and 'power exponentiation'.
}
  \item{hyper}{
The hyper parameters. Default is NULL. If not NULL, then must be a list with certain names. 
}
  \item{NewHyper}{
Vector of the names of the new hyper parameters from customized kernel function. The names of the hyper-parameters must have the format: xxxxxx.x, i.e. '6 digit' plus 'a dot' plus '1 digit'. This is required for both 'hyper' and 'NewHyper'
}
  \item{mean}{
Is the mean taken out when analysis? Default to be 0, which assumes the mean is zero. if assume mean is a constant, mean=1; if assume mean is a linear trend, mean='t'.
}
  \item{gamma}{
Power parameter used in powered exponential kernel function.
}
  \item{itermax}{
Number of maximum iteration in optimization function. Default to be 100. Normally the number of optimization steps is around 20. If reduce 'reltol', the iterations needed will be less.
}
  \item{reltol}{
Relative convergence tolerance. Smaller reltol means more accurate and slower to converge.
}
  \item{trace}{
The value of the objective function and the parameters is printed every trace'th iteration. Defaults to 0 which indicates no trace information is to be printed.
}
}
\details{
The most important function in the package, for fitting the GP model and store everything necessary for prediction. The optimization used in the function is 'nlminb'. Optimization might break down if the noise for the curve are too far away from normal. Jitter, LU decomposition and sparse matrix inverse are used to ensure the matrix inverse can always get an answer.

The names for the hyper parameters should be:"linear.a" for linear covariance function, "pow.ex.w", "pow.ex.v" for power exponential, "rat.qu.s", "rat.qu.a" for rational quadratic, "vv" for white noise. All hyper parameters should be in one list. 
}
\value{
\item{CovFun}{Covariance function type}
\item{fitted.mean}{Fitted value of training data}
\item{fitted.sd}{Standard deviation of the fitted value of training data}
\item{gamma}{Parameter used in powered exponential covariance function}
\item{hyper}{Hyper-parameter estimated from training data}
\item{I}{Variance of the estimated hyper-parameters}
\item{train.x}{Training covariates}
\item{train.y}{Training response, may be transformed, for prediction use only}
\item{train.yOri}{Original training response}
\item{Q}{Covariance matrix}
\item{inv}{Inverse of the covariance matrix}
\item{mean}{The mean assumed in the analysis}
\item{lrm}{'lm' object if mean is a linear regression. NULL otherwise.}
\item{conv}{0 means converge; 1 otherwise.}
\item{hyper0}{starting point of the hyper-parameters}
}

\references{
  Shi, J Q., and Choi, T. (2011), \emph{Gaussian Process Regression Analysis for Functional Data}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}

\seealso{
  \code{\link{gppredict}};  \code{\link{cov.linear}};  \code{\link{cov.pow.ex}};  \code{\link{cov.rat.qu}};  \code{\link{gpfr}}
}
\examples{
library(GPFDA)
library(MASS) ## used to generate data
hp <- list('pow.ex.w'=log(10),'linear.a'=log(10),'pow.ex.v'=log(5),
      'vv'=log(1))
c <- seq(0,1,len=40)
idx <- sort(sample(1:40,21))
X <- as.matrix(c[idx])
Y <- (mvrnorm(n=40,mu=c-c,Sigma=(cov.linear(hp,c)+cov.pow.ex(hp,c)))[,1]
      )*0.1+sin(c*6)
Y <- as.matrix(Y[idx])
x <- as.matrix(seq(0,1,by=0.03))
a <- gpr(X,Y,c('linear','pow.ex'))

## NOT RUN
## Further codes to provide preditions and plot can be found in demos, for example
## > demo('gpr_ex1')
## END
}

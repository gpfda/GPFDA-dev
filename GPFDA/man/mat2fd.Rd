\name{mat2fd}
\alias{mat2fd}

\title{
Create an fd object from a matrix
}
\description{
Easy setting up for creating an fd object
}
\usage{
mat2fd(mat,fdList=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
Input data, should be a matrix with ncol time points and nrow replications or samples.
}
  \item{fdList}{
A list with following items: `time': a sequence of time points default to be 100 points from 0 to 1; `nbasis': number of basis functions used in smoothing, default to be less or equal to 23; `norder': the order of the functional curves default to be 6, `bSpline': logical, if True, b-spline is used, otherwise use Fourier basis, default to be True; `Pen': default to be c(0,0), means that the penalty is on the second order derivative of the curve, since the weight for zero-th and first order derivatives of the curve are zero, `lambda':default to be 1e-4, the smoothing parameter for the penalty.
}
}
\details{
All items listed above have default values. If any item is required to change, add that item into the list, otherwise leave it as NULL. For example,if one only wants to change the number of basis functions, do: mat2fd{SomeMatrix,list(nbasis=21)}
}

\value{
\item{matfd}{An fd object}
}
\references{
  Shi, J Q., and Choi, T. (2011), \emph{Gaussian Process Regression Analysis for Functional Data}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}
\seealso{
  \code{\link{cov.linear}},\code{\link{xixj_sta}}
}
\examples{
ry=rnorm(20,sd=10)
y1=matrix(0,ncol=100,nrow=20)
for(i in 1:20)  y1[i,]=sin(seq(-1,pi,len=100))*ry[i]

y1fd=mat2fd(y1)
y1fd=mat2fd(y1,list(lambda=1))
}


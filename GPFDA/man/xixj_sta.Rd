\name{xixj_sta}
\alias{xixj_sta}

\title{
Stationary kernel function component.
}
\description{
Component of the distance to build a stationary kernel function or similar. 

\deqn{M=\sum {w_i*(x'_i-x_i^{T})}^{power}}

where \eqn{x_i} is the \eqn{i^{th}} column of the input matrix; \eqn{w_i} is the \eqn{i^{th}} element of the weight vector. Note that \eqn{x} and \eqn{x'} might be different.
}
\usage{
xixj_sta(mat,mat.new=NULL,w=NULL,power=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
Input data, could be a matrix or a vector.
}
  \item{mat.new}{
Second input data, could be a vector or a matrix. Default to be NULL. If NULL, mat.new=mat.
}
  \item{w}{
Weight to be add on each column of the matrix.
}
  \item{power}{
Argument 'power' X 2 will be the power to put on the distance. Default power is 1, which means \eqn{distance^2}. The range of the power to put on the distance is 0 to 2, thus argument 'power' is from 0 to 1.
}
}
\details{
If one wants to involve stationary kernel components in customized covariance matrix, this function will be used in derivatives of the kernel function. See examples in demo('co2').
}
\value{
\item{out}{A symmetric matrix used to build the linear kernel or similar}
}
\references{
  Shi, J Q., and Choi, T. (2011), \emph{Gaussian Process Regression Analysis for Functional Data}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}

\seealso{
  \code{\link{cov.linear}},
}



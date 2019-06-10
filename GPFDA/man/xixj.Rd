\name{xixj}
\alias{xixj}

\title{
Linear kernel function component.
}
\description{
Component to build a linear kernel function or similar. 

\deqn{M=\sum a_i*x'_i*x_i^{T}}

where \eqn{x_i} is the \eqn{i^{th}} column of the input matrix; \eqn{a_i} is the i'th element of the weight vector. Note that \eqn{x} and \eqn{x'} might be different. It is for non-stationary kernel functions.
}
\usage{
xixj(mat,mat.new=NULL,a=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
Input data, could be a matrix or a vector.
}
  \item{mat.new}{
Second input data, could be a vector or a matrix. Default to be NULL. If NULL, mat.new=mat.
}
  \item{a}{
Weight to be add on each column of the matrix.
}
}
\details{
When all 'a' are 1, this is simply mat\%*\%t(mat.new). If one wants to involve linear kernel components in customized covariance matrix, this function will be used in derivatives of the kernel function. See examples in demo('co2').
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
  \code{\link{cov.linear}},\code{\link{xixj_sta}}
}



\name{D2}
\alias{D2}

\title{
Second derivative of the likelihood
}
\description{
Computer the second derivative of the likelihood function with respect to one of the hyper-parameters, with first and second derivative of the kernel function given.
}
\usage{
D2(d1, d2, inv.Q, Alpha.Q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d1}{
First derivative of the kernel function with respect to the required hyper-parameter.
}
  \item{d2}{
Second derivative of the kernel function with respect to the required hyper-parameter.
}
  \item{inv.Q}{
Inverse matrix of the covariance matrix
}
  \item{Alpha.Q}{
  This is \emph{iQY(iQY)'-iQ}. where iQ is the inverse of the covariance matrix, Y is the response.
}
}
\details{
The function is to calculate the second derivative of the normal likelihood, using the first and second derivative of the kernel functions. The first and second derivative need to be pre-defined, for example of customized covariance function, see "demo('co2')".
}
\value{
\item{out}{A number}
}
\references{
  Shi, J Q., and Choi, T. (2011), \emph{Gaussian Process Regression Analysis for Functional Data}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}
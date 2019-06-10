\name{plot.gpfr}
\alias{plot.gpfr}

\title{
Plot Gaussian Process regression with functional mean for either training or predicting
}
\description{
Plot Gaussian Process with functional mean for training or predicting with 'gpfr' class object.
}
\usage{
\method{plot}{gpfr}(x,...,type=c('raw','fitted','prediction'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The 'gpr' object from either training or predicting of the Gaussian Process.
}
  \item{...}{
Other arguments from general 'plot' function, such as: 'axes', etc.
}
  \item{type}{
Function provides three types of plots: raw, fitted and prediction.
}
}

\author{
Jian Qing Shi & Yafeng Cheng
}

\seealso{
  \code{\link{gpfr}};  \code{\link{gpfrpred}};  \code{\link{plot}};  \code{\link{plot.gpr}}
}

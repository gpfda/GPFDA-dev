\name{plot.gpr}
\alias{plot.gpr}

\title{
Plot Gaussian Process training or predicting
}
\description{
Plot Gaussian Process training or predicting only for 'gpr' class object.
}
\usage{
\method{plot}{gpr}(x,...,fitted=F,col.no=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The 'gpr' object from either training or predicting of the Gaussian Process.
}
  \item{...}{
Other arguments from general 'plot' function, such as: 'axes', etc.
}
  \item{fitted}{
Plot fitted value or not. Default to be FALSE, which is to plot the predictions. 
}
  \item{col.no}{
Column number of the input matrix. If the input matrix has more than one columns, than one of them will be used in the plot. Default to be the first one.
}
}

\author{
Jian Qing Shi & Yafeng Cheng
}

\seealso{
  \code{\link{gppredict}};  \code{\link{gpr}};  \code{\link{plot}}
}

\name{betaPar}
\alias{betaPar}

\title{
Create an fdPar object
}
\description{
Easy setting up for create a fdPar object.
}
\usage{
betaPar(betaList=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betaList}{
A list contain following items: `rtime': range of time, default to be 0 and 1; `nbasis': number of basis functions used in smoothing, default to be less or equal to 19; norder: the order of the functional curves default to be 6;`bSpline': logical, if True, b-spline is used, otherwise use Fourier basis, default to be True; `Pen': default to be c(0,0);`lambda':default to be 1e4;`bivar':logical, if True, the bivariate basis will be calculated, otherwise normal basis, default to be False; `lambdas':the smoothing parameter for the penalty of the additional basis, default to be 1e4.
}
}
\details{
All items listed above have default values. If any item is required to change, add that item into the list, otherwise leave it as NULL. For example,if one only wants to change the number of basis functions, do: betaPar{list(nbasis=11)}
}

\value{
\item{betaPar}{An fdPar object}
}
\references{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
}
\author{
Jian Qing Shi & Yafeng Cheng
}
\seealso{
  \code{\link{cov.linear}},\code{\link{xixj_sta}}
}

\examples{
library(GPFDA)
beta1=betaPar()
beta2=betaPar(list(nbasis=7,lambda=0.01))
}

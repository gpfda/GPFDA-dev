\name{GPFDA-package}
\alias{GPFDA-package}
\alias{GPFDA}
\docType{package}
\title{
Gaussian Process in Functional Data Analysis
}
\description{
uses functional regression to be the mean function, and the Gaussian Process to be the covariance structure. 

\deqn{y_m(t)=mu_m(t)+tau_m(x)+epsilon_m(t)}

Where \eqn{m} is the \eqn{m^{th}} data or curve; \eqn{\mu_m} is from functional regression; and \eqn{\tau_m} is from Gaussian Process regression with mean 0 covariance matrix \eqn{k({\bf\theta})}.
}
\details{
\tabular{ll}{
Package: \tab GPFDA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-09-30\cr
License: \tab \cr
}
}
\author{
  Jian Qing Shi & Yafeng Cheng

  Maintainer: yafeng.cheng@ncl.ac.uk
}
\references{
  Shi, J Q., and Choi, T. (2011), \emph{Gaussian Process Regression Analysis for Functional Data}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
}
\keyword{functional, Gaussian Process}


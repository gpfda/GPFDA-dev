#'GPFDA: A package for Gaussian Process Regression for Functional Data Analysis
#'
#'Gaussian Process Regression for Functional Data Analysis
#'
#'@details The main functions of the package are:
#'
#'  \describe{ \item{gpr}{Gaussian process regression using stationary separable
#'  covariance kernels.} \item{nsgpr}{Gaussian process regression using
#'  nonstationary and/or nonseparable covariance kernels.} 
#'  \item{mgpr}{Multivariate Gaussian process -- regression for multivariate outputs.}
#'  \item{gpfr}{Functional regression model given by
#'  \deqn{y_m(t)=\mu_m(t)+\tau_m(x)+\epsilon_m(t),} where \eqn{m} is the
#'  \eqn{m}-th curve or surface; \eqn{\mu_m} is from functional regression;
#'  and \eqn{\tau_m} is from Gaussian Process regression with mean 0 covariance
#'  matrix \eqn{k(\bf \theta)}.}
#'  }
#'@author Jian Qing Shi, Yafeng Cheng, Evandro Konzen
#'@references Shi, J. Q., and Choi, T. (2011), ``Gaussian Process Regression
#'  Analysis for Functional Data'', CRC Press.
#'
#'@docType package
#'@name GPFDA
#'@keywords internal 
"_PACKAGE"

#'@useDynLib GPFDA, .registration=TRUE
#'@importFrom Rcpp sourceCpp
#'@import stats
NULL
#> NULL
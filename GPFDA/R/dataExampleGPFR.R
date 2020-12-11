#' Data simulated in the GPFR example
#'
#' A list containing training and test data simulated from a functional 
#' regression model. \cr \cr
#' In the training set, there are M=20 independent realisations and 
#' the functional response and the functional covariate are observed on a grid 
#' of n=20 time points.  \cr \cr
#' The test set includes a single realisation observed on a grid of n_new=100 
#' time points.  \cr \cr
#' Both training and test sets also have a scalar covariate.
#'
#' @format   A list with seven elements:
#' \describe{ \item{tt}{A vector of length 50} 
#' \item{response_train}{A (20 x 50) matrix}
#' \item{x_train}{A (20 x 50) matrix}
#' \item{scalar_train}{A (20 x 2) matrix}
#' \item{t_new}{A vector of length 100} 
#' \item{response_new}{A vector of length 100}
#' \item{x_new}{A vector of length 100}
#' \item{scalar_new}{A (1 x 2) matrix}
#' } 
#' @details Data used in the GPFR example, see vignette("gpfr").
"dataExampleGPFR"
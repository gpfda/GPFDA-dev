#' Data simulated in the MGPR example
#'
#' A list containing data simulated from a MGPR model. \cr \cr 
#' The dataset contains 30 realisations from a trivariate process. Each of the 
#' three functions is observed on 250 time points on [0,1].
#'
#' @format   A list with two elements:
#' \describe{ \item{input}{List of 3 numeric vectors, each one being the time 
#' points where the corresponding function is observed.} 
#' \item{response}{List of 3 matrices containing the observed 250 datapoints. 
#' Each column is an independent realisation.}
#' } 
#' @details Data used in the MGPR example, see vignette("mgpr").
"dataExampleMGPR"
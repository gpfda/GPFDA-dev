% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gp.functions6.R
\name{plotImage}
\alias{plotImage}
\title{Draw an image plot for a given two-dimensional input}
\usage{
plotImage(
  response,
  input,
  realisation = 1,
  n1,
  n2,
  main = " ",
  zlim = NULL,
  cex.axis = 1,
  cex.lab = 2.5,
  legend.cex.axis = 1,
  font.main = 2,
  cex.main = 2,
  legend.width = 2,
  mar = c(2.1, 2.1, 3.1, 6.1),
  oma = c(0, 1, 0, 0),
  nGrid = 200,
  enlarge_zlim = NULL
)
}
\arguments{
\item{response}{Data to be plotted (e.g. matrix of predictions)}

\item{input}{Matrix of two columns representing the input coordinates.}

\item{realisation}{Integer identifying which realisation should be plotted
(if there are multiple).}

\item{n1}{Number of datapoints in the first coordinate direction}

\item{n2}{Number of datapoints in the second coordinate direction}

\item{main}{Title for the plot}

\item{zlim}{Range of z-axis}

\item{cex.axis}{Graphical parameter}

\item{cex.lab}{Graphical parameter}

\item{legend.cex.axis}{Graphical parameter}

\item{font.main}{Graphical parameter}

\item{cex.main}{Graphical parameter}

\item{legend.width}{Graphical parameter}

\item{mar}{Graphical parameter}

\item{oma}{Graphical parameter}

\item{nGrid}{Dimension of output grid in each coordinate direction}

\item{enlarge_zlim}{Additional quantity to increase the range of zlim}
}
\value{
A plot
}
\description{
Draw an image plot for a given two-dimensional input
}
\examples{
## See examples in vignette:
# vignette("gpr_ex2", package = "GPFDA")
}

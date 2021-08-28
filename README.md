# GPFDA-dev
Developing version of the R package GPFDA

### Short history of the package
The package was initially created to implement the model proposed by Shi, J. Q., et al. "Mixed‐effects Gaussian process functional regression models with application to dose–response curve prediction." Statistics in Medicine 31.26 (2012): 3165-3177. (https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.4502)

The early stage of the functions was in Matlab. Later Yafeng Cheng developed an R version and published on CRAN. A corresponding user manual with theoretical details and examples is available at https://www.staff.ncl.ac.uk/j.q.shi/ps/gpfda.pdf. This manual corresponds to a previous version of the package and will be updated soon. 

In 06/2019, Evandro Konzen joins the project. He brings non-stationary kernels and efficient computation to the package. 

### Changes from 06/2019
* Improvement on computation efficiency using C++ functions.
* Addition of code for GP with nonstationary and/or nonseparable covariance structure -- see Konzen, E., Shi, J. Q. and Wang, Z. (2020) "Modeling Function-Valued Processes with Nonseparable and/or Nonstationary Covariance Structure" (https://arxiv.org/abs/1903.09981).
* Addition of code for multivariate (convolved) GP -- see Chapter 8 of Shi, J. Q., and Choi, T. (2011). "Gaussian process regression analysis for functional data". CRC Press. 
* Addition of vignettes to demonstrate the main functionalities of the package in many GPR and GPFR models.

An animation of a trivariate GP can be seen at https://ekonzen.shinyapps.io/CGP_N3/

#### Changes and additions made to version 4.0.0 (under development) of the package include:
* R methods (print, summary, predict) have been added to GPR, GPFR and MGPR models;
* In GPFR models, documentation has been improved, and bugs and error messages have been fixed;
* Most plots, which were base R plots, have been converted to ggplot2.
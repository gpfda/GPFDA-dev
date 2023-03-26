

### fitFuncReg is for the regression with functional response
fitFuncReg <- function(response,uReg=NULL,fxReg=NULL,fyList=NULL,uCoefList=NULL,
                  fxList=NULL,concurrent=TRUE,fxCoefList=NULL,time=NULL){
  
  y <- response
  uModel <- NULL;res <- NULL;fittedFR <- NULL;
  if(!is.null(uReg)){
    if(inherits(y, "fdata")){
      y <- y$data
    }
  }
  
  if(is.null(time))stop("Argument 'time' must be provided.")
  
  if(!is.null(time)){
    
    fyListSupplied <- ifelse(is.null(fyList), F, T)
    uCoefListSupplied <- ifelse(is.null(uCoefList), F, T)
    fxCoefListSupplied <- ifelse(is.null(fxCoefList), F, T)
    fxListSupplied <- ifelse(is.null(fxList), F, T)
    
    
    if(!fyListSupplied){
      cat('  fyList was not supplied. Default smoothness specifications will be applied to the functional response.','\n')
    }else{
      if(!is.list(fyList)){
        stop('  fyList must be either NULL or a list.','\n')
      }
    }
    
    if(uCoefListSupplied){
      if( !( is.list(uCoefList) & 
             sum(sapply(uCoefList, is.list))==length(uCoefList) ) ){
        stop('  uCoefList must be either NULL or a list of list(s).','\n')
      }
    }
    
    if(fxCoefListSupplied){
      if( !( is.list(fxCoefList) & 
             sum(sapply(fxCoefList, is.list))==length(fxCoefList) ) ){
        stop('  fxCoefList must be either NULL or a list of list(s).','\n')
      }
    }
    
    if(fxListSupplied){
      if( !( is.list(fxList) & 
             sum(sapply(fxList, is.list))==length(fxList) ) ){
        stop('  fxList must be either NULL or a list of list(s).','\n')
      }
    }
    
    if(is.null(fyList)){
      fyList <- list(time=time)
    }else{
      fyList$time <- time
    }
    
    if(is.null(uCoefList)){
      uCoefList <- list(list(rtime=range(time)))
    }else{
      uCoefList <- lapply(uCoefList,function(i) c(i,list(rtime=range(time))))
    }
    
    if(is.null(fxCoefList)){
      fxCoefList <- list(list(rtime=range(time)))
    }else{
      fxCoefList <- lapply(fxCoefList,function(i) c(i,list(rtime=range(time))))
    }
    
    if(is.null(fxList)){
      fxList <- list(list(time=time))
    }else{
      # fxList <- lapply(fxList,function(i) c(i,list(rtime=range(time))))
      fxList <- lapply(fxList,function(i) c(i,list(time=time)))
    }
    
    
    
  }
  if(inherits(y, "matrix")){
    ## define 'fd' object for y if y is a matrix
    y <- mat2fd(y,fyList)
    fyBasis <- y$basis
  }
  if(!inherits(y, "fd")){
    stop('class of response must be one of matrix, fd or fdata')
  }
  y_time <- seq(y$basis$rangeval[1],y$basis$rangeval[2],
                len=length(y$fdnames$time))
  
  if(is.null(uCoefList[[1]]$nbasis)){
    uCoefList <- lapply(uCoefList,function(i)
      c(i,list(nbasis=y$basis$nbasis)))
  }
  if(is.null(uCoefList[[1]]$norder)){
    uCoefList <- lapply(uCoefList,function(i)
      c(i,list(norder=c(fyList$norder,6)[1])))
  }
  if(is.null(uCoefList[[1]]$Pen)) {
    uCoefList <- lapply(uCoefList,function(i){
      if(!is.null(fyList$Pen)) c(i,list(Pen=fyList$Pen))
      if(is.null(fyList$Pen)) c(i,Pen=c(0,0))
    })}
  

  
  
  ## functional response with scalar covariates ------------------------
  
  if(!is.null(uReg)){
    ## define list of x
    if(!inherits(uReg, "matrix")){stop("'uReg' must be a matrix.")}
    x <- uReg
    nx <- ncol(x)
    lxList <- vector('list',length=nx)
    for(i in 1:nx) lxList[[i]] <- x[,i]
    
    ## define list of beta
    if(length(uCoefList)!=length(lxList)){
      if(uCoefListSupplied){
        cat('  The length of uCoefList is not equal to the number of scalar covariates entered in uReg. The first element of uCoefList will be applied to the functional regression coefficient of each scalar covariate.', '\n')
      }else{
        cat('  uCoefList was not supplied. Smoothness specifications in fyList will be applied to each scalar covariate.','\n')
      }
      betalist <- lapply(lxList,function(i){
        i <- betaPar(uCoefList[[1]])
      })
    }
    if(length(uCoefList)==length(lxList))
      betalist <- lapply(uCoefList,betaPar)
    
    
    #regression
    uModel <- fRegress(y, lxList, betalist)
    betaEstMat <- do.call('cbind',lapply(uModel$betaestlist,function(i)
      predict(i,y_time)))
    ml_fitted <- uReg%*%t(betaEstMat)
    
    if(inherits(response, "fd")){
      rawResponse <- eval.fd(y_time,response)
    }
    if(inherits(response, "matrix")){
      if(nrow(response)==nrow(ml_fitted)) residML <- response-ml_fitted
      if(nrow(response)==ncol(ml_fitted)) residML <- t(response)-ml_fitted
    }
    y <- residML
    fittedFRu <- list(ml_fitted)
    names(fittedFRu) <- paste0("fittedFRu")
    res <- c(res,list(y)); fittedFR <- c(fittedFR,fittedFRu)
  }
  
  
  ## functional response with functional covariates ------------------------
  
  mfTrainfd <- NULL
  fxListLambda <- NULL
  fxModel <- list(NULL)
  if(!is.null(fxReg)){

    
    y <- mat2fd(y,fyList)
    fyBasis <- y$basis
    
    ## set up list of 'fd' object for x
    if(inherits(fxReg, "matrix") | inherits(fxReg, "fd")){
      fxReg <- list(fxReg)
    }

    if(inherits(fxReg, "list")){
      if(length(unique(lapply(fxReg, function(x) class(x)[1])))!=1){
        stop('All functional covariates are expected to have the same class.')
      }
      if(unique(lapply(fxReg, function(x) class(x)[1]))=='matrix'){
        if(ncol(fxReg[[1]])!=length(y$fdnames$time)) fxReg <- lapply(fxReg,t)
        if(length(fxList)!=length(fxReg)){
          
          fxListLambda <- rep(NA, length(fxReg))
          for(j in 1:length(fxReg)){
            fxListLambda[j] <- fxList[[1]]$lambda
          }
          
          if(fxListSupplied){
            cat('  The length of fxList is not equal to the number of functional covariates entered in fxReg. The first element of fxList will be applied to each functional covariate.', '\n')
          }else{
            cat('  fxList was not supplied. Default smoothness specifications will be applied to each functional covariate.','\n')
          }
          
          fxReg <- lapply(fxReg,function(i){
            i <- mat2fd(i,fxList[[1]])
          })
          
        }
        if(length(fxList)==length(fxReg)){
          fxListLambda <- rep(NA, length(fxReg))
          for(j in 1:length(fxReg)){
            fxListLambda[j] <- fxList[[j]]$lambda
          }
          
          fxReg <- lapply(1:length(fxReg),function(i){
            mat2fd(fxReg[[i]], fxList[[i]])
          })
        }

      }
      
  
      ## set up list of fdPar object for beta
      if(length(fxCoefList)!=length(fxReg)){
        
        if(fxCoefListSupplied){
          cat('  The length of fxCoefList is not equal to the number of functional covariates entered in fxReg. The first element of fxCoefList will be applied to the functional regression coefficient of each functional covariate.', '\n')
        }else{
          cat('  fxCoefList was not supplied. Default smoothness specifications will be applied to the functional regression coefficient of each functional covariate.','\n')
        }
        
        if(!concurrent) {
          betalist <- lapply(fxReg,function(i) i=betaPar(list(bivar=TRUE)))
        }else{
          betalist <- lapply(fxReg,function(i) i=betaPar(fxCoefList[[1]]))
        }
        
      }
      if(length(fxCoefList)==length(fxReg)){
        if(!concurrent) betalist <- lapply(fxCoefList, function(i)
          betaPar(list(bivar=TRUE)))
        if(concurrent) betalist <- lapply(fxCoefList, betaPar)
      }
      
      ## regression
      fxModel <- NULL
      for(i in seq_along(fxReg)){
        if(is.matrix(y)) y <- mat2fd(y,fyList)
        x <- fxReg[[i]]
        if(concurrent){
          const <- rep(1,dim(x$coef)[2])
          xlist <- list(const=const,x=x)
          
          b1 <- betalist[[i]]
          bList <- list(const=b1,x=b1)
          mf <- fRegress(y,xlist,bList)
          fxModel <- c(fxModel,list(mf))
          
          # evaluate functional coefficients
          betaEstMat <- list(predict(mf$betaestlist$const,y_time)[,1],
                             x=predict(mf$betaestlist$x,y_time)[,1])
          mf_fitted <- apply(t(eval.fd(y_time,x)),1,function(i)
            i=i*betaEstMat[[2]]+betaEstMat[[1]])
          
          if(inherits(y, "fd")){
            rawResponse <- t(eval.fd(y_time,y))
            residMF <- rawResponse-t(mf_fitted)
          }
          if(inherits(y, "matrix")){
            if(nrow(y)==nrow(mf_fitted)) residMF <- rawResponse-mf_fitted
            if(nrow(y)==ncol(mf_fitted)) residMF <- t(rawResponse)-mf_fitted
          }
          y <- residMF
          fittedFRx_i <- list(t(mf_fitted))
          names(fittedFRx_i) <- paste0("fittedFRx", i)
          res <- c(res,list(y)); fittedFR <- c(fittedFR,fittedFRx_i)
        }
        
        if(!concurrent){
          bList <- list(betaPar(), betalist[[i]])
          mf <- linmod(x,y,bList)
          fxModel <- c(fxModel,list(mf))
          
          betaEstMat <- list(b0=eval.fd(y_time,mf$beta0estfd)[,1],
                             b1=eval.bifd(y_time,y_time,mf$beta1estbifd)[,1])
          mf_fitted <- apply(t(eval.fd(y_time,x)),1,function(i)
            i=i%*%betaEstMat[[2]]/length(y_time)^2+betaEstMat[[1]])
          
          if(inherits(y, "fd")){
            rawResponse <- t(eval.fd(y_time,y))
          }
          if(inherits(y, "matrix")){
            if(nrow(y)==nrow(mf_fitted)) residMF <- rawResponse-mf_fitted
            if(nrow(y)==ncol(mf_fitted)) residMF <- t(rawResponse)-mf_fitted
          }
          y <- residMF
          fittedFRx_i <- list(t(mf_fitted))
          names(fittedFRx_i) <- paste0("fittedFRx", i)
          res <- c(res,list(y)); fittedFR <- c(fittedFR,fittedFRx_i)
        }
      }
      mfTrainfd <- fxReg
    }
  }
  
  
  out <- list('model'=list(uModel=uModel,fxModel=fxModel),'resid'=y,
              'fyBasis'=fyBasis, 'fxListLambda'=fxListLambda,
              'mfTrainfd'=mfTrainfd,'fyl'=fyList,'fittedFR'=fittedFR)
  return(out)
}




#' Create an 'fd' object from a matrix
#'
#' Easy setting up for creating an 'fd' object
#'
#' @param mat Input data, should be a matrix with ncol time points and nrow
#'   replications or samples.
#' @param fdList A list with following items: \describe{ \item{time}{Sequence
#'   of time points (default to be 100 points from 0 to 1).}
#'   \item{nbasis}{Number of basis functions used in smoothing, default to be the minimum between 23 and 20\% of the number of time points. } \item{norder}{Order of the functional curves default
#'   to be 6.} \item{bSpline}{Logical, if TRUE (default), b-Spline basis is
#'   used; otherwise, Fourier basis is used.} \item{Pen}{Default to be c(0,0),
#'   meaning that the penalty is on the second order derivative of the curve,
#'   since the weight for zero-th and first order derivatives of the curve are
#'   set to zero.} \item{lambda}{Smoothing parameter for the penalty. Default to
#'   be 1e-4.} }
#'
#' @details All items listed above have default values. If any item is required
#'   to change, add that item into the list; otherwise, leave it as NULL. For
#'   example, if one only wants to change the number of basis functions, do:
#'
#'   \code{mat2fd(SomeMatrix,list(nbasis=21))}
#'
#' @references Ramsay, J., and Silverman, B. W. (2006),
#   ``Functional Data Analysis'', 2nd ed., Springer, New York.
#' @return An 'fd' object
#' @export
#' @import fda
#' @examples
#' require(fda)
#' require(fda.usc)
#' nrep <- 20   # number of replications
#' n <- 100     # number of time points
#' input <- seq(-1, pi, length.out=n) # time points
#' ry <- rnorm(nrep, sd=10)
#' y <- matrix(NA, ncol=n, nrow=nrep)
#' for(i in 1:nrep)  y[i,] <- sin(2*input)*ry[i]
#'
#' plot.fdata(fdata(y,input))
#'
#' yfd <- mat2fd(y, list(lambda=0.01))
#' plot(yfd)
#'
#' yfd <- mat2fd(y, list(lambda=0.00001))
#' plot(yfd)
#'
mat2fd <- function(mat, fdList=NULL){
  fl <- list(time=seq(0,1,len=ncol(mat)),nbasis=min(as.integer(ncol(mat)/5),23),
             norder=6,bSpline=TRUE,Pen=c(0,0),lambda=1e-4)
  nbasis <- c(fdList$nbasis,fl$nbasis)[1]
  norder <- c(fdList$norder,fl$norder)[1]
  lambda <- c(fdList$lambda,fl$norder)[1]
  bSpline <- c(fdList$bSpline,fl$bSpline)[1]
  time <- list(a=fdList$time,b=fl$time)
  time <- time[[which(unlist(lapply(time,is.null))^2==0)[1]]]
  if(!bSpline) fl$Pen <- c(c(0,(2*pi/diff(range(time)))^2,0))
  Pen <- list(a=fdList$Pen,b=fl$Pen)
  Pen <- Pen[[which(unlist(lapply(Pen,is.null))^2==0)[1]]]
  if(bSpline)  basis <- create.bspline.basis(range(time),nbasis,norder)
  if(!bSpline) basis <- create.fourier.basis(range(time),nbasis,
                                              diff(range(time)))
  Par <- vec2Lfd(Pen,range(time))
  matfd <- smooth.basisPar(time,t(mat),basis,Lfdobj=Par,lambda)$fd
  return(matfd)
}

# Create an fdPar object
#
# Easy setting up for creating an fdPar object.
#
# @param betaList
# \describe{
# \item{rtime}{Range of time, default to be 0 and 1}
#  \item{nbasis}{Number of basis functions used in smoothing, default to be
#   less or equal to 19.}
#  \item{bSpline}{Logical, if TRUE (default), b-Spline basis is
#   used; otherwise, Fourier basis is used.}
#  \item{Pen}{Default to be c(0,0),
#   meaning that the penalty is on the second order derivative of the curve,
#   since the weight for zero-th and first order derivatives of the curve are
#   set to zero.}
#  \item{lambda}{Smoothing parameter for the penalty. Default to be 1e4.}
#  \item{bivar}{Logical. Used for non-concurrent models. If TRUE, bivariate basis is used;
#  if FALSE (default), normal basis is used}
#  \item{lambdas}{Smoothing parameter for the penalty of the additional basis.
#  Default to 1.}
# }
#
# @details All items listed above have default values. If any item is required
#   to change, add that item into the list, otherwise leave it as NULL. For
#   example, if one only wants to change the number of basis functions, do:
#   betaPar{list(nbasis=11)}
# @references  Ramsay, J., and Silverman, B. W. (2006),
#   ``Functional Data Analysis'', 2nd ed., Springer, New York.
# @return A list
# @export
# @import fda
# @examples
# require(fda)
# beta1 <- betaPar()
# beta2 <- betaPar(list(nbasis=7,lambda=0.01))
betaPar <- function(betaList=NULL){
  bl <- list(rtime=c(0,1),nbasis=19,norder=4,bSpline=TRUE,Pen=c(0,0),lambda=1e4,
             bivar=FALSE,lambdas=1)
  nbasis <- c(betaList$nbasis,bl$nbasis)[1]
  norder <- c(betaList$norder,bl$norder)[1]
  lambda <- c(betaList$lambda,bl$lambda)[1]
  bSpline <- c(betaList$bSpline,bl$bSpline)[1]
  bivar <- c(betaList$bivar,bl$bivar)[1]
  rtime <- list(a=betaList$rtime,b=bl$rtime)
  rtime <- rtime[[which(unlist(lapply(rtime,is.null))^2==0)[1]]]
  if(!bSpline) fl$Pen <- c(c(0,(2*pi/diff(rtime))^2,0))
  Pen <- list(a=betaList$Pen,b=bl$Pen)
  Pen <- Pen[[which(unlist(lapply(Pen,is.null))^2==0)[1]]]
  if(bSpline)  basis <- create.bspline.basis(rtime,nbasis,norder)
  if(!bSpline) basis <- create.fourier.basis(rtime,nbasis,diff(rtime))
  Par <- vec2Lfd(Pen,rtime)
  if(!bivar){
    betaPar <- fdPar(basis, Par, lambda)
  }else{
    lambdas <- c(betaList$lambdas,bl$lambdas)
    betaPar <- bifdPar(bifd(matrix(0,nbasis,nbasis), basis, basis),
                       Par, Par, lambda, lambdas)
  }
  return(betaPar)
}


### response is expected to be matrices with ncol replications and nrow
### observations
### Data is expected to be a list with matrices
repgp.loglikelihood <- function(hyper.p,response,Data,Cov,gamma=1,nu=1.5,
                                time=NULL,...){

  single <- rep(0,ncol(response))
  for(i in 1:ncol(response)){
    old_input <- as.matrix(do.call('cbind',lapply(Data,function(j) j=j[,i])))
    single[i] <- gp.loglikelihood2(hyper.p,input=old_input,
                                   response=response[,i,drop=F],
                                   Cov=Cov,gamma=gamma,nu=nu,...)
  }
  out <- sum(single)
  return(out)
}



### response is expected to be matrices with ncol replications and nrow
### observations
### Data is expected to be a list with matrices
repgp.Dloglikelihood <- function(hyper.p, response,Data,Cov,gamma=1,nu=1.5,
                                 time=NULL,...){
  
  single <- matrix(0,ncol=length(hyper.p),nrow=ncol(response))
  for(i in 1:nrow(single)){
    old_input <- as.matrix(do.call('cbind',lapply(Data,function(j) j=j[,i])))
    single[i,] <- gp.Dlikelihood2(hyper.p,input=old_input,
                                  response=response[,i,drop=F],
                                  Cov=Cov,gamma=gamma,nu=nu,...)
  }
  out <- apply(single,2,sum)
}



#' Gaussian process functional regression (GPFR) model
#'
#' Use functional regression (FR) model for the mean structure and Gaussian
#' Process (GP) for the covariance structure. \cr \cr Let 'n' be the number of
#' time points 't' of functional objects and 'nrep' the number of independent
#' replications in the sample.
#'
#' @param response Response data. It can be an 'fd' object or a matrix with
#'   'nrep' rows and 'n' columns.
#' @param time Input 't' of functional objects. It is a numeric vector of
#'   length 'n'.
#' @param uReg Scalar covariates for the FR model. It should be a matrix with
#'   'nrep' rows.
#' @param fxReg Functional covariates for the FR model. It can be a matrix with
#'   'nrep' rows  and 'n' columns, an 'fd' object, or a list of matrices or 'fd'
#'   objects.
#' @param gpReg Covariates in the GP model. It should be a matrix, a numeric vector,
#' an 'fd' object, a list of matrices or a list of 'fd' objects.
#' @param fyList A list to control the smoothing of response.
#' @param uCoefList A list to control the smoothing of the regression
#'   coefficient function of the scalar covariates in the FR model.
#' @param fxList A list to control the smoothing of functional covariates in the
#'   FR model.
#' @param fxCoefList A list to control the smoothing of the regression
#'   coefficient function of functional covariates in the functional concurrent
#'   model.
#' @param concurrent Logical. If TRUE (default), concurrent functional
#'   regression will be carried out; otherwise, the full functional regression
#'   will be carried out.
#' @param Cov Covariance function(s) to use. Options are: 'linear', 'pow.ex',
#'   'rat.qu', and 'matern'. Default to 'power.ex'.
#' @param gamma Power parameter used in powered exponential kernel function. It
#'   must be 0<gamma<=2.
#' @param nu Smoothness parameter of the Matern class. It must be a positive
#'   value.
#' @param hyper Vector of initial hyperparameters. Default to NULL.
#' @param NewHyper Vector of names of new hyperparameters from the customized
#'   kernel function.
#' @param useGradient Logical. If TRUE, first derivatives will be used in the
#'   optimization.
#' @param rel.tol Relative tolerance passed to nlminb(). Default to be 1e-10.
#' @param trace.iter Print the processing of iterations of optimization.
#' @param fitting Logical. If TRUE, GPFR fitting is carried out on the supplied data. Default to TRUE.
#'
#' @importFrom  fda.usc is.fdata
#' @importFrom  fda is.fd
#' @details \code{fyList} is a list with the following items:
#'
#'   \itemize{ \item \code{nbasis}: number of basis functions used in
#'   smoothing, default to be the minimum between 23 and 20\% of the number of time points. \item \code{norder}:
#'   order of the functional curves; default to be 6. \item \code{bSpline}:
#'   logical. If TRUE (default), B-splines basis is used; otherwise, Fourier
#'   basis is used. \item \code{Pen}: default to be c(0,0), meaning that the
#'   penalty is only applied to the second order derivative of the curve, with
#'   no penalty for the zero-th and first order derivatives of the curve. \item
#'   \code{lambda}: smoothing parameter for the penalty, default to be 1e-4. }
#'   \code{fxList} is similar to \code{fyList}. However, it is a list of lists
#'   to allow for different specifications for each functional covariate if
#'   there are multiple ones.
#'
#'   \code{uCoefList} and \code{fxCoefList} are similar to
#'   each other. Each one is expected to be a list of lists. If a list of one
#'   element is provided, then the items of this element are applied to each of
#'   the functional coefficients of scalar covariates and of functional covariates,
#'   respectively.
#'
#'   \itemize{ 
#'   \item \code{nbasis}: nnumber of basis functions used in smoothing, default to be
#'   less than or equal to 19. \item  \code{norder}: order of the functional
#'   curves; default to be 6. \item \code{bSpline}: logical. If TRUE (default),
#'   B-splines basis is used; otherwise, Fourier basis is used. \item
#'   \code{Pen}: default to be c(0,0). \item \code{lambda}: smoothing parameter
#'   for the penalty, default to be 1e4. \item \code{bivar}:logical. Used for
#'   non-concurrent models; if TRUE, bivariate basis will be used; if FALSE
#'   (default), normal basis will be used; see details in
#'   \code{\link[fda]{bifdPar}}. \item \code{lambdas}: smoothing parameter for
#'   the penalty of the additional basis, default to be 1. } Note that all items
#'   have default settings.
#'
#' @references \itemize{ \item Ramsay, J., and Silverman, B. W. (2006),
#'   ``Functional Data Analysis'', 2nd ed., Springer, New York. \item Shi, J.
#'   Q., and Choi, T. (2011), ``Gaussian Process Regression Analysis for
#'   Functional Data'', CRC Press. }
#'
#' @return A list containing: \describe{ 
#'   \item{fitted.mean}{Fitted values obtained by GPFR model}
#'   \item{fitted.sd}{Standard deviation of the fitted values obtained by GPFR model}
#'   \item{FRmodelsList}{List of fitted FR models used for the mean function}
#'   \item{residFR}{Residuals obtained in FR model which were used in GPR model.} 
#'   \item{CovFun}{Covariance kernel used in GPR model} \item{gamma}{Parameter 'gamma' used
#'   in Gaussian process with powered exponential kernel} \item{nu}{Parameter
#'   'nu' used in Gaussian process with Matern kernel} 
#'   \item{hyper}{Marginal ML estimates of hyperparameters of GPR model}
#'   \item{rawResponse}{Raw response data}
#'   \item{uTrain}{Training scalar covariates used in FR model}
#'   \item{fxTrain}{Training functional covariates used in FR model}
#'   \item{gpTrain}{Training covariates used in GPR model}
#'   \item{time}{Input time 't'} 
#'   \item{...}{Other objects used for 'print', 'summary', and 'predict' methods} 
  # \item{FRfit}{Estimates of coefficients of functional objects in FR model}
  # \item{iuuL}{Inverse of cross-product of scalar covariates. This is needed for predictions.}
  # \item{iuuF}{Inverse of cross-product of functional covariates. This is needed for predictions.}
  # \item{fittedFR}{Fitted values from the FR model}
  # \item{fyList}{List containing specifications for the smoothing of response.}
  # \item{mfTrainfd}{List of 'fd' objects used in FR model with functional covariates}
#'   }
#' @export
#'
#' @examples
#' ## See examples in vignette:
#' # vignette("gpfr", package = "GPFDA")
gpfr <- function(response, time, uReg=NULL, fxReg=NULL, gpReg=NULL,
                 fyList=NULL, uCoefList=NULL,
                 fxList=NULL, fxCoefList=NULL,
                 concurrent=TRUE,
                 Cov='pow.ex', gamma=2, nu=1.5, hyper=NULL, NewHyper=NULL, 
                 useGradient=T, rel.tol=1e-10, trace.iter=5, fitting=TRUE){
  
  
  if(!is.matrix(response) & !is.fd(response)){
    stop("'response' must be either a matrix or an 'fd' object.")
  }
  
  if(is.numeric(gpReg)) gpReg <- as.matrix(gpReg)
  if(is.matrix(gpReg)) gpReg <- list(gpReg)
  if(is.list(gpReg)) col.no <- length(gpReg)
  # if(is.matrix(gpReg)) col.no <- 1
  if(!is.matrix(gpReg) & !is.list(gpReg)){
    cat('No gpReg found, doing functional regression only')
    col.no <- 1
  }
  if(is.null(hyper)){
    hyper <- list()
    if(any(Cov=='linear')){
      hyper$linear.a <- rnorm(col.no)
      hyper$linear.i <- log(1)
    }
    if(any(Cov=='pow.ex')){
      hyper$pow.ex.v <- runif(1,-1,1)
      hyper$pow.ex.w <- (-abs(rnorm(col.no)))
    }
    if(any(Cov=='matern')){
      hyper$matern.v <- log(1)
      hyper$matern.w <- rep(log(1), col.no)
    }
    if(any(Cov=='rat.qu')){
      hyper$rat.qu.w <- rnorm(col.no)
      hyper$rat.qu.s <- runif(1,0.01,1)
      hyper$rat.qu.a <- runif(1,0.01,1)
    }
    hyper$vv <- sample(x=c(0.2,0.5),1)
    hyper.nam <- names(hyper)
    
    if(!is.null(NewHyper)){
      hyper.nam <- c(hyper.nam,NewHyper)
      nh.length <- length(NewHyper)
      for(i in 1:nh.length){
        hyper <- c(hyper,runif(1,-1,1))
      }
      names(hyper) <- hyper.nam
    }
  }
  
  #################################################
  
  y <- rawResponse <- response

  model <- fitFuncReg(response=response,uReg=uReg,fxReg=fxReg,fyList=fyList,
                 uCoefList=uCoefList,fxList=fxList,concurrent=concurrent,
                 fxCoefList=fxCoefList,time=time)
  
  # Inverse of t(uReg)%*%uReg
  iuuL <- chol2inv(chol(crossprod(uReg)))
  
  fittedFR <- model$fittedFR
  
  ## convert fd/fdata class to matrix
  resid <- model$resid # residuals
  Data <- gpReg
  if(!inherits(resid, "matrix")){
    stop("Residuals from functional regression are not in a 'matrix'.")
  }
  ftime <- model$fyl$time
  if(!is.null(ftime)) time <- ftime
  if(is.null(ftime) & is.null(time)) stop("Input 'time' must be supplied.")
  
  
  Data_classes <- lapply(Data, function(x) class(x)[1])
  
  # if(unique(unlist(lapply(Data,class)))[1]=='fdata'){
  if(unique(unlist(Data_classes))[1]=='fdata'){
    Data <- lapply(Data,function(i) i=t(i$data))
  }
  
  if(inherits(resid, "matrix")){
    resid <- t(resid)
  }
  if(unique(unlist(Data_classes)[1]=='matrix')){
    Data <- lapply(Data,t)
  }
  
  if(inherits(Data,'fd')){
    Data <- eval.fd(time,Data)
  }
  if(unique(unlist(Data_classes)=='fd')){
    Data <- lapply(Data,function(i) (eval.fd(time,i)))
  }
  
  
  ### this is an approximation of iuu for functional regression
  iuuF <- NULL
  if(!is.null(fxReg)){
    if(concurrent){
      iuuF <- lapply(model$model$fxModel,function(i){
        BetaBasis=i$betaestlist$x$fd$basis
        xBasis=i$xfdlist$x$basis
        out=inprod(BetaBasis,xBasis)
        C=i$xfdlist$x$coefs
        out=out%*%C/(ncol(out)*nrow(C))
        out=tcrossprod(out)
        i=out
      })
    }
  }
  
  
  stepsize <- nrow(resid)  #training data size
  tdsize <- as.integer(stepsize/2)  #choose half data for training
  
  
  sample_idx <- 1:tdsize*2-1
  residFR <- resid;
  gpTrain <- Data ## keep original data before reducing the dimension
  resid <- resid[sample_idx,]
  if(ncol(Data[[1]])!=ncol(resid)){
    Data <- lapply(Data, t)
  }
  Data <- lapply(Data,function(i) i=i[sample_idx,])
  
  init0 <- unlist(hyper)
  
  if("pow.ex"%in%Cov & is.null(gamma)){
    stop("Argument 'gamma' must be informed for pow.ex kernel")
  }
  if("matern"%in%Cov & is.null(nu)){
    stop("Argument 'nu' must be informed for matern kernel")
  }
  
  if("matern"%in%Cov & !is.null(nu)){
    if("matern"%in%Cov & !(nu%in%c(3/2, 5/2)) & useGradient){
      useGradient <- F
      cat("Gradient was not used.
      For Matern kernel, the gradient is only available
          if either nu=3/2 or nu=5/2.
      For other values of 'nu', useGradient is
          automatically set to FALSE.")
    }}
  
  cat('\n')
  cat('   Start optimization: ','\n')
  cat('Iteration; Marginal log-lik; hyperparameters: ','\n')
  
  if(!useGradient){repgp.Dloglikelihood <- NULL}
  pp <- nlminb(init0,repgp.loglikelihood,repgp.Dloglikelihood,
               response=resid,Data=Data,Cov=Cov,gamma=gamma,nu=nu,
               control=list(trace=trace.iter,rel.tol=rel.tol))
  
  cat('    Optimization has been finished.','\n','\n')
  hpMMLE <- pp[[1]]
  
  names(hpMMLE) <- names(init0)
  pp.df <- data.frame(hpMMLE=hpMMLE,ppN=substr(names(init0),1,8))
  names(pp.df) <- c('hpMMLE','ppN')
  hpMMLE <- split(pp.df$hpMMLE,pp.df$ppN)
  
  mean <- t(Reduce('+',fittedFR)) ## mean from functional regression model
  
  fitted <- fitted.sd <- NULL
  n <- length(Cov)
  hyper.cg <- hpMMLE
  if(fitting){
    fitted <- fitted.sd <- matrix(NA, ncol=ncol(residFR),
                                  nrow=nrow(residFR))
    for(i in 1:ncol(resid)){
      dr <- as.matrix(do.call('cbind',lapply(gpTrain,
                                             function(j) j=j[,i])))
      yy <- as.matrix(residFR[,i])
      
      CovList <- vector('list',n)
      for(k in 1:n) CovList[k] <- list(paste0('cov.',Cov[k]))
      
      CovL <- lapply(CovList,function(j){
        f <- get(j)
        if(j=='cov.pow.ex')
          return(f(hyper.cg,dr,dr,gamma=gamma))
        if(j=='cov.matern')
          return(f(hyper.cg,dr,dr,nu=nu))
        if(!(j%in%c('cov.pow.ex', 'cov.matern')))
          return(f(hyper.cg,dr,dr))
      }  )
      if(length(CovL)==1)
        Q <- CovL[[1]]
      if(length(CovL)>1)
        Q <- Reduce('+',CovL)
      
      Q <- Q+diag(exp(hyper.cg$vv),dim(Q)[1])
      invQ <- chol2inv(chol(Q))
      QR <- invQ%*%yy
      AlphaQ <- QR%*%t(QR)-invQ
      yfit <- (Q-diag(exp(hyper.cg$vv),dim(Q)[1]))%*%invQ%*%(yy)+mean[,i]
      s2 <- exp(hyper.cg$vv)*rowSums(
        (Q-diag(exp(hyper.cg$vv),dim(Q)[1]))*t(invQ))
      fitted[,i] <- yfit
      fitted.sd[,i] <- sqrt(s2)
      if(i%%5==0) cat(paste0('Fitting ',i,'th curve;','\n'))
    }
    
    cat(paste0('GPFR model fitted to all curves.','\n'))
  }
  
  ##### Save objects for print and summary methods
  
  FRfit <- list()
  
  yBasisType <- model$fyBasis$type
  yNumBasis <- model$fyBasis$nbasis
  yRange <- model$fyBasis$rangeval
  yLambda <- model$fyl$lambda
  
  FRfit$yBasisType <- yBasisType
  FRfit$yNumBasis <- yNumBasis
  FRfit$yRange <- yRange
  FRfit$yLambda <- yLambda
  
  uCoefs <- NULL
  if(!is.null(model$model$uModel[[1]])){
    
    numScalCovariates <- ncol(uReg)
    uCoefs <- vector("list", numScalCovariates)
    uCoeffsBasisType <- rep(NA, numScalCovariates)
    uCoeffsNumBasis <- rep(NA, numScalCovariates)
    uCoeffsLambda <- rep(NA, numScalCovariates)
    uCoeffsRange <- vector("list", numScalCovariates)
    for(j in 1:numScalCovariates){
      coefs_uj <- c(model$model$uModel$betaestlist[[j]]$fd$coefs)
      names(coefs_uj) <- model$model$uModel$betaestlist[[j]]$fd$fdnames[[1]]
      uCoefs[[j]] <- coefs_uj
      uCoeffsBasisType[j] <- model$model$uModel$betalist[[j]]$fd$basis$type
      uCoeffsNumBasis[j] <- model$model$uModel$betalist[[j]]$fd$basis$nbasis
      uCoeffsRange[[j]] <- model$model$uModel$betalist[[j]]$fd$basis$rangeval
      uCoeffsLambda[j] <- model$model$uModel$betaestlist[[j]]$lambda
      
    }
    
    FRfit$uCoeffs <- uCoefs
    FRfit$uCoeffsBasisType <- uCoeffsBasisType
    FRfit$uCoeffsNumBasis <- uCoeffsNumBasis
    FRfit$uCoeffsRange <- uCoeffsRange
    FRfit$uCoeffsLambda <- uCoeffsLambda

  }
  
  if(!is.null(model$model$fxModel[[1]])){
    
    numFuncCovariates <- length(model$model$fxModel)
    fxCoeffs <- vector("list", numFuncCovariates)
    fxCoeffsBasisType <- rep(NA, numFuncCovariates)
    fxCoeffsNumBasis <- rep(NA, numFuncCovariates)
    fxCoefsRange <- vector("list", numFuncCovariates)
    fxCoeffsLambda <- rep(NA, numFuncCovariates)
    fxRegBasisType <- rep(NA, numFuncCovariates)
    fxRegNumBasis <- rep(NA, numFuncCovariates)
    fxRegRange <- vector("list", numFuncCovariates)
    fxRegLambda <- rep(NA, numFuncCovariates)
    for(j in 1:numFuncCovariates){
      coefs_fxj <- c(model$model$fxModel[[j]]$betaestlist$const$fd$coefs)
      names(coefs_fxj) <- model$model$fxModel[[j]]$betaestlist$x$fd$basis$names
      fxCoeffs[[j]] <- coefs_fxj
      fxCoeffsBasisType[j] <- model$model$fxModel[[j]]$betaestlist$x$fd$basis$type
      fxCoeffsNumBasis[j] <- model$model$fxModel[[j]]$betaestlist$x$fd$basis$nbasis
      fxCoefsRange[[j]] <- model$model$fxModel[[j]]$betaestlist$x$fd$basis$rangeval
      
      fxCoeffsLambda[j] <- model$model$fxModel[[j]]$betaestlist$x$lambda
      
      
      fxRegBasisType[j] <- model$mfTrainfd[[j]]$basis$type
      fxRegNumBasis[j] <- model$mfTrainfd[[j]]$basis$nbasis
      fxRegRange[[j]] <- model$mfTrainfd[[j]]$basis$rangeval
      fxRegLambda[j] <- model$fxListLambda[j]
    }
    
    FRfit$fxCoeffs <- fxCoeffs
    FRfit$fxCoeffsBasisType <- fxCoeffsBasisType
    FRfit$fxCoeffsNumBasis <- fxCoeffsNumBasis
    FRfit$fxCoefsRange <- fxCoefsRange
    FRfit$fxCoeffsLambda <- fxCoeffsLambda
    
    FRfit$fxRegBasisType <- fxRegBasisType
    FRfit$fxRegNumBasis <- fxRegNumBasis
    FRfit$fxRegRange <- fxRegRange
    
    FRfit$fxRegLambda <- fxRegLambda
  }
  
#####

  out <- list(
    'fitted.mean'=fitted,
    'fitted.sd'=fitted.sd,
    'FRmodelsList'=model$model,
    'residFR'=residFR,
    'CovFun'=Cov, 'gamma'=gamma,'nu'=nu,
    'hyper'=hpMMLE, 
    'rawResponse'=rawResponse,
    'uTrain'=uReg,
    'fxTrain'=fxReg,
    'gpTrain'=gpTrain,
    'time'=time,         
    'FRfit'=FRfit, 
    'mfTrainfd'=model$mfTrainfd,
    'iuuL'=iuuL,'iuuF'=iuuF,
    'fittedFR'=fittedFR,
    'fyList'=fyList)
  
  
  class(out) <- 'gpfr'
  
  return(out)
}


#' Prediction of GPFR model
#'
#' Make predictions for test input data based on the GPFR model learnt by the
#' 'gpfr' function. Both Type I and Type II predictions can be made.
#'
#' @param object An object of class 'gpfr' obtained by the the 'gpfr' function.
#' @param testInputGP Test input data for the GP prediction. It must be a numeric
#' vector, a matrix (for the case of multiple functional covariates) or an 'fd' object.
#' @param testTime Test time points for prediction. If NULL, default settings
#'   will be applied.
#' @param uReg Scalar covariates data of a new batch for the FR model.
#' @param fxReg Functional covariates data of a new batch for the FR model. It 
#' must be a numeric vector, a matrix (for the case of multiple functional covariates) 
#' or an 'fd' object.
#' @param gpReg Input data for the GP part used for Type I prediction. It must
#'   be a list of three items. The names of the items must be 'response',
#'   'input', and 'time'. The item 'response' is the observed response for a new
#'   batch; 'input' is the observed functional covariates for a new
#'   batch,;'time' is the observed time for the previous two. If NULL (default),
#'   Type II prediction is carried out.
#' @param GPpredict Logical. If TRUE (default), GPFR prediction is carried out;
#'   otherwise only predictions based on the FR model is carried out.
#' @param ... Not used here.
#'
#' @importFrom  fda.usc is.fdata
#'
#' @details If 'gpReg' is provided, then Type I prediction is made. Otherwise,
#'   Type II prediction is made.
#' @return A list containing: \describe{ \item{ypred.mean}{The mean values of
#'   the prediction.} \item{ypred.sd}{The standard deviation of the
#'   predictions.} \item{predictionType}{Prediction type if  GPFR prediction is
#'   carried out.} \item{train}{All items trained by 'gpfr'.} }
#'
#' @references \itemize{ \item Ramsay, J., and Silverman, B. W. (2006),
#'   ``Functional Data Analysis'', 2nd ed., Springer, New York. \item Shi, J.
#'   Q., and Choi, T. (2011), ``Gaussian Process Regression Analysis for
#'   Functional Data'', CRC Press. }
#' @export
#'
#' @examples
#' ## See examples in vignette:
#' # vignette("gpfr", package = "GPFDA")
predict.gpfr <- function(object, testInputGP, testTime=NULL, uReg=NULL, fxReg=NULL,
                        gpReg=NULL, GPpredict=TRUE, ...){
  
  if(!inherits(object, "gpfr")){
    stop("'object' must be of type 'gpfr'")
  }
  
  train <- object
  
  if(is.null(gpReg)){
    predictionType <- 2
  }else{
    predictionType <- 1
  }
  model <- train$FRmodelsList
  if(is.null(model$uModel) & !is.null(uReg)){
    cat('    model with scalar variable is not found, ignoring uReg','\n')
    uReg <- NULL
  }
  if(is.null(model$fxModel[[1]]) & !is.null(fxReg)){
    cat('    model with functional variable is not found, ignoring uReg','\n')
    fxReg <- NULL
  }
  if(!is.null(model$uModel) & is.null(uReg)){
    stop("     'uReg' must be supplied for the model with scalar covariate(s) ','\n'")
    uReg <- NULL
  }
  if(!is.null(model$fxModel[[1]]) & is.null(fxReg)){
    stop("     'fxReg' must be supplied for the model with functional covariate(s) ','\n'")
    fxReg <- NULL
  }
  
  if(!is.null(model$uModel) | !is.null(model$fxModel))
    rtime <- c(model$uModel$yhatfdobj$fd$basis$rangeval,
               model$fxModel[[1]]$yhatfdobj$basis$rangeval,
               model$fxModel[[1]]$yhatfdobj$fd$basis$rangeval)[1:2]
  if(is.null(model$uModel) & is.null(model$fxModel)) rtime <- c(0,1)
  
  if(inherits(testInputGP, "numeric")){
    testInputGP <- as.matrix(testInputGP)
  }
  if(is.matrix(testInputGP)){
    test <- testInputGP
    test <- t(test)
    if(is.null(testTime)) time <- seq(rtime[1],rtime[2],len=col(test))
    if(!is.null(testTime)) time <- testTime
  }
  
  if(inherits(testInputGP, "fd")){
    if(is.null(testTime)) time <- seq(rtime[1],rtime[2],
                                      len=testInputGP$fdnames$time)
    if(!is.null(testTime)) time <- testTime
    test <- eval.fd(time,testInputGP)
  }
  testtime <- time
  if(is.null(uReg)) uRegList <- NULL
  if(!is.null(uReg)) uRegList <- list(f=uReg)
  timeList <- list(f=time)
  if(is.null(fxReg)) fRegList <- NULL
  if(!is.null(fxReg))  fRegList <- list(f=fxReg)
  
  ml_var <- 0
  mf_var <- 0
  
  if(!is.null(gpReg) & !inherits(gpReg,'list')){
    cat("Type I prediction was expecting 'gpReg' to be a list with a response and
        an input. Type II prediction will be made instead.")
    gpReg <- NULL
  }
  type <- 2
  
  if(!is.null(gpReg) & inherits(gpReg,'list')){
    type <- 1
    nl <- names(gpReg)
    if(sum(c('response','input','time')%in%names(gpReg))!=3)
      stop("Check the names in the 'gpReg' list. These names must be 'response',
           'input' and 'time'.")
    if(!1%in%dim(as.matrix(gpReg$response)))
      stop("Check the dimension of the 'response'.")
    gplReg <- uReg
    gpfReg <- gpReg$input
    gpresp <- gpReg$response
    gptime <- gpReg$time
    yhat_ml_list <- vector('list',length=2)
    if(!is.null(uReg)) uRegList <- list(f=uReg,gp=gplReg)
    if(!is.null(time)) timeList <- list(f=time,gp=gptime)
    if(!is.null(fxReg)) fRegList <- list(f=fxReg,gp=gpfReg)

  }
  
  ### predict uModel for fitFuncReg
  if(is.null(uRegList)) yhat_ml <- 0
  if(!is.null(uRegList)){
    for(ii in seq_along(uRegList)){
      uReg <- uRegList[[ii]]
      time <- timeList[[ii]]
      
      if(!is.null(uReg)){
        if(is.vector(uReg)) uReg <- t(as.matrix(uReg))
        if(is.matrix(uReg)){
          if(length(model$uModel$betaestlist)!=ncol(uReg))
            stop('dimension of uReg does not match the model')
          if(length(model$uModel$betaestlist)==ncol(uReg)){
            x <- model$uModel$betaestlist
            for(i in 1:ncol(uReg))
              x[[i]] <- as.matrix(uReg[,i])
          }
        }
        if(inherits(uReg,"list")){
          if(unique(unlist(lapply(uReg, function(x) class(x)[1])))=='fd')
            x <- lapply(uReg,function(i) t(i$coefs))
          if(unique(unlist(lapply(uReg, function(x) class(x)[1])))=='matrix')
            x <- uReg
        }
        betalist <- lapply(model$uModel$betaestlist,function(i){
          i <- predict(i,time)
        })
        if(ii==1){
          ml_var <- do.call('cbind',x)%*%train$iuuL%*%t(do.call('cbind',x))
        }
        for(i in 1:length(x)){
          x[[i]] <- x[[i]]%*%t(betalist[[i]])
        }
        if(ii==1) yhat_ml <- t(Reduce('+',x))
        if(ii==2) gpyhat_ml <- t(Reduce('+',x))
      }
      
    }
  }
  
  ### predict fxModel for fitFuncReg
  if(is.null(fRegList)) yhat_mf <- gpyhat_mf <- matrix(0,ncol=1,nrow=1)
  if(!is.null(fRegList)){
    for(ii in seq_along(fRegList)){
      fxReg <- fRegList[[ii]]
      time <- timeList[[ii]]

      if(!is.null(fxReg)){
        
        #   #   #   #   #   #   
        
        if(is.numeric(fxReg)){
          fxReg <- list(fxReg)
        }
        if(is.list(fxReg)){
          if(is.numeric(fxReg[[1]])){
            fxReg <- lapply(fxReg, function(x) t(x))
          }
        }
        
        # if(is.numeric(fxReg)){fxReg <- list(as.matrix(fxReg))}
        
        #   #   #   #   #   #   
        
        if(is.fdata(fxReg)) fxReg <- list(t(fxReg$data))
        if(is.matrix(fxReg) | is.fd(fxReg)) fxReg <- list((fxReg))
        if(unique(unlist(lapply(fxReg, function(x) class(x)[1])))=='fd'){
          fxReg <- lapply(fxReg,function(i) eval.fd(time,i))
        }
        if(unique(unlist(lapply(fxReg, function(x) class(x)[1])))=='matrix'){
          fxReg <- lapply(fxReg,t)
        }
        if(unique(unlist(lapply(fxReg,ncol)))>1){
          if(length(fxReg)!=1 & ii==1){
            stop('new samples of functional covariates for functional
                 regression are having wrong dimensions')}
          if(length(fxReg)!=1 & ii==2){
            stop('input functional covariates for Gaussian process are having
                 wrong dimensions')}
          if(length(fxReg)==1){
            ftmp <- vector('list',length=ncol(fxReg[[1]]))
            for(i in seq_along(ftmp)) ftmp[[i]] <- as.matrix(fxReg[[1]][,i])
            fxReg <- ftmp
            rm(ftmp)
          }
        }
        
        fxReg <- lapply(fxReg,function(i){
          i=cbind(matrix(1,nrow=nrow(fxReg[[1]])),i)
        })
        ## find beta and multiply with the fx
        if(!'bifd'%in%unlist(lapply(model$fxModel[[1]], function(x) class(x)[1]))){
          fbeta <- lapply(model$fxModel,function(i){
            intcept=predict(i$betaestlist[[1]],time)
            slope=predict(i$betaestlist[[2]],time)
            i=cbind(intcept,slope)
          })
          
          if(ii==1){
            mf_var <- rep(0,length(fbeta))
            for(i in seq_along(fbeta)){
              iuuTime <- seq(train$mfTrainfd[[i]]$basis$rangeval[1],
                             train$mfTrainfd[[i]]$basis$rangeval[2],
                             len=nrow(fxReg[[i]]))
              Basis <- train$mfTrainfd[[i]]$basis
              fx <- as.matrix(fxReg[[i]][,2])
              fx <- t(
                eval.fd(seq(
                  iuuTime[1],iuuTime[2],
                  len=train$FRmodelsList$fxModel[[1]]$betaestlist$x$fd$basis$nbasis),
                  smooth.basis(iuuTime,fx,Basis)$fd))
              mf_var[i] <- fx%*%train$iuuF[[i]]%*%t(fx)
            }
          }
          # for(i in seq_along(fbeta)){
            fxReg[[i]] <- apply(fxReg[[i]],2,function(j){
              j=as.matrix(fbeta[[i]][,1])+as.matrix(fbeta[[i]][,2])*j
            })
          # }
          if(ii==1) yhat_mf <- Reduce('+',fxReg)
          if(ii==2) gpyhat_mf <- Reduce('+',fxReg)
        }
        
        if('bifd'%in%unlist(lapply(model$fxModel[[1]], function(x) class(x)[1]))){
          fbeta <- lapply(model$fxModel,function(i){
            intcept=eval.fd(time,i$beta0estfd)
            slope=eval.bifd(time,time,i$beta1estbifd)
            i=cbind(intcept,slope)
          })
          for(i in seq_along(fbeta)){
            fxReg[[i]] <- apply(fxReg[[i]],2,function(j){
              j <- as.matrix(fbeta[[i]][,1]) +
                as.matrix(fbeta[[i]][,-1])%*%as.matrix(j)
            })}
          if(i==1) yhat_mf <- Reduce('+',fxReg)
          if(i==2) gpyhat_mf <- Reduce('+',fxReg)
        }
        
      }
      
      
    }
  }
  #
  f.mean <- yhat_ml+apply(yhat_mf,1,sum)
  f.var <- apply((train$rawResponse-t(train$residFR))^2,2,sum)/(nrow(
    train$rawResponse)-1)
  train$fyList$time <- train$time
  f.var <- mat2fd(t(as.matrix(f.var)),train$fyList)
  f.var <- abs(eval.fd(testtime,f.var))
  # Fypredup <- f.mean + 1.96*sqrt(f.var)
  # Fypredlo <- f.mean - 1.96*sqrt(f.var)
  if(!GPpredict) {
    # return(list(ypred=cbind(f.mean,Fypredup,Fypredlo),
    #                            unclass(train)))
    result <- c(list(ypred.mean = f.mean,
                     ypred.sd = sqrt(f.var),
                     testTime=testtime),
                unclass(train))
  }else{
    
    ## Type I prediction
    if(type==1){
      
      gpresp <- as.matrix(gpReg$response)
      if(nrow(gpresp)==1) gpresp <- t(gpresp)
      ygpobs <- as.matrix(gpresp)-gpyhat_ml-apply(gpyhat_mf,1,sum)
      y_gppred <- predict.gpr(object=NULL,hyper=train$hyper,
                             input=as.matrix(gpReg$input), Y=as.matrix(ygpobs),
                             inputNew=t(as.matrix(test)), Cov=train$CovFun,
                             gamma=train$gamma, nu=train$nu)
      ygppred <- y_gppred$pred.mean
      s2 <- ((y_gppred$pred.sd)^2-exp(train$hyper$vv))%*%(1 + ml_var)
      #+sum(mf_var))
      ypred <- yhat_ml+apply(yhat_mf,1,sum) + ygppred
      ## fd regression plus gp regression
      cat('  Type I predictions calculated.','\n')
    }
    
    ## Type II prediction
    if(is.null(gpReg)|type==2){
      
      fitted <- matrix(0,ncol=ncol(train$residFR),nrow=ncol(test))
      fitted.var <- fitted
      
      for(i in 1:ncol(train$residFR)){
        input <- do.call('cbind',lapply(train$gpTrain,function(j) j=j[,i]))
        y_gppred <- predict.gpr(object=NULL,inputNew=t(as.matrix(test)),
                               hyper=train$hyper,input=as.matrix(input),
                               Y=as.matrix(train$residFR[,i]),
                               Cov=train$CovFun,gamma=train$gamma, nu=train$nu)
        ygppred <- y_gppred$pred.mean
        s2 <- ((y_gppred$pred.sd)^2-exp(train$hyper$vv))
        #%*%(1 + ml_var+sum(mf_var))
        fitted[,i] <- yhat_ml+apply(yhat_mf,1,sum) + ygppred
        fitted.var[,i] <- s2
      }
      ypred <- as.matrix(apply(fitted,1,mean))
      s2 <- as.matrix(apply(fitted.var,1,mean))+apply(fitted^2,1,mean)-ypred^2
      cat('  Type II predictions calculated.','\n')
    }
    
    result <- c(list(ypred.mean=ypred,
                     ypred.sd=sqrt(s2),
                     testTime=testtime,
                     predictionType=predictionType),
                unclass(train))

  }
  
  class(result) <- 'gpfr'
  return(result)
}









#' Print method for 'gpfr' objects
#'
#' @param x A 'gpfr' object.
#' @param ... Not used here.
#'
#' @return Returns a summary of an object of class 'gpfr'.
#' @export
#'
#' @seealso \link[GPFDA]{gpfr}
print.gpfr <- function(x, ...) {
  
  if(!inherits(x,"gpfr")){
    stop("'x' must be of type gpfr")
  }
  
  object <- x
  
  cat(paste0("=============================================================", "\n"))
  
  cat(paste0("Smoothness specifications for the functional response: \n \n"))
  
  df <- as.data.frame(matrix(NA, 1, 4))
  names(df) <- c("nbasis", "type", "lambda", "range")
  range_j <- object$FRfit$yRange
  df[1, ] <- c(object$FRfit$yNumBasis,
               object$FRfit$yBasisType,
               object$FRfit$yLambda,
               paste0("[", range_j[1], ",", range_j[2], "]"))
  
  print(df, row.names = FALSE)
  
  cat(paste0("-----------------------------", "\n"))
  
  cat("\n")
  cat(paste0("FR model for the mean function:", "\n"))
  
  numScalarCovariates <- length(object$FRfit$uCoeffs)
  numFunctionalCovariates <- length(object$FRfit$fxCoeffs)
  
  cat("\n")
  cat(paste0("Number of scalar covariates:  ",  numScalarCovariates,"\n"))
  cat(paste0("Number of functional covariates:  ",  numFunctionalCovariates,"\n"))
  
  cat("\n")
  
  if(!is.null(object$FRfit$uCoeffs)){
    cat(paste0("-----------------------------", "\n"))
    cat(paste0("Smoothness specifications for the functional regression coefficient of
    each scalar covariate: \n \n"))
    
    df <- as.data.frame(matrix(NA, numScalarCovariates, 4))
    names(df) <- c("nbasis", "type", "lambda", "range")
    for(j in 1:numScalarCovariates){
      range_j <- object$FRfit$uCoeffsRange[[j]]
      df[j, ] <- c(object$FRfit$uCoeffsNumBasis[j],
                   object$FRfit$uCoeffsBasisType[j],
                   object$FRfit$uCoeffsLambda[j],
                   paste0("[", range_j[1], ",", range_j[2], "]"))
    }
    print(df)
  }
  
  
  
  
  if(!is.null(object$FRfit$fxRegBasisType)){
    cat("\n")
    cat(paste0("-----------------------------", "\n"))
    
    cat(paste0("Smoothness specifications for each functional covariate: \n \n"))
    
    df <- as.data.frame(matrix(NA, numFunctionalCovariates, 4))
    names(df) <- c("nbasis", "type", "lambda", "range")
    for(j in 1:numFunctionalCovariates){
      range_j <- object$FRfit$fxRegRange[[j]]
      df[j, ] <- c(object$FRfit$fxRegNumBasis[j],
                   object$FRfit$fxRegBasisType[j],
                   object$FRfit$fxRegLambda[j],
                   paste0("[", range_j[1], ",", range_j[2], "]"))
    }
    print(df)
    
  }
  
  if(!is.null(object$FRfit$fxCoeffs)){
    cat("\n")
    cat(paste0("-----------------------------", "\n"))
    
    cat(paste0("Smoothness specifications for the functional regression coefficient of each functional covariate: \n \n"))
    
    
    df <- as.data.frame(matrix(NA, numFunctionalCovariates, 4))
    names(df) <- c("nbasis", "type", "lambda", "range")
    for(j in 1:numFunctionalCovariates){
      range_j <- object$FRfit$fxCoefsRange[[j]]
      df[j, ] <- c(object$FRfit$fxCoeffsNumBasis[j],
                   object$FRfit$fxCoeffsBasisType[j],
                   object$FRfit$fxCoeffsLambda[j],
                   paste0("[", range_j[1], ",", range_j[2], "]"))
    }
    print(df)
  }
  
  
  cat("\n")
  cat(paste0("=============================================================", "\n"))
  cat(paste0("GPR model for the covariance function:", "\n"))
  
  cat("\n")
  cat(paste0("Covariance kernel: ", object$CovFun, "\n"))
  if(object$CovFun=="matern"){
    cat(paste0("nu = ", object$nu, "\n"))
  }else if(object$CovFun=="pow.ex"){
    cat(paste0("gamma = ", object$gamma, "\n"))
  }
  
  
}


#' Summary method for 'gpfr' objects
#'
#' @param object A 'gpfr' object.
#' @param ... Not used here.
#' 
#' @return Some fitting results of the GPFR model.
#' @export
#'
#' @seealso \link[GPFDA]{gpfr}
summary.gpfr <- function(object, ...){
  
  if(!inherits(object,"gpfr")){
    stop("'object' must be of type gpfr")
  }
  
  cat(paste0("=============================================================", "\n"))
  cat(paste0("-----------------------------", "\n"))
  cat(paste0("FR model for the mean function:", "\n"))
  
  numScalarCovariates <- length(object$FRfit$uCoeffs)
  numFunctionalCovariates <- length(object$FRfit$fxCoeffs)
  
  cat("\n")
  cat(paste0("Number of scalar covariates:  ",  numScalarCovariates,"\n"))
  cat(paste0("Number of functional covariates:  ",  numFunctionalCovariates,"\n"))
  
  
  if(!is.null(object$FRfit$uCoeffs)){
    cat("\n")
    cat(paste0("-----------------------------", "\n"))
    cat(paste0("Estimates of coefficients of the functional coefficient for each
    scalar covariate:", "\n"))
    
    print(object$FRfit$uCoeffs)
    
  }
  
  
  if(!is.null(object$FRfit$fxCoeffs)){
    cat("\n")
    cat(paste0("-----------------------------", "\n"))
    
    
    
  }
  
  
  
  
  cat("\n")
  cat(paste0("=============================================================", "\n"))
  cat(paste0("GPR model for the covariance function:", "\n"))
  
  cat("\n")
  cat(paste0("Covariance kernel: ", object$CovFun, "\n"))
  if(object$CovFun=="matern"){
    cat(paste0("nu = ", object$nu, "\n"))
  }else if(object$CovFun=="pow.ex"){
    cat(paste0("gamma = ", object$gamma, "\n"))
  }
  
  cat("\n")
  cat(paste0("Marginal ML estimates of hyperparameters:", "\n"))
  print(unlist(object$hyper))
  
  cat("\n")
  
  
}


#' Plot GPFR model for either training or prediction
#'
#' @param x Plot GPFR for training or prediction from a given object of 'gpfr'
#'   class.
#' @param type Required type of plots. Options are: 'raw',
#'   'meanFunction', 'fitted' and 'prediction'.
#' @param ylab Title for the y axis.
#' @param xlab Title for the x axis.
#' @param ylim Graphical parameter. If NULL (default), it is chosen automatically.
#' @param realisations Index vector identifying which training realisations
#' should be plotted. If NULL (default), all training realisations are plotted.
#' For predictions, 'realisations' should be '0' if no training realisation
#' is to be plotted.
#' @param alpha Significance level used for 'fitted' or 'prediction'. Default is 0.05.
#' @param colourTrain Colour for training realisations when 'type' is set to
#' 'prediction' and 'realisations' is positive.
#' @param colourNew Colour for predictive mean for the new curve when 'type' is
#' set to 'prediction'.
#' @param lwd Graphical parameter.
#' @param cex.lab Graphical parameter.
#' @param cex.axis Graphical parameter.
#' @param cex.main Graphical parameter.
#' @param main Title of the plot.
#' @param ... Not used here.
#' @import  ggplot2
#' @importFrom reshape2 melt
#' @return A plot.
#' @export
#' @examples
#' ## See examples in vignette:
#' # vignette("gpfr", package = "GPFDA")
plot.gpfr <- function (x, type=c('raw','meanFunction','fitted','prediction'),
                       ylab='y', xlab='t', ylim=NULL, realisations=NULL,
                       alpha=0.05,
                       colourTrain="red", colourNew="blue", lwd=0.5,
                       cex.lab=10, cex.axis=10, cex.main=15, main=NULL, ...){
  
  if(!inherits(x,"gpfr")){
    stop("'object' must be of type gpfr")
  }
  
  obj <- x
  
  if(missing(type)) {
    type <- 'fitted'
    message("Argument 'type' was not supplied. Fitted GPFR model will be plotted.")
  }
  
  if((length(type)!=1L) |
     ((length(type)=1L)&(!type%in%c('raw','meanFunction','fitted','prediction')))){
    stop("Argument 'type' must a character object of length one. It must be either
         'raw', 'meanFunction', 'fitted' or 'prediction'.")
  }
  
  
  z <- stats::qnorm(1-alpha/2)
  
  if(is.null(realisations)){
    idx <- 1:ncol(obj$fitted.mean)
  }else{
    idx <- realisations
  }
  numRealis <- length(idx)
  
  
  if(type=='raw'){
    df <- data.frame(t=obj$time, y=t(obj$rawResponse[idx,,drop=F]))
    meltdf <- reshape2::melt(df, id="t")
    out <- ggplot(meltdf, aes(x=t,y=value,colour=variable,group=variable)) +
      geom_line(lwd=lwd) + theme(legend.position="none")
    if(numRealis<6){
      out <- out + geom_point(shape = 21, fill = "black", size = 0.2, stroke = 1, aes(color = variable))
    }
    if(is.null(main)){ main <- "Raw data"}
  }
  
  #############################################
  if(type=='meanFunction'){
    
    df <- data.frame(t=obj$time, y=t(obj$rawResponse[idx,,drop=F]))
    meltdf <- reshape2::melt(df, id="t")
    out <- ggplot(meltdf, aes(x=t,y=value,colour=variable,group=variable)) +
      geom_line(lwd=lwd, linetype = "dashed") + theme(legend.position="none")
    if(numRealis<6){
      out <- out + geom_point(shape = 21, fill = "black", size = 0.2, stroke = 1,
                              aes(color = variable))
    }
    
    fittedFRtotal <- t(Reduce("+", obj$fittedFR)[idx, ])
    
    df2 <- data.frame(t=obj$time,  y=fittedFRtotal)
    meltdf2 <- reshape2::melt(df2, id="t")
    out <- out + geom_line(data = meltdf2)
    
    if(is.null(main)){ main <- "Mean function fit"}
  }
  
  #############################################
  if(type=='fitted'){
    
    df <- data.frame(t=obj$time, y=obj$fitted.mean[,idx,drop=F])
    meltdf <- reshape2::melt(df, id="t")
    out <- ggplot(meltdf, aes(x=t,y=value,colour=variable,group=variable)) +
      geom_line(lwd=lwd) + theme(legend.position="none")
    
    df2 <- data.frame(t=obj$time,
                      y=t(obj$rawResponse[idx,,drop=F]))
    meltdf2 <- reshape2::melt(df2, id="t")
    out <- out + geom_point(data=meltdf2, shape = 21, fill = "black", 
                            size = 0.2, stroke = 1,
                            aes(color = variable))
    
    yCIlo <- (obj$fitted.mean-obj$fitted.sd*z)[,idx[1:numRealis]]
    yCIup <- (obj$fitted.mean+obj$fitted.sd*z)[,idx[1:numRealis]]
    
    out <- out + geom_ribbon(aes(ymin=c(yCIlo), ymax=c(yCIup)), alpha=0.1, linewidth=0)
    
    if(is.null(main)){ main <- "GPFR fit"}
  }
  
  
  #############################################
  if(type=='prediction'){
    
    if(is.null(ylim)){
      if(identical(idx, 0)){
        ylim <- range(c(obj$ypred.mean - z*obj$ypred.sd),
                      c(obj$ypred.mean + z*obj$ypred.sd))
        # ylim <- range(c(obj$ypred[,2], obj$ypred[,3], obj$ypred[,1]))
      }else{
        ylim <- range(c(obj$rawResponse[idx,]),
                      c(obj$ypred.mean - z*obj$ypred.sd),
                      c(obj$ypred.mean + z*obj$ypred.sd))
      }
    }
    
    
    if(identical(idx, 0)){
      
      df <- data.frame(t=obj$testTime,
                        y=obj$ypred.mean[,1])
      meltdf <- reshape2::melt(df, id="t")
      
      out <- ggplot(data=meltdf, aes(x=t,y=value,colour=colourNew,group=variable),
                    lwd=lwd, linetype="dashed",
                    color=colourNew) +
        geom_line(lwd=lwd, colour=colourNew) + theme(legend.position="none") +
        geom_point(data=meltdf, shape = 21, fill = "black",
                   size = 1, stroke = 1, color = colourNew)
      
      
      yCIlo <- c(obj$ypred.mean - z*obj$ypred.sd)
      yCIup <- c(obj$ypred.mean + z*obj$ypred.sd)
      
      df2 <- data.frame(t=obj$testTime, yCIlo=yCIlo, yCIup=yCIup)
      out <- out + geom_ribbon(data=df2, aes(x=t, ymin=c(yCIlo), ymax=c(yCIup)),
                               alpha=0.3, linewidth=0, inherit.aes=FALSE)
      
    }else{
      
      df <- data.frame(t=obj$time, y=t(obj$rawResponse[idx,,drop=F]))
      meltdf <- reshape2::melt(df, id="t")
      out <- ggplot(meltdf, aes(x=t,y=value,colour=colourTrain,group=variable)) +
        geom_line(lwd=lwd) + theme(legend.position="none")
      
      if(numRealis<6){
        out <- out + geom_point(shape = 21, fill = "black", size = 0.2, stroke = 1,
                                aes(color = colourTrain))
      }
      
      df2 <- data.frame(t=obj$testTime,
                        y=obj$ypred.mean[,1])
      meltdf2 <- reshape2::melt(df2, id="t")
      out <- out + geom_line(data=meltdf2, lwd=lwd, linetype="dashed",
                             color=colourNew) +
        geom_point(data=meltdf2, shape = 21, fill = "black",
                   size = 1, stroke = 1, color = colourNew)
      
      yCIlo <- c(obj$ypred.mean - z*obj$ypred.sd)
      yCIup <- c(obj$ypred.mean + z*obj$ypred.sd)
      
      df3 <- data.frame(t=obj$testTime, yCIlo=yCIlo, yCIup=yCIup)
      out <- out + geom_ribbon(data=df3, aes(x=t, ymin=c(yCIlo), ymax=c(yCIup)),
                               alpha=0.3, linewidth=0, inherit.aes=FALSE)
      
    }
    
    if(is.null(main)){
      if(obj$predictionType==1){
        main <- "Type I prediction"
      }else{
        main <- "Type II prediction"
      }
    }
    
  }
  #############################################
  
  out <- out + labs(title = main) + theme(plot.title = element_text(hjust = 0.5))
  
  if(!is.null(cex.main)){
    out <- out + theme(plot.title = element_text(size=cex.main))
  }
  if(!is.null(xlab)){
    out <- out + labs(x = xlab)
  }
  if(!is.null(xlab)){
    out <- out + labs(y = ylab)
  }
  if(!is.null(ylim)){
    out <- out + coord_cartesian(ylim = ylim)
  }
  
  if(!is.null(cex.lab)){
    out <- out + theme(axis.title=element_text(size=cex.lab))
  }
  if(!is.null(cex.axis)){
    out <- out + theme(axis.text=element_text(size=cex.axis))
  }
  
  suppressWarnings(print(out))
  
}



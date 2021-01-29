

### main1 is the function for the case of scalar response variable. 
### This function will do lm, flm, gp, or combinations 
### of the three, depends on the avaliable data.
main1 <- function(response,uReg=NULL,fxReg=NULL,fxList=NULL,fxCoefListScalarResp=NULL){ 
  
  y <- response
  
  ## multivariate linear regression
  ml <- NULL
  if (!is.null(uReg)){
     #  if(class(uReg)!='matrix') stop('The covariates for scalar multivariate 
     #  linear regression are expected to store in a matrix, not other format')
    ml <- lm(y~uReg)
    resid_ml <- as.matrix(y-resid(ml))
    y <- resid_ml
  }
  
  ## functional regression
  temp <- list(NULL)
  if(!is.null(fxReg)){
    if(class(fxReg)=='matrix'|class(fxReg)=='fd')
      fxReg <- list(fxReg)
    if(class(fxReg)=='list'){
      if(length(unique(unlist(lapply(fxReg,class))))!=1) 
        stop('functional covariates are expected to have same class')
      if(unique(unlist(lapply(fxReg,class)))=='matrix'){
        temp <- list(NULL)
        res <- list(y)
        
        ## set up functional variable for fx
        if(length(fxList)!=length(fxReg)){
          cat('  The length of fxList is not equal to the number of functional covariates entered in fxReg','\n')
          
          if(length(fxList)==0){
            cat('     Default fxList is applied','\n')
            fxList <- lapply(fxReg,function(i){
              i <- min(as.integer(ncol(i)/5),10)
              names(i) <- list('nx_basis')
              return(i)
            })
          } 
          else{
            cat('     The first element of fxList will be applied to all functional covariates.','\n')
            fxList <- lapply(fxReg,function(i){
              i <- fxList[[1]]
              return(i)
            })
          }
          
        }
          
        ## set up functional parameter for fbeta
        if(length(fxCoefListScalarResp)!=length(fxReg)){
          cat('       Length of fxCoefListScalarResp is not equal to the length of list of  functional covariates','\n')
          if(length(fxCoefListScalarResp)==0){
            cat('     Default fxCoefListScalarResp is applied','\n')
            fxCoefListScalarResp <- fxList
            fxCoefListScalarResp <- lapply(fxCoefListScalarResp,function(i){
              names(i) <- 'nbeta_basis'
              return(i)
            })
          }
          else{
            cat('     First element of fxCoefListScalarResp is applied to all items','\n')
            fxCoefListScalarResp <- lapply(fxReg,function(i){
              i <- fxCoefListScalarResp[[1]]
              return(i)
            })
          }
        }
          
        for(i in seq_along(fxReg)){
        ## functional regression with scalar response and functional covariates 
        ## using fdata.
          mf <- freg1(y,fxReg[[i]],fxList[[i]],fxCoefListScalarResp[[i]])
          res <- list(res,as.matrix(resid(mf)))
          temp <- c(temp,list(mf))
          if(is.null(temp[[1]])) temp <- temp[-1]
          y <- res[[length(res)]]
        } 
      }
      if(unique(unlist(lapply(fxReg,class)))=='fd'){
        res <- list(y)
        ## set up functional parameter for fbeta
        if(length(fxCoefListScalarResp)!=length(fxReg)){
          cat('     Length of fxCoefListScalarResp is not equal to the length of list of 
              functional covariates','\n')
          if(length(fxCoefListScalarResp)==0){
            cat('     Default fxCoefListScalarResp is applied','\n')
            fbeta[[i]] <- betaPar()
          }
          else{
            cat('     First element of fxCoefListScalarResp is applied to all items','\n')
            fbeta <- lapply(fxReg,function(i){
              i <- betaPar(fxCoefListScalarResp[[1]])
              return(i)
            })
          }
        }
        if(length(fxCoefListScalarResp)==length(fxReg))
          fbeta <- lapply(fxCoefListScalarResp,betaPar)
          
        for(i in seq_along(fxReg)){          
        ## functional regression with scalar response and functional covariates 
        ## using 'fd'.
          mf <- freg2(y,fxReg[[i]],fbeta[[i]])
          res <- list(res,as.matrix(resid(mf)))
          temp <- c(temp,list(mf))
          if(is.null(temp[[1]])) temp <- temp[-1]
          y <- res[[length(res)]]
        } 
      }
    }
  }
  out <- list(model=list(ml=ml,mf=temp),res=y)
  return(out)
}

#' @importFrom  fda.usc fdata
#' @importFrom  fda.usc fregre.basis
#' @importFrom  fda create.bspline.basis
freg1 <- function(y,fx,nx_bas,nbeta_bas){
  ### works for scalar response and a single functional covariate;
  ### y is expected to be a column matrix; 
  ### fx is expected to be a matrix.
  if (nrow(fx)!=nrow(y)) stop('unequal sample size for response and functional
                              covariates')
  fx <- fdata(fx,argvals=seq(0,1,len=ncol(fx)))
  tt <- fx[['argvals']]
  basis1 <- create.bspline.basis(rangeval=range(tt),nbasis=nx_bas)
  basis2 <- create.bspline.basis(rangeval=range(tt),nbasis=nbeta_bas)
  fm <- fregre.basis(fx,y,basis1,basis2)
  return(fm)
}

freg2 <- function(y,xlist,fbeta){
  ### works for scalar response and a single functional covariate;
  ### y is expected to be a column matrix; 
  ### fx is expected to be an 'fd' object;
  ### fbeta is expected to be an 'fd' onject.
  fm <- fRegress(y[,1], xlist, fbeta)
  return(fm)
}



main2 <- function(response,uReg=NULL,fxReg=NULL,fyList=NULL,uCoefList=NULL,
               fxList=NULL,concurrent=TRUE,fxCoefList=NULL,time=NULL){
  ### main2 is for the regression with functional response, 
  ### This function will do lm, flm, gp, or conbinations 
  ### of the three, depends on the avaliable data.
  y <- response
  ml <- NULL;res <- NULL;fittedFM <- NULL;
  if(!is.null(uReg)){
    if(class(y)[1]=='fdata') y <- y$data
  }
  if(!is.null(time)){
    if(!is.null(fyList)) fyList$time <- time
    if(is.null(fyList)) fyList <- list(time=time)
    if(!is.null(uCoefList)) uCoefList <- lapply(uCoefList,function(i) 
      c(i,list(rtime=range(time))))
    if(is.null(uCoefList)) uCoefList <- list(list(rtime=range(time)))
    if(!is.null(fxCoefList)) fxCoefList <- lapply(fxCoefList,function(i) 
      c(i,list(rtime=range(time))))
    if(is.null(fxCoefList)) fxCoefList <- list(list(rtime=range(time)))
    if(!is.null(fxList)) fxList <- lapply(fxList,function(i) 
      c(i,list(rtime=range(time))))
    if(is.null(fxList)) fxList <- list(list(time=time))
  }
  
  if(class(y)[1]=='matrix'){
    ## define 'fd' object for y if y is a matrix
    y <- mat2fd(y,fyList)
  }
  if(class(y)[1]!='fd'){
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
    
  if(!is.null(uReg)){
    ## define list of x 
    if(class(uReg)[1]!='matrix') stop('class of uReg is expected to be matrix')
    x <- uReg
    nx <- ncol(x)
    lxList <- vector('list',length=nx)
    for(i in 1:nx) lxList[[i]] <- x[,i]
    
    ## define list of beta
    if(length(uCoefList)!=length(lxList)){
      cat('  The length of uCoefList is not equal to the number of scalar covariates entered in uReg.','\n')
      if(length(uCoefList)==0){
        cat('  Default uCoefList is applied.','\n')
        betalist <- lapply(lxList,function(i){
          i <- betaPar()
        })
      }
      
      if(length(uCoefList)>0){
        cat('  The first element of uCoefList will be applied to the functional regression coefficients of all scalar covariates.', '\n')
        betalist <- lapply(lxList,function(i){
          i <- betaPar(uCoefList[[1]])
        })
      }
    }
    if(length(uCoefList)==length(lxList))
      betalist <- lapply(uCoefList,betaPar)
    
    
    #regression
    ml <- fRegress(y, lxList, betalist)
    betaEstMat <- do.call('cbind',lapply(ml$betaestlist,function(i) 
      predict(i,y_time)))
    ml_fitted <- uReg%*%t(betaEstMat)
    
    if(class(response)[1]=='fd') y_raw <- eval.fd(y_time,response)
    if(class(response)[1]=='matrix'){
      if(nrow(response)==nrow(ml_fitted)) residML <- response-ml_fitted
      if(nrow(response)==ncol(ml_fitted)) residML <- t(response)-ml_fitted
    }
    y <- residML
    res <- c(res,list(y)); fittedFM <- c(fittedFM,list(ml_fitted))
  }
  mfTrainfd <- NULL
  ## functional response with functional covariates
  temp <- list(NULL)
  if(!is.null(fxReg)){
    y <- mat2fd(y,fyList)
    
    ## set up list of 'fd' object for x
    if(class(fxReg)[1]=='matrix' | class(fxReg)[1]=='fd')
      fxReg <- list(fxReg)
    if(class(fxReg)[1]=='list'){
      if(length(unique(unlist(lapply(fxReg,class))))!=1) 
        stop('functional covariates are expected to have same class')
      if(unique(unlist(lapply(fxReg,class)))=='matrix'){
        if(ncol(fxReg[[1]])!=length(y$fdnames$time)) fxReg <- lapply(fxReg,t)
        if(length(fxList)!=length(fxReg)){
          cat('     The length of fxList is not equal to the number of functional covariates entered in fxReg.','\n')
          fxReg <- lapply(fxReg,t)
          if(length(fxList)==0){
            cat('     Default fxList is applied','\n')
            fxReg <- lapply(fxReg,mat2fd)
          }
          if(length(fxList)>0){
            cat('     The first element of fxList will be applied to all functional covariates.','\n')
            
            fxReg <- lapply(fxReg,function(i){
              i <- mat2fd(i,fxList[[1]])
            })
          }
        }
        if(length(fxList)==length(fxReg)){
          fxReg <- lapply(1:length(fxReg),function(i){
            mat2fd(fxReg[[i]],fxList[[i]])
          })
        }
      }
      
      
      ## set up list of fdPar object for beta
      if(length(fxCoefList)!=length(fxReg)){
        cat('     The length of fxCoefList list is not equal to the number of functional covariates entered in fxReg.','\n')
        if(length(fxCoefList)==0){
          cat('     Default fxCoefList is applied','\n')
          if(1-concurrent){
            betalist <- lapply(fxReg,function(i) i=betaPar(list(bivar=TRUE)))
          }else{
            betalist <- lapply(fxReg,function(i) i=betaPar())
          }
        }
        if(length(fxCoefList)>0){
          cat('     The first element of fxCoefList will be applied to the functional regression coefficients of all functional covariates.','\n')
          if(1-concurrent) {
            betalist <- lapply(fxReg,function(i) i=betaPar(list(bivar=TRUE)))
          }else{
            betalist <- lapply(fxReg,function(i) i=betaPar(fxCoefList[[1]]))
            }
        }
      }
      if(length(fxCoefList)==length(fxReg)){
        if(1-concurrent) betalist <- lapply(fxCoefList, function(i) 
          betaPar(list(bivar=TRUE)))
        if(concurrent) betalist <- lapply(fxCoefList, betaPar)
      }
        
      ## regression
      temp <- NULL
      for(i in seq_along(fxReg)){
        if(is.matrix(y)) y <- mat2fd(y,fyList)
        x <- fxReg[[i]]
        if(concurrent){
          const <- rep(1,dim(x$coef)[2])
          xlist <- list(const=const,x=x)
          
          b1 <- betalist[[i]]
          bList <- list(const=b1,x=b1)
          mf <- fRegress(y,xlist,bList) 
          temp <- c(temp,list(mf))
          
          # evaluate functional coefficients
          betaEstMat <- list(predict(const=mf$betaestlist$const,y_time)[,1],
                             x=predict(const=mf$betaestlist$x,y_time)[,1])
          mf_fitted <- apply(t(eval.fd(y_time,x)),1,function(i) 
            i=i*betaEstMat[[2]]+betaEstMat[[1]])
          
          if(class(y)=='fd') y_raw <- t(eval.fd(y_time,y))
          if(class(y)=='matrix'){
            if(nrow(y)==nrow(mf_fitted)) residMF <- y_raw-mf_fitted
            if(nrow(y)==ncol(mf_fitted)) residMF <- t(y_raw)-mf_fitted
          }
          y <- residMF
          res <- c(res,list(y)); fittedFM <- c(fittedFM,list(mf_fitted))
        }
        
        if(1-concurrent){
          bList <- list(betaPar(), betalist[[i]])
          mf <- linmod(x,y,bList) 
          temp <- c(temp,list(mf))
          
          betaEstMat <- list(b0=eval.fd(y_time,mf$beta0estfd)[,1],
                             b1=eval.bifd(y_time,y_time,mf$beta1estbifd)[,1])
          mf_fitted <- apply(t(eval.fd(y_time,x)),1,function(i) 
            i=i%*%betaEstMat[[2]]/length(y_time)^2+betaEstMat[[1]])
          
          if(class(y)[1]=='fd') y_raw <- t(eval.fd(y_time,y))
          if(class(y)[1]=='matrix'){
            if(nrow(y)==nrow(mf_fitted)) residMF <- y_raw-mf_fitted
            if(nrow(y)==ncol(mf_fitted)) residMF <- t(y_raw)-mf_fitted
          }
          y <- residMF
          res <- c(res,list(y)); fittedFM <- c(fittedFM,list(mf_fitted))
        }
      }
      mfTrainfd <- fxReg
    }
  }
  
  
  out <- list(model=list(ml=ml,mf=temp),res=y,resList=res,
              'mfTrainfd'=mfTrainfd,fyl=fyList,'fittedFM'=fittedFM)
  return(out)
}


main3 <- function(response,uReg){
  ### this is for the case that the response are longitudinal data, while the 
  ### covariates are scalars. 
  
  y <- response
  ml <- NULL
  if(!is.null(uReg)){
    ml <- lm(y~uReg)
    y <- as.matrix(resid(ml))
  }

  out <- list(model=list(ml=ml,fl=NULL),res=y)
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
#'   \item{nbasis}{Number of basis functions used in smoothing, default to be
#'   less or equal to 23.} \item{norder}{Order of the functional curves default
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
mat2fd <- function(mat,fdList=NULL){
  fl <- list(time=seq(0,1,len=ncol(mat)),nbasis=min(as.integer(ncol(mat)/5),23),
             norder=6,bSpline=TRUE,Pen=c(0,0),lambda=1e-4)
  nbasis <- c(fdList$nbasis,fl$nbasis)[1]
  norder <- c(fdList$norder,fl$norder)[1]
  lambda <- c(fdList$lambda,fl$norder)[1]
  bSpline <- c(fdList$bSpline,fl$bSpline)[1]
  time <- list(a=fdList$time,b=fl$time)
  time <- time[[which(unlist(lapply(time,is.null))^2==0)[1]]]
  if(1-bSpline) fl$Pen <- c(c(0,(2*pi/diff(range(time)))^2,0))
  Pen <- list(a=fdList$Pen,b=fl$Pen)
  Pen <- Pen[[which(unlist(lapply(Pen,is.null))^2==0)[1]]]
  if(bSpline)  basis <- create.bspline.basis(range(time),nbasis,norder)
  if(1-bSpline) basis <- create.fourier.basis(range(time),nbasis,
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
  if(1-bSpline) fl$Pen <- c(c(0,(2*pi/diff(rtime))^2,0))
  Pen <- list(a=betaList$Pen,b=bl$Pen)
  Pen <- Pen[[which(unlist(lapply(Pen,is.null))^2==0)[1]]]
  if(bSpline)  basis <- create.bspline.basis(rtime,nbasis,norder)
  if(1-bSpline) basis <- create.fourier.basis(rtime,nbasis,diff(rtime))
  Par <- vec2Lfd(Pen,rtime)
  if(1-bivar){
    betaPar <- fdPar(basis, Par, lambda)
  }else{
    lambdas <- c(betaList$lambdas,bl$lambdas)
    betaPar <- bifdPar(bifd(matrix(0,nbasis,nbasis), basis, basis),
                    Par, Par, lambda, lambdas)
  } 
  return(betaPar)  
}


repgp.loglikelihood <- function(hyper.p,response,Data,Cov,gamma=1,nu=1.5,
                                time=NULL,...){
  ### response is expected to be matrices with ncol replications and nrow 
  ### observations
  ### Data is expected to be a list with matrices
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


fisherinfo <- function(pp.cg,X,Y,Cov,gamma,nu){
  n <- length(Cov)
  CovList <- vector('list',n)
  for(i in 1:n) CovList[i] <- list(paste0('cov.',Cov[i]))
  CovL <- lapply(CovList,function(j){
    f <- get(j)
    if(j=='cov.pow.ex')
      return(f(pp.cg,X,X,gamma=gamma))
    if(j=='cov.matern')
      return(f(pp.cg,X,X,nu=nu))
    if(!(j%in%c('cov.pow.ex', 'cov.matern')))
      return(f(pp.cg,X,X))
  }  )
  if(length(CovL)==1)
    Q <- CovL[[1]]
  if(length(CovL)>1)
    Q <- Reduce('+',CovL)
  
  response <- as.matrix(Y)
  X <- as.matrix(X)
  Q <- Q+diag(exp(pp.cg$vv),dim(Q)[1])
  invQ <- chol2inv(chol(Q))
  
  QR <- invQ%*%response
  AlphaQ <- QR%*%t(QR)-invQ
  
  
  D2fx <- lapply(seq_along(pp.cg),function(i){
    Dp <- pp.cg[i]
    name.Dp <- names(Dp)
    f <- get(paste0('D2',name.Dp))
    if(name.Dp%in%c('pow.ex.w','pow.ex.v') )
      D2para <- f(pp.cg,X,gamma=gamma,inv.Q=invQ,Alpha.Q=AlphaQ)
    if(name.Dp%in%c('matern.w','matern.v') )
      D2para <- f(pp.cg,X,nu=nu,inv.Q=invQ,Alpha.Q=AlphaQ)
    if(name.Dp%in%c('linear.a'))
      D2para <- f(hyper=pp.cg, input=X, Alpha.Q=AlphaQ)
    if(name.Dp%in%c('linear.i'))
      D2para <- f(hyper=pp.cg, inv.Q=invQ, Alpha.Q=AlphaQ)
    if(!name.Dp%in%c('pow.ex.w','pow.ex.v','matern.w','matern.v',
                     'linear.a','linear.i'))
      D2para <- f(hyper=pp.cg, input=X, inv.Q=invQ, Alpha.Q=AlphaQ)
    return(D2para)
  })
  names(D2fx) <- names(pp.cg)
  II <- abs(-1/(unlist(D2fx)*dim(X)[1]))
  return(II)
}

#' @importFrom  fda.usc is.fdata
#' @importFrom  fda is.fd
gpfrtrain <- function(response,uReg=NULL,fxReg=NULL,fyList=NULL,uCoefList=NULL,
                      fxList=NULL,fxCoefListScalarResp=NULL,concurrent=TRUE,
                      fxCoefList=NULL,gpReg=NULL,hyper=NULL,
                      Cov,gamma=2,nu=1.5,useGradient=T,time=NULL,
                      rel.tol=1e-10,
                      trace.iter=5,fitting=FALSE){
  y <- y_raw <- response
  if(is.vector(y)) model <- 'main1'
  if(is.matrix(y)){
    if(ncol(y)==1 | nrow(y)==1) model <- 'main1'
    else model <- 'main2'
  }
  if(is.fd(y) | is.fdata(y)) model  <- 'main2'
  
  if(is.data.frame(y)) model <- 'main3'
  ModelType <- model
  
  if(model=='main1'){
    model <- main1(response=response,uReg=uReg,fxReg=fxReg,fxList=fxList,
                   fxCoefListScalarResp=fxCoefListScalarResp)
  }
  if(model=='main2'){
    model <- main2(response=response,uReg=uReg,fxReg=fxReg,fyList=fyList,
                   uCoefList=uCoefList,fxList=fxList,concurrent=concurrent,
                   fxCoefList=fxCoefList,time=time)  
  }
  iuuL <- NULL
  iuuL <- chol2inv(chol(crossprod(cbind(uReg))))
  
  fittedFM <- model$fittedFM
  
  ## convert fd/fdata class to matrix
  response <- model$res; Data <- gpReg
  if(class(response)[1]!='matrix') stop('expecting matrix residual from 
                                        functinoal regression')
  ftime <- model$fyl$time
  if(!is.null(ftime)) time <- ftime
  if(is.null(ftime) & is.null(time)) stop('expecting input time')
  
  
  
  if(unique(unlist(lapply(Data,class)))[1]=='fdata'){
    Data <- lapply(Data,function(i) i=t(i$data))
  }

  if(class(response)[1]=='matrix') response <- t(response)
  if(unique(unlist(lapply(Data,class))[1]=='matrix')){
    Data <- lapply(Data,t)
  }
  
  if(class(Data)=='fd')
    Data <- (eval.fd(time,Data))
  if(unique(unlist(lapply(Data,class))=='fd')){
    Data <- lapply(Data,function(i) (eval.fd(time,i)))
  }
  
  
  ### this is an approximation of iuu for functional regression
  iuuF <- NULL
  if(ModelType=='main2' & !is.null(fxReg)){
    if(concurrent){
      iuuF <- lapply(model$model$mf,function(i){
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
 
  
  stepsize <- nrow(response)  #training data size
  tdsize <- as.integer(stepsize/2)  #choose half data for training
  
  
  sample_idx <- 1:tdsize*2-1
  response_raw <- response; 
  Data_raw <- Data ## keep original data before reducing the dimension
  response <- response[sample_idx,]
  if(ncol(Data[[1]])!=ncol(response)){
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
  

  cat('   Start optimization: ','\n')

  if(!useGradient){repgp.Dloglikelihood <- NULL}
  pp <- nlminb(init0,repgp.loglikelihood,repgp.Dloglikelihood,
               response=response,Data=Data,Cov=Cov,gamma=gamma,nu=nu, 
               control=list(trace=trace.iter,rel.tol=rel.tol))

  cat('    Optimization has been finished.','\n','\n')
  pp.cg <- pp[[1]]
  
  names(pp.cg) <- names(init0)
  pp.df <- data.frame(pp.cg=pp.cg,pp.N=substr(names(init0),1,8))
  names(pp.df) <- c('pp.cg','pp.N')
  pp.cg <- split(pp.df$pp.cg,pp.df$pp.N)
  
  allbat <- vector('list',length=ncol(response))
  for(i in 1:ncol(response)){
    allbat[[i]] <- cbind(response[,i],do.call('cbind',
                                              lapply(Data,function(j) j=j[,i])))
  }
  
  Qlist <- lapply(allbat,function(l) fisherinfo(pp.cg=pp.cg,
                                                X=as.matrix(l[,-1]),
                                                Y=l[,1],Cov=Cov,
                                                gamma=gamma,nu=nu))
  II <- abs(-1/apply(do.call('cbind',Qlist),1,sum))
  # cat("Fisher's information done",'\n','\n')
  
  mean <- t(Reduce('+',fittedFM)) ## mean from functional regression model
  
  fitted <- fitted.sd <- NULL
  n <- length(Cov)
  hyper.cg <- pp.cg
  if(fitting){
    fitted <- fitted.sd <- matrix(0,ncol=ncol(response_raw),
                                  nrow=nrow(response_raw))
    for(i in seq_along(allbat)){
      dr <- as.matrix(do.call('cbind',lapply(Data_raw,
                                             function(j) j=j[,i])))
      yy <- as.matrix(response_raw[,i])
      
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
    
    cat(paste0('GPFR model fitted to all the curves.','\n'))
  }
    
  result <- list('hyper'=pp.cg,'I'=II, 'modellist'=model$model,'CovFun'=Cov,
                 'gamma'=gamma,'nu'=nu,'init_resp'=y_raw,meanFM=mean,
                 'resid_resp'=response_raw,'fitted.mean'=fitted,
                 'fitted.sd'=fitted.sd,'ModelType'=ModelType,'lTrain'=uReg,
                 'fTrain'=fxReg,'mfTrainfd'=model$mfTrainfd,'gpTrain'=Data_raw,
                 'time'=time,'iuuL'=iuuL,'iuuF'=iuuF,
                 'fittedFM'=fittedFM,'fyList'=fyList)
  class(result) <- 'gpfr'
  return(result)
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
#' @param fyList A list to control the smoothing of response.
#' @param uCoefList A list to control the smoothing of the regression
#'   coefficient function of the scalar covariates in the FR model.
#' @param fxList A list to control the smoothing of functional covariates in the
#'   FR model.
#' @param concurrent Logical. If TRUE (default), concurrent functional
#'   regression will be carried out; otherwise, the full functional regression
#'   will be carried out.
#' @param fxCoefList A list to control the smoothing of the regression
#'   coefficient function of functional covariates in the functional concurrent
#'   model.
### @param fxCoefListScalarResp A list to control the smoothing of functional covariates in
###   the FR model with scalar response and functional covariates.
#' @param gpReg Covariates in the GP model. It should be a matrix, a numeric vector, 
#' an 'fd' object, a list of matrices or a list of 'fd' objects.
#' @param hyper Vector of initial hyperparameters. Default to NULL.
#' @param NewHyper Vector of names of new hyperparameters from the customized
#'   kernel function.
#' @param Cov Covariance function(s) to use. Options are: 'linear', 'pow.ex',
#'   'rat.qu', and 'matern'. Default to 'power.ex'.
#' @param gamma Power parameter used in powered exponential kernel function. It
#'   must be 0<gamma<=2.
#' @param nu Smoothness parameter of the Matern class. It must be a positive
#'   value.
#' @param useGradient Logical. If TRUE, first derivatives will be used in the
#'   optimization.
#' @param rel.tol Relative tolerance passed to nlminb(). Default to be 1e-10.
#' @param trace.iter Print the processing of iterations of optimization.
#' @param fitting Logical. If TRUE, fitting is carried out. Default to FALSE.
#'
#' @details \code{fyList} is a list with the following items:
#'
#'   \itemize{ \item \code{time}: a sequence of time points; default to be 100
#'   points from 0 to 1. \item \code{nbasis}: number of basis functions used in
#'   smoothing, default to be less than or equal to 23. \item \code{norder}:
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
#'   \itemize{ \item \code{rtime}: range of time, default to be c(0,1). \item
#'   \code{nbasis}: nnumber of basis functions used in smoothing, default to be
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
#' @return A list containing: \describe{ \item{hyper}{Estimated hyperparameters}
#'   \item{I}{A vector of estimated standard deviation of hyperparameters}
#'   \item{modellist}{List of FR models fitted before Gaussian process}
#'   \item{CovFun}{Covariance function used} \item{gamma}{Parameter 'gamma' used
#'   in Gaussian process with powered exponential kernel} \item{nu}{Parameter
#'   'nu' used in Gaussian process with Matern kernel} \item{init_resp}{Raw
#'   response data} \item{resid_resp}{Residual after the fitted values from FR
#'   models have been taken out} \item{fitted}{Fitted values}
#'   \item{fitted.sd}{Standard deviation of the fitted values}
#'   \item{ModelType}{The type of the model applied in the function.}
#'   \item{lTrain}{Training scalar covariates for the FR model}
#'   \item{fTrain}{Training functional covariates for the FR model}
#'   \item{mfTrainfd}{List of 'fd' objects from training data for FR model with
#'   functional covariates} \item{gpTrain}{Training data for Gaussian Process}
#'   \item{time}{Input time 't'} \item{iuuL}{Inverse of covariance matrix for
#'   uReg} \item{iuuF}{Inverse of covariance matrix for fxReg}
#'   \item{fittedFM}{Fitted values from the FR model} \item{fyList}{fyList
#'   object used} }
#' @export
#'
#' @examples
#' ## See examples in vignette:
#' # vignette("gpfr", package = "GPFDA")
gpfr <- function(response, time=NULL, uReg=NULL, fxReg=NULL,
                 fyList=NULL, uCoefList=NULL, 
                 fxList=NULL, concurrent=TRUE, fxCoefList=NULL,
                 gpReg=NULL, hyper=NULL, NewHyper=NULL, Cov='pow.ex', gamma=2, nu=1.5, 
                 useGradient=T, rel.tol=1e-10, trace.iter=5, fitting=FALSE){
  
  
  if(!is.matrix(response) & !is.fd(response)){
    stop("'response' must be either a matrix or an 'fd' object.")
  }
  fxCoefListScalarResp <- NULL # models with scalar response not implemented yet
  
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
  
  a1 <- gpfrtrain(response=response,uReg=uReg,fxReg=fxReg,gpReg=gpReg,
                  fyList=fyList,uCoefList=uCoefList,fxList=fxList,
                  fxCoefList=fxCoefList,
                  fxCoefListScalarResp=fxCoefListScalarResp,
                  hyper=hyper,
                  Cov=Cov,gamma=gamma,nu=nu,useGradient=useGradient,
                  fitting=fitting,time=time,
                  rel.tol=rel.tol,trace.iter=trace.iter,
                  concurrent=concurrent)
  return(a1)
}


#' Prediction of GPFR model
#'
#' Make predictions for test input data based on the GPFR model learnt by the
#' 'gpfr' function. Both Type I and Type II predictions can be made.
#'
#' @param train An object of class 'gpfr' obtained by the the 'gpfr' function.
#' @param testInputGP Test input data for the GP prediction. It must be a numeric 
#' vector, a matrix or an 'fd' object.
#' @param testTime Test time points for prediction. If NULL, default settings
#'   will be applied.
#' @param uReg Scalar covariates data of a new batch for the FR model.
#' @param fxReg Functional covariates data of a new batch for the FR model.
#' @param gpReg Input data for the GP part used for Type I prediction. It must
#'   be a list of three items. The names of the items must be 'response',
#'   'input', and 'time'. The item 'response' is the observed response for a new
#'   batch; 'input' is the observed functional covariates for a new
#'   batch,;'time' is the observed time for the previous two. If NULL (default),
#'   Type II prediction is carried out.
#' @param GPpredict Logical. If TRUE (default), GPFR prediction is carried out;
#'   otherwise only predictions based on the FR model is carried out.
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
gpfrPredict <- function(train, testInputGP, testTime=NULL, uReg=NULL, fxReg=NULL,
                     gpReg=NULL, GPpredict=TRUE){
  if(class(train)!='gpfr'){
    stop("The argument 'train' is expected to be an object of class 'gpfr' ",'\n')
  }
  
  if(is.null(gpReg)){
    predictionType <- 2
  }else{
    predictionType <- 1
  }
  model <- train$modellist
  if(is.null(model$ml) & !is.null(uReg)){
    cat('    model with scalar variable is not found, ignoring uReg','\n')
    uReg <- NULL
  }
  if(is.null(model$mf[[1]]) & !is.null(fxReg)){
    cat('    model with functional variable is not found, ignoring uReg','\n')
    fxReg <- NULL
  }
  if(!is.null(model$ml) & is.null(uReg) & train$ModelType=='main1'){
    stop('    expecting input variable for model with scalar variable. ','\n')
    uReg <- NULL
  }
  if(!is.null(model$ml) & is.null(uReg) & train$ModelType=='main2'){
    stop('    expecting input variable for model with scalar variable. ','\n')
    uReg <- NULL
  }
  if(!is.null(model$mf[[1]]) & is.null(fxReg)){
    stop('    expecting input variable for model with functional variable. ',
         '\n')
    fxReg <- NULL
  }
  
  if(!is.null(model$ml) | !is.null(model$mf))
    rtime <- c(model$ml$yhatfdobj$fd$basis$rangeval,
               model$mf[[1]]$yhatfdobj$basis$rangeval,
            model$mf[[1]]$yhatfdobj$fd$basis$rangeval)[1:2]
  if(is.null(model$ml) & is.null(model$mf)) rtime <- c(0,1)
  
    if(class(testInputGP)=='numeric'){
    testInputGP <- as.matrix(testInputGP)
  }
  if(is.matrix(testInputGP)){
    test <- testInputGP
    test <- t(test)
    if(is.null(testTime)) time <- seq(rtime[1],rtime[2],len=col(test))
    if(!is.null(testTime)) time <- testTime
  }
  
  if(class(testInputGP)[1]=='fd'){
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
  
  if(!is.null(gpReg)& class(gpReg)!='list'){
    cat('Type I prediction is expecting gpReg to be a list with a response and 
        an input. Do Type II prediction instead')
    gpReg <- NULL
  }
  type <- 2
  if(!is.null(gpReg) & class(gpReg)=='list'){
    type <- 1
    nl <- names(gpReg)
    if(sum(c('response','input','time')%in%names(gpReg))!=3)
      stop('check the name of the gpReg list. there must be one "response", 
           one "input" and one "time"')
    if(!1%in%dim(as.matrix(gpReg$response)))
      stop('check the diemsion of the "response", it must be a column or 
           row matrix, or a vector')
    gplReg <- uReg
    gpfReg <- gpReg$input
    gpresp <- gpReg$response
    gptime <- gpReg$time
    yhat_ml_list <- vector('list',length=2)
    if(!is.null(uReg)) uRegList <- list(f=uReg,gp=gplReg)
    if(!is.null(time)) timeList <- list(f=time,gp=gptime)
    if(!is.null(fxReg)) fRegList <- list(f=fxReg,gp=gpfReg)
    # else stop('expecting Test input')
  }
  
  ### predict ml for main2
  if(is.null(uRegList)) yhat_ml <- 0
  if(!is.null(uRegList)){
    for(ii in seq_along(uRegList)){
      uReg <- uRegList[[ii]]
      time <- timeList[[ii]]
      
      if (train$ModelType=='main2'){
        if(!is.null(uReg)){
          if(is.vector(uReg)) uReg <- t(as.matrix(uReg))
          if(is.matrix(uReg)){
            if(length(model$ml$betaestlist)!=ncol(uReg))
              stop('dimension of uReg does not match the model')
            if(length(model$ml$betaestlist)==ncol(uReg)){
              x <- model$ml$betaestlist
              for(i in 1:ncol(uReg))
                x[[i]] <- as.matrix(uReg[,i])
            }
          }
          if(class(uReg)[1]=='list'){
            if(unique(unlist(lapply(uReg,class)))=='fd')
              x <- lapply(uReg,function(i) t(i$coefs))
            if(unique(unlist(lapply(uReg,class)))=='matrix')
              x <- uReg
          }
          betalist <- lapply(model$ml$betaestlist,function(i){
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
  }

  ### predict mf for main2
  if(is.null(fRegList)) yhat_mf <- gpyhat_mf <- matrix(0,ncol=1,nrow=1)
  if(!is.null(fRegList)){
    for(ii in seq_along(fRegList)){
      fxReg <- fRegList[[ii]]
      time <- timeList[[ii]]
      if (train$ModelType=='main2'){
        if(!is.null(fxReg)){
          if(is.fdata(fxReg)) fxReg <- list(t(fxReg$data))
          if(is.matrix(fxReg) | is.fd(fxReg)) fxReg <- list((fxReg))
          if(unique(unlist(lapply(fxReg,class)))=='fd'){
            fxReg <- lapply(fxReg,function(i) eval.fd(time,i))
          }
          if(unique(unlist(lapply(fxReg,class)))=='matrix'){
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
          if(!'bifd'%in%unlist(lapply(model$mf[[1]],class))){
            fbeta <- lapply(model$mf,function(i){
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
                    len=train$modellist$mf[[1]]$betaestlist$x$fd$basis$nbasis),
                             smooth.basis(iuuTime,fx,Basis)$fd))
                mf_var[i] <- fx%*%train$iuuF[[i]]%*%t(fx)
              }
            }
            fxReg[[i]] <- apply(fxReg[[i]],2,function(j){
              j=as.matrix(fbeta[[i]][,1])+as.matrix(fbeta[[i]][,2])*j
            })
            if(ii==1) yhat_mf <- Reduce('+',fxReg)
            if(ii==2) gpyhat_mf <- Reduce('+',fxReg)
          }
          
          if('bifd'%in%unlist(lapply(model$mf[[1]],class))){
            fbeta <- lapply(model$mf,function(i){
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
  }
#   
  f.mean <- yhat_ml+apply(yhat_mf,1,sum)
  f.var <- apply((train$init_resp-t(train$resid_resp))^2,2,sum)/(nrow(
    train$init_resp)-1)
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
    y_gppred <- gprPredict(train=FALSE,hyper=train$hyper, 
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

    fitted <- matrix(0,ncol=ncol(train$resid_resp),nrow=ncol(test))
    fitted.var <- fitted
    
    for(i in 1:ncol(train$resid_resp)){
      input <- do.call('cbind',lapply(train$gpTrain,function(j) j=j[,i]))
      y_gppred <- gprPredict(train=FALSE,inputNew=t(as.matrix(test)),
                            hyper=train$hyper,input=as.matrix(input), 
                            Y=as.matrix(train$resid_resp[,i]), 
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
  
  
  # ypredup <- ypred + 1.96*sqrt(s2)
  # ypredlo <- ypred - 1.96*sqrt(s2)
  # 
  # CI <- cbind(ypred, ypredup, ypredlo)
  
  # result <- c(list(ypred=CI, testtime=time, predtime=testtime,
  #                  ypred.mean=ypred,
  #                  ypred.sd=sqrt(s2)),predictionType=predictionType,
  #             unclass(train))
  
  result <- c(list(ypred.mean=ypred,
                 ypred.sd=sqrt(s2),
                 testTime=testtime, 
                 predictionType=predictionType),
                 unclass(train))
  
  }
  

  class(result) <- 'gpfr'
  return(result)
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
#' @param mar Graphical parameter passed to par().
#' @param oma Graphical parameter passed to par().
#' @param cex.lab Graphical parameter passed to par().
#' @param cex.axis Graphical parameter passed to par().
#' @param cex.main Graphical parameter passed to par().
#' @param ... Other graphical parameters passed to plot().
#' @importFrom graphics polygon
#' @importFrom graphics matpoints
#' @importFrom graphics matlines
#' @importFrom graphics lines
#' @importFrom grDevices rgb
#' @importFrom fda matplot
#' @return A plot.
#' @export
#'
#' @examples
#' ## See examples in vignette:
#' # vignette("gpfr", package = "GPFDA")
plot.gpfr <- function (x, type=c('raw','meanFunction','fitted','prediction'),
                       ylab='y', xlab='t', ylim=NULL, realisations=NULL, 
                       alpha=0.05,
                       colourTrain=2, colourNew=4,
                       mar=c(4.5,5.1,2.2,0.8), oma=c(0,0,1,0), 
                       cex.lab=1.5, cex.axis=1, cex.main=1.5, ...){
  obj <- x
  if(!type%in%c('raw','meanFunction','fitted','prediction')){
    stop('type must be one of the raw, fitted or prediction')
  }
  
  z <- stats::qnorm(1-alpha/2)
  
  if(is.null(realisations)){
    idx <- 1:ncol(obj$fitted.mean)
  }else{
    idx <- realisations
  }
  numRealis <- length(idx)
  
  old <- par(mar=mar, oma=oma, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)
  
  if(type=='raw'){
    if(is.null(ylim)){
      ylim <- range(obj$init_resp[idx,])
    }
    matplot(obj$time,t(obj$init_resp[idx,,drop=F]),type='l',lwd=2,lty=3,
            main='Raw data',xlab=xlab,ylab=ylab, ylim=ylim)
    if(numRealis<6){
      matpoints(obj$time,t(obj$init_resp[idx,,drop=F]),pch=4,cex=1,lty=3)
    }
  }
  
  if(type=='meanFunction'){
    if(is.null(ylim)){
      ylim <- range(c(obj$init_resp[idx,], c(obj$modellist$ml$yhatfdobj$y[,idx])))
    }
    matplot(obj$time,t(obj$init_resp[idx,,drop=F]),type='l',lwd=2,lty=3,
            main='Mean function fit', xlab=xlab, ylab=ylab, ylim=ylim)
    if(numRealis<6){
      matpoints(obj$time,t(obj$init_resp[idx,,drop=F]),pch=4,cex=1,lty=3)
    }
    for(i in 1:numRealis){
      ii <- idx[i]
      lines(obj$modellist$ml$yhatfdobj$argvals[,1],
            obj$modellist$ml$yhatfdobj$y[,ii], lwd=2, lty=1, col=ii)
    }
    
    # lines(obj$modellist$ml$yhatfdobj, lwd=2, lty=1)
    # matplot(obj$time,t(obj$init_resp[idx,,drop=F]),type='p',pch=4,lty=3,cex=1,lwd=1,
    #         main='Mean function',xlab=xlab,ylab=ylab)
    # lines(obj$modellist$ml$yhatfdobj, lwd=2)
  }
  
  if(type=='fitted'){
    if(is.null(ylim)){
      ylim <- range(c(obj$init_resp[idx,], 
                      c((obj$fitted.mean-obj$fitted.sd*z)[,idx]),
                      c((obj$fitted.mean+obj$fitted.sd*z)[,idx])
                    ))
    }
    matplot(obj$time,t(obj$init_resp[idx,,drop=F]),type='p',pch=4,lty=3,cex=0,lwd=1,
            main='GPFR fit',xlab=xlab,ylab=ylab, ylim=ylim)
    for(i in 1:numRealis){
      ii <- idx[i]
      polygon(c(obj$time, rev(obj$time)), 
              c((obj$fitted.mean-obj$fitted.sd*z)[,ii], 
                rev((obj$fitted.mean+obj$fitted.sd*z)[,ii])), 
              col = rgb(127,127,127,80, maxColorValue = 255), border = NA)
    }
    matlines(obj$time,obj$fitted.mean[,idx,drop=F],type='l',pch=4,lwd=1.5, lty=1)
    matpoints(obj$time,t(obj$init_resp[idx,,drop=F]),pch=4,lty=3,cex=1,lwd=1)
    
  }
  
  if(type=='prediction'){
    
    if(obj$predictionType==1){
      main <- "Type I prediction"
    }else{
      main <- "Type II prediction"
    }
    if(is.null(ylim)){
      if(identical(idx, 0)){
        ylim <- range(c(obj$ypred.mean - z*obj$ypred.sd),
                      c(obj$ypred.mean + z*obj$ypred.sd))
        # ylim <- range(c(obj$ypred[,2], obj$ypred[,3], obj$ypred[,1]))
      }else{
        ylim <- range(c(obj$init_resp[idx,]),
                      c(obj$ypred.mean - z*obj$ypred.sd),
                      c(obj$ypred.mean + z*obj$ypred.sd))
      }
    }
    if(identical(idx, 0)){
      plot(NA, main=main,xlab=xlab,ylab=ylab, ylim=ylim, 
           xlim = range(obj$time))
    }else{
      matplot(obj$time,t(obj$init_resp[idx,,drop=F]),type='l',lwd=2,lty=3,col=colourTrain,
              main=main,xlab=xlab,ylab=ylab, ylim=ylim)
        if(numRealis<6){
          matpoints(obj$time,t(obj$init_resp[idx,,drop=F]),pch=4,cex=1,lty=3,col=colourTrain)
        }
    }
    
    polygon(c(obj$testTime, rev(obj$testTime)), 
            c(obj$ypred.mean + z*obj$ypred.sd, rev(c(obj$ypred.mean - z*obj$ypred.sd))), 
            col = rgb(127,127,127,100, maxColorValue = 255), border = NA)
    lines(obj$testTime,obj$ypred.mean[,1],col=colourNew,lwd=2)
    
  }
  
  par(old)
  
}


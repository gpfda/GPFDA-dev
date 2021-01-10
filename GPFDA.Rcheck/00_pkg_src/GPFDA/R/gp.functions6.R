
######################## Covariance functions ############################ 

#' @title Calculate a covariance matrix
#' @description Evaluates one of the following covariance functions at input 
#' vectors t and t': 
#' \itemize{
#' \item Powered exponential
#' \item Rational quadratic
#' \item Matern
#' \item Linear
#' }
#'
#' @param hyper The hyperparameters. It must be a list with certain names. See 
#' details.
#' @param input The covariate t. It must be either a matrix, where each column
#'   represents a covariate, or a vector if there is only one covariate.
#' @param input.new The covariate t'. It also must be a vector or a matrix. 
#' If NULL (default), 'input.new' will be set to be equal to 'input' and the 
#' function will return a squared, symmetric covariance matrix.
#' @param gamma Power parameter used in powered exponential kernel function. It
#'   must be 0<gamma<=2. Default to 2, which gives the squared exponential 
#'   covariance function.
#' @param nu Smoothness parameter of the Matern class. It must be a positive
#'   value.
#'
#' @details The names for the hyperparameters should be: 
#' \itemize{
#' \item "pow.ex.v" and "pow.ex.w" (powered exponential);
#' \item "rat.qu.v", "rat.qu.w" and "rat.qu.a" (rational quadratic);
#' \item "matern.v" and "matern.w" (Matern);
#' \item "linear.i" and "linear.a" (linear);
#' \item "vv" (Gaussian white noise).
#' }
#' @references Shi, J. Q., and Choi, T. (2011), ``Gaussian Process Regression
#'   Analysis for Functional input'', CRC Press.
#' @return A covariance matrix
#' @name covMat
NULL


#' @rdname covMat
#' @export
cov.pow.ex <- function(hyper, input, input.new=NULL, gamma=2){
  
  hyper <- lapply(hyper,exp)
  input <- as.matrix(input)
  Q <- ncol(input)
  if(Q==1){
    A <- as.matrix(hyper$pow.ex.w)
  }else{
    A <- diag(hyper$pow.ex.w)
  }
  
  if (is.null(input.new)) {
    distMat <- distMatSq(input=input, A=A, power=gamma)
    exp.v.power <- hyper$pow.ex.v*exp(-distMat)
  }else{
    input.new <- as.matrix(input.new)
    distMat <- distMat(input=input, inputNew=input.new, A=A, power=gamma)
    exp.v.power <- hyper$pow.ex.v*exp(-distMat)
  }

  return(exp.v.power)
}


#' @rdname covMat
#' @export
cov.rat.qu <- function(hyper, input, input.new=NULL){
  
  hyper <- lapply(hyper,exp)
  input <- as.matrix(input)
  Q <- ncol(input)
  if(Q==1){
    A <- as.matrix(hyper$rat.qu.w)
  }else{
    A <- diag(hyper$rat.qu.w)
  }
  
  if (is.null(input.new)) {
    distMat <- distMatSq(input=input, A=A, power=2)
    covratqu <- hyper$rat.qu.v*(1 + distMat)^(-hyper$rat.qu.a)
  }else{
    input.new <- as.matrix(input.new)
    distMat <- distMat(input=input, inputNew=input.new, A=A, power=2)
    covratqu <- hyper$rat.qu.v*(1 + distMat)^(-hyper$rat.qu.a)
  }
  
  return(covratqu)
}

#' @rdname covMat
#' @export
cov.matern <- function(hyper, input, input.new=NULL, nu){
  
  input <- as.matrix(input)
  cc <- exp(hyper$matern.v)
  Q <- ncol(input)
  if(Q==1){
    A <- as.matrix(exp(hyper$matern.w))
  }else{
    A <- diag((exp(hyper$matern.w)))
  }
  
  if (is.null(input.new)) {
    
    if(nu==3/2){
      distmat <- sqrt(distMatSq(input=input, A=A, power=2))
      covmatern <- cc*(1+sqrt(3)*distmat)*exp(-sqrt(3)*distmat)
    }else{
      if(nu==5/2){
        distmat <- sqrt(distMatSq(input=input, A=A, power=2))
        covmatern <- cc*(1+sqrt(5)*distmat+(5/3)*distmat^2)*
          exp(-sqrt(5)*distmat)
      }else{
        covmatern <- CovMaternCppSq(input=input, cc=cc, A=A, nu=nu)
      }
    }

  }else{
    input.new <- as.matrix(input.new)
    if(nu==3/2){
      distmat <- sqrt(distMat(input=input, inputNew=input.new, A=A, power=2))
      covmatern <- cc*(1+sqrt(3)*distmat)*exp(-sqrt(3)*distmat)
    }else{
      if(nu==5/2){
        distmat <- sqrt(distMat(input=input, inputNew=input.new, A=A, power=2))
        covmatern <- cc*(1+sqrt(5)*distmat+(5/3)*distmat^2)*
          exp(-sqrt(5)*distmat)
      }else{
        covmatern <- CovMaternCpp(input=input, inputNew=input.new, cc=cc, A=A, nu=nu)   
      }
    }
  }
  
  return(covmatern)
}


#' @rdname covMat
#' @export
cov.linear <- function(hyper, input, input.new=NULL){
  
  hyper <- lapply(hyper,exp)
  input <- as.matrix(input)
  Q <- ncol(input)
  if(Q==1){
    A <- as.matrix(hyper$linear.a)
  }else{
    A <- diag(hyper$linear.a)
  }
  
  if (is.null(input.new)) {
    cov.lin <- distMatLinearSq(input=input, A=A)
  }else{
    input.new <- as.matrix(input.new)
    cov.lin <- distMatLinear(input=input, inputNew=input.new, A=A)
  }
  
  cov.lin <- hyper$linear.i + cov.lin
  
  return(cov.lin)
}


#' Gaussian Process regression
#'
#' Gaussian Process regression for a single or multiple independent
#' realisations.
#'
#' @param input Input covariates. It must be either a matrix, where each column
#'   represents a covariate, or a vector if there is only one covariate.
#' @param response Response data. It should be a matrix, where each column is a
#'   realisation. It can be a vector if there is only one realisation.
#' @param Cov Covariance function(s) to use. Options are: 'linear', 'pow.ex',
#'   'rat.qu', and 'matern'. Default to 'power.ex'.
#' @param m If Subset of Data is to be used, m denotes the subset size and
#'   cannot be larger than the total sample size. Default to NULL.
#' @param hyper The hyperparameters. Default to NULL. If not NULL, then it must
#'   be a list with appropriate names.
#' @param NewHyper Vector of names of the new hyperparameters of the customized
#'   kernel function. These names must have the format: xxxxxx.x, i.e. '6 digit'
#'   followed by 'a dot' followed by '1 digit'. This is required for both
#'   'hyper' and 'NewHyper'
#' @param meanModel Type of mean function. It can be \describe{ \item{0}{Zero
#'   mean function} \item{1}{Constant mean function to be estimated}
#'   \item{'t'}{Linear model for the mean function} \item{'avg'}{The average
#'   across replications is used as the mean function. This is only used if
#'   there are more than two realisations observed at the same input coordinate
#'   values.} } Default to 0. If argument 'mu' is specified, then 'meanModel'
#'   will be set to 'userDefined'.
#' @param mu Mean function specified by the user. It must be a vector. Its
#'   length must be the same as the sample size, that is, nrow(response).
#' @param gamma Power parameter used in powered exponential kernel function. It
#'   must be 0<gamma<=2.
#' @param nu Smoothness parameter of the Matern class. It must be a positive
#'   value.
#' @param useGradient Logical. If TRUE, first derivatives will be used in the
#'   optimization.
#' @param iter.max Maximum number of iterations allowed. Default to 100. If
#'   'rel.tol' is reduced, then the number of iterations needed will be less.
#' @param rel.tol Relative convergence tolerance. Default to 8e-10. Smaller
#'   rel.tol means higher accuracy and more time to converge.
#' @param trace The value of the objective function and the parameters is
#'   printed every trace'th iteration. Defaults to 0 which indicates no trace
#'   information is to be printed.
#' @param nInitCandidates Number of initial hyperparameter vectors. The
#'   optimization starts with the best.
#'
#' @details The most important function of the package. It fits the GPR model
#'   and stores everything necessary for prediction. The optimization used in
#'   the function is 'nlminb'. The names for the hyperparameters should be:
#'   "linear.a" for linear covariance function, "pow.ex.w", "pow.ex.v" for power
#'   exponential, "rat.qu.s", "rat.qu.a" for rational quadratic, "matern.w",
#'   "matern.v" for Matern, "vv" for variance of Gaussian white noise. All
#'   hyperparameters should be in one list.
#'
#' @return A list containing: \describe{ \item{hyper}{Hyperparameters vector
#'   estimated from training data} \item{var.hyper}{ Variance of the estimated
#'   hyperparameters} \item{fitted.mean }{Fitted values for the training data }
#'   \item{fitted.sd }{Standard deviation of the fitted values for the training 
#'   data}
#'   \item{train.x }{ Training covariates} \item{train.y }{ Training response}
#'   \item{ train.yOri}{Original training response } \item{train.DataOri }{
#'   Original training covariates} \item{idxSubset }{Index vector identifying
#'   which observations were selected if Subset of Data was used.} \item{
#'   CovFun}{ Covariance function type} \item{ gamma}{Parameter used in powered
#'   exponential covariance function } \item{nu }{Parameter used in Matern
#'   covariance function } \item{Q}{Covariance matrix } \item{mean}{Mean
#'   function } \item{meanModel}{Mean model used} \item{meanLinearModel}{'lm'
#'   object if mean is a linear regression. NULL otherwise.} \item{conv}{An
#'   integer. 0 means converge; 1 otherwise. } \item{hyper0}{Starting point of
#'   the hyperparameters vector.} }
#'
#' @references Shi, J. Q., and Choi, T. (2011), ``Gaussian Process Regression
#'   Analysis for Functional Data'', CRC Press.
#'
#' @export
#' @examples
#' ## See examples in vignettes:
#' 
#' # vignette("gpr_ex1", package = "GPFDA")
#' # vignette("gpr_ex2", package = "GPFDA")
#' # vignette("co2", package = "GPFDA")
gpr <- function(input, response, Cov='pow.ex', m = NULL, hyper=NULL, 
                NewHyper=NULL, meanModel=0, mu=NULL, gamma=2, nu=1.5, 
                useGradient=T, iter.max=100, rel.tol=8e-10, trace=0, 
                nInitCandidates = 1000){
  
  
  if("pow.ex"%in%Cov & is.null(gamma)){
    stop("Argument 'gamma' must be informed for pow.ex kernel")
  }
  if("matern"%in%Cov & is.null(nu)){
    stop("Argument 'nu' must be informed for matern kernel")
  }
  
  if("matern"%in%Cov & !is.null(nu)){
    if("matern"%in%Cov & !(nu%in%c(3/2, 5/2)) & useGradient){
      useGradient <- F
      warning("Gradient was not used. For Matern kernel, the gradient is only 
              available if either nu=3/2 or nu=5/2. For other values of 'nu', 
              useGradient is automatically set to FALSE.")
    }}
  input <- as.matrix(input)
  
  dimData <- ncol(input)
  response <- as.matrix(response)
  n <- nrow(response)
  nrep <- ncol(response)
  
  if(!is.null(mu)){
    if(!length(mu)==n){
      stop("'mu' defined by the user must have the same length as the response 
           variable.")
    }
    mu <- matrix(rep(mu, nrep), ncol=nrep, byrow=F)
    meanModel <- 'userDefined'
  }
  
  y.original <- response
  input.original <- input
  
  idxSubset <- NULL
  if(!is.null(m)){
    if(m>n){stop("m cannot be bigger than n.")}
    idxSubset <- sort(sample(x=1:n, size=m, replace=F))
    response <- response[idxSubset,,drop=F]
    input <- input[idxSubset,,drop=F]
    if(!is.null(mu)){
      mu <- mu[idxSubset,,drop=F]
    }
  }

  if(is.null(hyper)){
    
    ## lower bounds for candidates
    hyper=list()
    if(any(Cov=='linear')){
      hyper$linear.a <- rep(log(1e-4), dimData)
      hyper$linear.i <- log(1e-4)
    }
    if(any(Cov=='pow.ex')){
      hyper$pow.ex.v <- log(1e-4)
      hyper$pow.ex.w <- rep(log(1e-4), dimData)
    }
    if(any(Cov=='matern')){
      hyper$matern.v <- log(1e-4)
      hyper$matern.w <- rep(log(1e-4), dimData)
    }
    if(any(Cov=='rat.qu')){
      hyper$rat.qu.a <- log(1e-4)
      hyper$rat.qu.v <- log(1e-4)
      hyper$rat.qu.w <- rep(log(1e-4), dimData)
    }
    hyper$vv=log(1e-4)
    hyper_low <- unlist(hyper)
    
    if(!is.null(NewHyper)){
      hyper.nam <- c(names(hyper_low), NewHyper)
      for(i in 1:length(NewHyper)){
        hyper_low <- c(hyper_low, log(1e-4))
      }
      names(hyper_low) <- hyper.nam
    }
    
    ## upper bounds for candidates
    hyper=list()
    if(any(Cov=='linear')){
      hyper$linear.a <- rep(log(1e4), dimData)
      hyper$linear.i <- log(1e4)
    }
    if(any(Cov=='pow.ex')){
      hyper$pow.ex.v <- log(1e4)
      hyper$pow.ex.w <- rep(log(1e4), dimData)
    }
    if(any(Cov=='matern')){
      hyper$matern.v <- log(1e4)
      hyper$matern.w <- rep(log(1e4), dimData)
    }
    if(any(Cov=='rat.qu')){
      hyper$rat.qu.a <- log(1e4)
      hyper$rat.qu.v <- log(1e4)
      hyper$rat.qu.w <- rep(log(1e4), dimData)
    }
    hyper$vv=log(1e4)
    hyper_upp <- unlist(hyper)
    
    if(!is.null(NewHyper)){
      hyper.nam <- c(names(hyper_upp), NewHyper)
      for(i in 1:length(NewHyper)){
        hyper_upp <- c(hyper_upp, log(1e4))
      }
      names(hyper_upp) <- hyper.nam
    }
    
    if(length(hyper_upp)!=length(hyper_low)){
      stop("hyper_upp and hyper_low must have the same dimension")
    }
    
    hyper <- hyper_low
    hyper.nam <- names(hyper_low)
    hp.name <- hyper.nam
  }
  
  if(!is.null(hyper)){
    hyper <- hyper[substr(names(hyper),1,6)%in%c(Cov,'vv')]
  }  
  hp.name <- names(unlist(hyper))
  
  if(meanModel==0) {response <- response; mu <- 0}
  if(meanModel==1) {
    mu <- mean(response)
    response <- as.matrix(response-mu)
  }
  meanLinearModel <- NULL
  if(meanModel=='t') {
    trend <- data.frame(yyy=c(response), xxx=rep(c(input), nrep))
    meanLinearModel <- lm(yyy~xxx, data=trend)
    response <- matrix(resid(meanLinearModel), nrow=nrow(response), byrow=F)
    mu <- matrix(fitted(meanLinearModel), nrow=nrow(response), byrow=F)
  }
  
  if(meanModel=='avg') {
    if(nrep<3){
      stop('Mean function can only be the average across replications when
           there are more than two replications.')
    }
    mu <- apply(response, 1, mean)
    mu <- matrix(rep(mu, nrep), ncol=nrep, byrow=F)
    response <- response - mu
  }
  
  #### Try a number of hp vector and start with the best
  candidates <- matrix(0, nInitCandidates, length(hyper_upp))
  for(iCand in 1:nInitCandidates){
    candidates[iCand,] <- runif(n=length(hyper_upp), min=hyper_low, 
                                max=hyper_upp)
  }
  colnames(candidates) <- hp.name
  cat(c('\n','--------- Initialising ---------- \n'))
  resCand <- apply(candidates, 1, function(x){
    gp.loglikelihood2(hyper.p=x, input=input, response=response,
                      Cov=Cov, gamma=gamma, nu=nu)})
  
  best_init <- candidates[which.min(resCand),]
  
  trace <- round(trace)
  if(trace>0){
    cat(c('iter:  -loglik:',hp.name,'\n'), sep='     ')
  }
  if(!useGradient){gp.Dlikelihood2 <- NULL}
  CG0 <- nlminb(start=best_init, objective=gp.loglikelihood2, 
                gradient=gp.Dlikelihood2,
                input=input, response=response, Cov=Cov, gamma=gamma, nu=nu,
                control=list(iter.max=iter.max, rel.tol=rel.tol, trace=trace))
  
  # if(trace!=F&CG0$convergence==0)
  #   cat('\n','    optimization finished. Converged.','\n')
  # if(trace!=F&CG0$convergence==1)
  #   cat('\n','    optimization finished. Failed Converge.','\n')
  if(trace>0){
    cat('\n','    optimization finished.','\n')
  }
  CG <- CG0[[1]]
  names(CG) <- hp.name
  CG.df <- data.frame(CG=CG,CG.N=substr(hp.name,1,8))
  names(CG.df) <- c('CG','CG.N')
  hyper.cg <- split(CG.df$CG,CG.df$CG.N)
  
  nkernels <- length(Cov)
  CovList <- vector('list',nkernels)
  for(i in 1:nkernels){
    CovList[i] <- list(paste0('cov.',Cov[i]))
  }
  CovL <- lapply(CovList,function(j){
    f <- get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper.cg, input=input, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper.cg, input=input, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){
      return(f(hyper=hyper.cg, input=input))
    }
  })
  if(length(CovL)==1)
    Q <- CovL[[1]]
  if(length(CovL)>1)
    Q <- Reduce('+',CovL)
  
  diag(Q) <- diag(Q)+exp(hyper.cg$vv)
  invQ <- chol2inv(chol(Q))
  
  if("matern"%in%Cov){
    if(nu%in%c(3/2, 5/2)){
      calcVarHyperPar <- T
    }else{
      calcVarHyperPar <- F
    }
  }else{
    calcVarHyperPar <- T
  }
  
  if(calcVarHyperPar){
    
    D2fxList <- vector('list',nrep)
    for(irep in 1:nrep){
      
      QR <- invQ%*%as.matrix(response[,irep])
      AlphaQ <- QR%*%t(QR)-invQ
      
      D2fx <- lapply(seq_along(hyper.cg), function(i){
        Dp <- hyper.cg[i]
        name.Dp <- names(Dp)
        f <- get(paste0('D2',name.Dp))
        if(name.Dp%in%c('pow.ex.w','pow.ex.v') )
          D2para <- f(hyper=hyper.cg, input=input, gamma=gamma, inv.Q=invQ, 
                      Alpha.Q=AlphaQ)
        if(name.Dp%in%c('matern.w','matern.v') )
          D2para <- f(hyper=hyper.cg, input=input, nu=nu, inv.Q=invQ, 
                      Alpha.Q=AlphaQ)
        if(!name.Dp%in%c('pow.ex.w','pow.ex.v','matern.w','matern.v') & 
           !name.Dp%in%c('linear.a','linear.i'))
          D2para <- f(hyper=hyper.cg, input=input, inv.Q=invQ, Alpha.Q=AlphaQ)
        if(name.Dp%in%c('linear.a'))
          D2para <- f(hyper=hyper.cg, input=input, Alpha.Q=AlphaQ)
        if(name.Dp%in%c('linear.i'))
          D2para <- f(hyper=hyper.cg, inv.Q=invQ, Alpha.Q=AlphaQ)
        return(D2para)
      })
      names(D2fx) <- names(hyper.cg)
      D2fx <- unlist(D2fx)
      D2fxList[[irep]] <- D2fx
    }
    
    D2fx <- Reduce('+', D2fxList)
    
    var_hyper <- (-1/(unlist(D2fx)*dim(input)[1]))
  }else{
    var_hyper <- NULL
  }

  fitted <- (Q-diag(exp(hyper.cg$vv),dim(Q)[1]))%*%invQ%*%(response)+mu
  fitted.var <- exp(hyper.cg$vv)*rowSums(
    (Q-diag(exp(hyper.cg$vv),dim(Q)[1]))*t(invQ))
  result <- list('hyper'=hyper.cg, 'var.hyper'=var_hyper,
                 'fitted.mean'=fitted,
                 fitted.sd=sqrt(fitted.var),
                 'train.x'=input, 'train.y'=response,
                 'train.yOri'=y.original, 
                 'train.DataOri'=input.original, 
                 'idxSubset'=idxSubset, 'CovFun'=Cov, 'gamma'=gamma, 'nu'=nu, 
                 'Q'=Q, 'inv'=invQ, 'mu'=mu, 'meanModel'=meanModel, 
                 'meanLinearModel'=meanLinearModel,
                 conv=CG0$convergence, 'hyper0'=hyper)
  class(result) <- 'gpr'
  return(result)
}





#' Prediction using Gaussian Process
#'
#' @inheritParams gpr
#' @param train A 'gpr' object obtained from 'gpr' function. Default to NULL. If
#'   NULL, learning is done based on the other given arguments; otherwise,
#'   prediction is made based on the trained model of class gpr'.
#' @param input.new Test input covariates.  It must be either a matrix, where
#'   each column represents a covariate, or a vector if there is only one
#'   covariate.
#' @param noiseFreePred Logical. If TRUE, predictions will be noise-free.
#' @param Y  Training response. It should be a matrix, where each column is a
#'   realisation. It can be a vector if there is only one realisation.
#' @param mSR Subset size m if Subset of Regressors method is used for
#'   prediction. It must be smaller than the total sample size.
#'
#' @return A list containing  \describe{ \item{pred.mean}{Mean of predictions}
#'   \item{pred.sd}{Standard deviation of predictions} \item{newdata}{Test input 
#'   data}
#'   \item{noiseFreePred}{Logical. If TRUE, predictions are noise-free.}
#'   \item{...}{Objects of 'gpr' class. } }
#' @export
#' @examples
#' ## See examples in vignettes:
#' 
#' # vignette("gpr_ex1", package = "GPFDA")
#' # vignette("gpr_ex2", package = "GPFDA")
#' # vignette("co2", package = "GPFDA")
gprPredict <- function(train=NULL, input.new=NULL, noiseFreePred=F, hyper=NULL, 
                         input=NULL, Y=NULL, mSR=NULL,
                         Cov=NULL, gamma=NULL, nu=NULL, meanModel=0, mu=0){
  input.new <- as.matrix(input.new)
  if(mu==1 & !is.null(Y)){
    mu <- mean(Y)
  }
  if(!is.null(input)) input <- as.matrix(input)
  if(!is.null(Y)) Y <- as.matrix(Y-mu)
  if("pow.ex"%in%Cov & is.null(gamma)){
    stop("Argument 'gamma' must be informed for pow.ex kernel")
  }
  if("matern"%in%Cov & is.null(nu)){
    stop("Argument 'nu' must be informed for matern kernel")
  }
  if(class(train)=='gpr'){
    hyper <- train$hyper
    input <- train$train.x
    Y <- train$train.y
    Cov <- train$CovFun
    gamma <- train$gamma
    nu <- train$nu
    mu <- train$mu
    meanModel <- train$meanModel
    meanLinearModel <- train$meanLinearModel
  }
  
  nrep <- ncol(Y)
  
  if(is.null(train)){
    train <- gpr(input=input, response=Y, Cov=Cov, hyper=hyper, 
                 gamma=gamma, nu=nu)
  }
  if(is.null(input.new)) input.new <- input

  nkernels <- length(Cov)


if(is.null(mSR)){
  
  CovList <- vector('list',nkernels)
  for(i in 1:nkernels){
    CovList[i] <- list(paste0('cov.',Cov[i]))
  }
  
  CovL1 <- lapply(CovList,function(j){
    f <- get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper, input=input, input.new=input.new, 
                                 gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper, input=input, input.new=input.new, 
                                 nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){
      return(f(hyper=hyper, input=input, input.new=input.new))}
  })
  Q1 <- Reduce('+',CovL1)
  Q1 <- t(Q1)
  
  CovL <- lapply(CovList,function(j){
    f <- get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper, input=input, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper, input=input, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){
      return(f(hyper=hyper, input=input))}
  })
  Q <- Reduce('+',CovL)
  diag(Q) <- diag(Q)+exp(hyper$vv)
  
  for(i in 1:nkernels) CovList[i] <- list(paste0('diag.',Cov[i]))
  CovLn <- lapply(CovList,function(j){
    f <- get(j)
    f(hyper=hyper, input=input.new)
  })
  Qstar <- Reduce('+',CovLn)

  invQ <- chol2inv(chol(Q))
  
  QQ1 <- invQ%*%t(Q1)
  QR <- invQ%*%Y

  if(meanModel=='t'){
    newtrend <- data.frame(xxx=input.new[,1])
    mu <- predict(meanLinearModel, newdata=newtrend)
    mu <- Q1%*%QR + matrix(rep(mu, nrep), ncol=nrep, byrow=F)
  }else{
    if(meanModel=='avg'){
      if(!all(input.new[,1]%in%input[,1])){
        stop("For predicting response values at input locations input.new, 
             the mean function can only be 'avg' across replications if
             all input.new are included in input.")
      }
      mu <- Q1%*%QR + mu
    }else{
      mu <- Q1%*%QR + mu
    }
  }

  if(noiseFreePred){
    sigma2 <- Qstar - as.matrix(diag(Q1%*%invQ%*%t(Q1)))
  }else{
    sigma2 <- Qstar - as.matrix(diag(Q1%*%invQ%*%t(Q1)))+exp(hyper$vv)
  }
  pred.sd <- sqrt(sigma2[,1])
  
}else{ # mSR set
  
  
  if(mSR>=n){stop("mSR must be smaller than n.")}
  n <- nrow(input)
  idx <- sort(sample(x=1:n, size=mSR, replace=F))
    
  CovList <- vector('list',nkernels)
  for(i in 1:nkernels){
    CovList[i]=list(paste0('cov.',Cov[i]))
  }

  Cov_m_ns <- lapply(CovList,function(j){
    f <- get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper, input=input[idx,,drop=F], 
                                 input.new=input.new, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper, input=input[idx,,drop=F], 
                                 input.new=input.new, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){
      return(f(hyper=hyper, input=input[idx,,drop=F], input.new=input.new))}
  })
  K_m_nstar <- Reduce('+',Cov_m_ns)

  Cov_mm <- lapply(CovList,function(j){
    f <- get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper, input=input[idx,,drop=F], 
                                 gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper, input=input[idx,,drop=F], 
                                 nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){
      return(f(hyper=hyper, input=input[idx,,drop=F]))}
  })
  K_mm <- Reduce('+',Cov_mm)

  Cov_m_n <- lapply(CovList,function(j){
    f <- get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper, input=input[idx,,drop=F], 
                                 input.new=input, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper, input=input[idx,,drop=F], 
                                 input.new=input, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){
      return(f(hyper=hyper, input=input[idx,,drop=F], input.new=input))}
  })
  K_m_n <- Reduce('+',Cov_m_n)

  toBeInverted <- K_m_n%*%t(K_m_n) + exp(hyper$vv)*K_mm
  diag(toBeInverted) <- diag(toBeInverted) + 1e-8
  invTerm <- chol2inv(chol(toBeInverted))
  
  centYpred <- t(K_m_nstar)%*%invTerm%*%K_m_n%*%Y
  varYpred <- exp(hyper$vv)*t(K_m_nstar)%*%invTerm%*%K_m_nstar
  range(varYpred)
  
  
  if(meanModel=='t'){
    newtrend <- data.frame(xxx=input.new[,1])
    mu <- predict(meanLinearModel,newdata=newtrend)
    mu <- centYpred + matrix(rep(mu, nrep), ncol=nrep, byrow=F)
  }else{
    if(meanModel=='avg'){
      if(!all(input.new[,1]%in%input[,1])){
        stop("For predicting response values at input locations input.new, 
             the mean function can only be 'avg' across replications if
             all input.new are included in input.")
      }
      mu <- centYpred + mu
    }else{
      mu <- centYpred + mu
    }
  }
  
  if(noiseFreePred){
    sigma2 <- varYpred
  }else{
    sigma2 <- varYpred + exp(hyper$vv)
  }
  
  sigma2
  if(all(sigma2>0)){
    pred.sd <- sqrt(sigma2[,1])
  }else{
    pred.sd <- NA
    warning("Subset of Regressors give negative predictive variance.")
  }
  
}

  result <- c(list('pred.mean'=mu,
                'pred.sd'=pred.sd,
                'newdata'=input.new), 
              'noiseFreePred'=noiseFreePred,
              unclass(train))
  class(result) <- 'gpr'
  
  return(result)
}






################################# tools #####################################
rmse <- function(t, a){ 
  y <- sqrt(sum((a-t)^2)/length(t))
  return(y)
}

########################### derivatives #####################################

Dloglik.linear.a <- function(hyper, input, AlphaQ){
  Dlinear.aj <- sapply(1:ncol(input),function(i){
    Xi <- as.matrix(input[,i])
    A1 <- as.matrix(1)
    DPsi <- exp(hyper$linear.a[i])*distMatLinearSq(input=Xi, A=A1)
    res <- sum(diag(AlphaQ%*%DPsi))
    return(res)
  })
  return(Dlinear.aj)
}

Dloglik.linear.i <- function(hyper, Alpha, invQ){
  constMat <- matrix(exp(hyper$linear.i), nrow=nrow(invQ), ncol=ncol(invQ))
  Dfa0 <- 0.5*sum(diag(  (Alpha%*%t(Alpha) - invQ)%*%constMat ))
  return(Dfa0)
}

Dloglik.pow.ex.w <- function(hyper, input, AlphaQ, gamma){
  Dpow.ex.wj=sapply(1:ncol(input),function(i){
    Xi <- as.matrix(input[,i])
    A1 <- as.matrix(1)
    DPsi <- -cov.pow.ex(hyper=hyper, input=input, gamma=gamma)*
      distMatSq(input=Xi, A=A1, power=gamma)*exp(hyper$pow.ex.w[i])
    res <- 0.5*sum(diag(AlphaQ%*%DPsi))
    return(res)
  })
  return(Dpow.ex.wj)
}

Dloglik.pow.ex.v <- function(hyper, input, AlphaQ, gamma){
  DDpow.ex.v <- cov.pow.ex(hyper,input, gamma=gamma)
  Dpow.ex.v <- 0.5*sum(diag(AlphaQ%*%DDpow.ex.v))
  return(Dpow.ex.v)
}

#######################

Dloglik.matern.v <- function(hyper, input, AlphaQ, nu){
  DDmatern.v <- cov.matern(hyper=hyper, input=input, input.new=NULL, nu=nu)
  Dmatern.v <- 0.5*sum(diag(AlphaQ%*%DDmatern.v))
  return(Dmatern.v)
}

Dloglik.matern.w <- function(hyper, input, AlphaQ, nu){
  
  input <- as.matrix(input)
  cc <- exp(hyper$matern.v)
  dimData <- ncol(input)
  if(dimData==1){
    A <- as.matrix(exp(hyper$matern.w))
  }else{
    A <- diag((exp(hyper$matern.w)))
  }
  
  dist_C <- distMatSq(input=input, A=A, power=2)
  dist_D <- sqrt(dist_C)
  
  Dmatern.wj <- sapply(1:ncol(input),function(i){
    Xi <- as.matrix(input[,i])
    A1 <- as.matrix(1)
    dist_i <- distMatSq(input=Xi, A=A1, power=2)
    if(nu==3/2){
      DPsi <- -1.5*cc*exp(hyper$matern.w[i])*dist_i*exp(-sqrt(3)*dist_D)
    }
    if(nu==5/2){
      DPsi <- cc*exp(hyper$matern.w[i])*dist_i*exp(-sqrt(5)*dist_D)*(5/6)*
        (-1-sqrt(5)*dist_D)
    }
    res <- 0.5*sum(diag(AlphaQ%*%DPsi))
    return(res)
  })
  return(Dmatern.wj)
}


Dloglik.rat.qu.w <- function(hyper, input, AlphaQ){

  dimData <- ncol(input)
  d1list <- vector("list", dimData)
  for(i in 1:dimData){
    covmatrix <- cov.rat.qu(hyper=hyper, input=input)
    
    Xi <- as.matrix(input[,i])
    A1 <- as.matrix(1)
    hyperAdj <- hyper
    hyperAdj$rat.qu.v <- 0
    covmatrixAdj_a <- cov.rat.qu(hyper=hyperAdj, input=input)
    hyperAdj$rat.qu.a <- 0
    covmatrixAdj_0 <- cov.rat.qu(hyper=hyperAdj, input=input)
    
    d1 <- -exp(hyper$rat.qu.a)*exp(hyper$rat.qu.v)*covmatrixAdj_a*
      covmatrixAdj_0*distMatSq(input=Xi, A=A1, power=2)*exp(hyper$rat.qu.w[i])
    d1list[[i]] <- d1
  }

  Drat.qu.wj <- sapply(1:dimData, function(i){sum(0.5*AlphaQ*d1list[[i]])} )
  return(Drat.qu.wj)
}

Dloglik.rat.qu.a <- function(hyper, input, AlphaQ){
  
  covmatrix <- cov.rat.qu(hyper,input)
  
  hyper <- lapply(hyper,exp)
  input <- as.matrix(input)
  Q <- ncol(input)
  if(Q==1){
    A <- as.matrix(hyper$rat.qu.w)
  }else{
    A <- diag(hyper$rat.qu.w)
  }
  v.power <- distMatSq(input=input,A=A,power=2)
  log_term <- log( 1 + v.power )
  
  DDrat.qu.a <- log_term*covmatrix*(-hyper$rat.qu.a)
  Drat.qu.a <- 0.5*sum(diag(AlphaQ%*%DDrat.qu.a))

  return(Drat.qu.a)
}


Dloglik.rat.qu.v <- function(hyper, input, AlphaQ){
  DDrat.qu.v <- cov.rat.qu(hyper,input)
  Drat.qu.v <- 0.5*sum(diag(AlphaQ%*%DDrat.qu.v))
  return(Drat.qu.v)
}



Dloglik.vv <- function(hyper, Alpha, invQ){
  Dfvv <- 0.5*sum(diag(Alpha%*%t(Alpha) - invQ))*exp(hyper$vv)
  return(Dfvv)
}


D2linear.a <- function(hyper, input, Alpha.Q){
  D2linear.aj <- sapply(1:ncol(input),function(i){
    Xi <- as.matrix(input[,i])
    A1 <- as.matrix(1)
    DPsi <- exp(hyper$linear.a[i])*distMatLinearSq(input=Xi, A=A1)
    res <- sum(diag(Alpha.Q%*%DPsi))
    return(res)
  })
  return(D2linear.aj)
}


D2linear.i <- function(hyper, inv.Q, Alpha.Q){
  constMat <- matrix(exp(hyper$linear.i), nrow=nrow(inv.Q), ncol=ncol(inv.Q))
  D2flinear.i <- D2(constMat, constMat,inv.Q,Alpha.Q)
  return(D2flinear.i)
}



D2pow.ex.w <- function(hyper,input,gamma,inv.Q,Alpha.Q){
  
  dimData <- ncol(input)
  d1list <- vector("list", dimData)
  d2list <- vector("list", dimData)
  for(i in 1:dimData){
    covmatrix <- cov.pow.ex(hyper=hyper, input=input, gamma=gamma)
    Xi <- as.matrix(input[,i])
    A1 <- as.matrix(1)
    d1 <- -covmatrix*distMatSq(input=Xi, A=A1, power=gamma)*exp(hyper$pow.ex.w[i])
    d1list[[i]] <- d1
    d2 <- covmatrix*
      (distMatSq(input=Xi, A=A1, power=2*gamma)*exp(2*hyper$pow.ex.w[i]) - 
         distMatSq(input=Xi, A=A1, power=gamma)*exp(hyper$pow.ex.w[i]))
    d2list[[i]] <- d2
  }

  D2pow.ex.wj <- sapply(1:dimData, function(i){
    D2(d1=d1list[[i]], 
       d2=d2list[[i]], 
       inv.Q=inv.Q, Alpha.Q=Alpha.Q)} )

  return(D2pow.ex.wj)
}


D2pow.ex.v <- function(hyper, input, gamma, inv.Q, Alpha.Q){
  DDpow.ex.v <- cov.pow.ex(hyper=hyper, input=input, input.new=NULL, 
                           gamma=gamma)
  D2pow.ex.v <- D2(DDpow.ex.v, DDpow.ex.v, inv.Q, Alpha.Q)  
  return(D2pow.ex.v)
}

D2matern.v <- function(hyper, input, nu, inv.Q, Alpha.Q){
  DDmatern.v <- cov.matern(hyper=hyper, input=input, input.new=NULL, nu=nu)
  D2matern.v <- D2(DDmatern.v, DDmatern.v, inv.Q, Alpha.Q)  
  return(D2matern.v)
}



D2matern.w <- function(hyper, input, nu, inv.Q, Alpha.Q){
  
  input <- as.matrix(input)
  cc <- exp(hyper$matern.v)
  dimData <- ncol(input)
  if(dimData==1){
    A <- as.matrix(exp(hyper$matern.w))
  }else{
    A <- diag((exp(hyper$matern.w)))
  }

  dist_C <- distMatSq(input=input, A=A, power=2)
  dist_D <- sqrt(dist_C)

  d1list <- vector("list", dimData)
  d2list <- vector("list", dimData)
  
  for(i in 1:dimData){
    Xi <- as.matrix(input[,i])
    A1 <- as.matrix(1)
    dist_i <- distMatSq(input=Xi, A=A1, power=2)
    
    if(nu==3/2){
      d1list[[i]] <- -1.5*cc*exp(hyper$matern.w[i])*dist_i*exp(-sqrt(3)*dist_D)
      d2list[[i]] <- d1list[[i]]*(1 - 0.5*sqrt(3)*dist_C^(-0.5))*
        exp(hyper$matern.w[i])*dist_i
      diag(d2list[[i]]) <- 0
    }
    
    if(nu==5/2){
      d1list[[i]] <- cc*exp(hyper$matern.w[i])*dist_i*exp(-sqrt(5)*dist_D)*
        (5/6)*(-1-sqrt(5)*dist_D)
      d2list[[i]] <- (-5/6)*cc*dist_i*exp(hyper$matern.w[i])*
        exp(-sqrt(5)*dist_D)*(
        1+sqrt(5)*dist_D - 2.5*exp(hyper$matern.w[i])*dist_i)
    }
    
  }
  D2matern.wj <- sapply(1:dimData, function(i){
    D2(d1=d1list[[i]], 
       d2=d2list[[i]], 
       inv.Q=inv.Q, Alpha.Q=Alpha.Q)} )
  
  return(D2matern.wj)
  
}


######################################################################

# covmatrixAdj_a will use unit signal variance v and usual alpha parameter
# covmatrixAdj_0 will use unit signal variance v and alpha=1 (ie, log(alpha)=0, 
# exp(log(1)) = exp(0))
D2rat.qu.w <- function(hyper, input, inv.Q, Alpha.Q){
  
  dimData <- ncol(input)
  d1list <- vector("list", dimData)
  d2list <- vector("list", dimData)
  for(i in 1:dimData){
    covmatrix <- cov.rat.qu(hyper=hyper, input=input)
    
    Xi <- as.matrix(input[,i])
    A1 <- as.matrix(1)

    hyperAdj <- hyper
    hyperAdj$rat.qu.v <- 0
    covmatrixAdj_a <- cov.rat.qu(hyper=hyperAdj, input=input)
    hyperAdj$rat.qu.a <- 0
    covmatrixAdj_0 <- cov.rat.qu(hyper=hyperAdj, input=input)
    
    d1 <- -exp(hyper$rat.qu.a)*exp(hyper$rat.qu.v)*covmatrixAdj_a*
      covmatrixAdj_0*distMatSq(input=Xi, A=A1, power=2)*exp(hyper$rat.qu.w[i])
    d1list[[i]] <- d1
    
    # calc D2rat.qu
    d2 <- -exp(hyper$rat.qu.a)*exp(hyper$rat.qu.v)*covmatrixAdj_a*
      covmatrixAdj_0*((-exp(hyper$rat.qu.a)-1)*covmatrixAdj_0*
                        distMatSq(input=Xi, A=A1, power=4)*
                        exp(2*hyper$rat.qu.w[i]) + 
        distMatSq(input=Xi, A=A1, power=2)*exp(hyper$rat.qu.w[i])
    )
    
    d2list[[i]] <- d2
  }
  
  D2rat.qu.wj <- sapply(1:dimData, function(i){
    D2(d1=d1list[[i]], 
       d2=d2list[[i]], 
       inv.Q=inv.Q, Alpha.Q=Alpha.Q)} )
  
  return(D2rat.qu.wj)
}


D2rat.qu.a <- function(hyper, input, inv.Q, Alpha.Q){
  
  covmatrix <- cov.rat.qu(hyper, input)

  hyper <- lapply(hyper, exp)
  input <- as.matrix(input)
  Q <- ncol(input)
  if(Q==1){
    A <- as.matrix(hyper$rat.qu.w)
  }else{
    A <- diag(hyper$rat.qu.w)
  }
  v.power <- distMatSq(input=input, A=A, power=2)
  log_term <- log( 1 + v.power )
  
  d1 <- log_term*covmatrix*(-hyper$rat.qu.a)
  d2 <- (d1 + covmatrix)*(-hyper$rat.qu.a)*log_term

  D2rat.qu.a <- D2(d1,d2,inv.Q,Alpha.Q)
  
  return(D2rat.qu.a)
}

D2rat.qu.v <- function(hyper, input, inv.Q, Alpha.Q){
  covmatrix <- cov.rat.qu(hyper, input)
  D2rat.qu.v <- D2(covmatrix, covmatrix, inv.Q, Alpha.Q)  
  return(D2rat.qu.v)
}


D2vv <- function(hyper, input, inv.Q, Alpha.Q){
  D2fvv <- D2(diag(exp(hyper$vv),dim(input)[1]), 
              diag(exp(hyper$vv),dim(input)[1]),inv.Q,Alpha.Q)
  return(D2fvv)
}


#' Second derivative of the likelihood
#'
#' Calculate the second derivative of the likelihood function with respect to
#' one of the hyperparameters, given the first and second derivative of the
#' kernel with respect to that hyperparameter.
#'
#' @param d1 First derivative of the kernel function with respect to the required
#'  hyperparameter.
#' @param d2 Second derivative of the kernel function with respect to the
#'  required hyperparameter.
#' @param inv.Q Inverse of covariance matrix Q.
#' @param Alpha.Q  This is alpha * alpha'- invQ, where invQ is the inverse of the
#'  covariance matrix Q, and alpha = invQ * Y, where Y is the response.
#'
#' @details The function calculates the second derivative of the log-likelihood,
#'  using the first and second derivative of the kernel functions.
#' @return A number.
#'
#' @references Shi, J. Q., and Choi, T. (2011), ``Gaussian Process Regression
#'  Analysis for Functional Data'', CRC Press.
#' @examples
#' ## This function is used in the vignette 'co2':
#' # vignette("co2", package = "GPFDA")
#' @export
D2 <- function(d1, d2, inv.Q, Alpha.Q){
  Aii <- t(d1)%*%inv.Q%*%d1
  al <- Alpha.Q+inv.Q
  return(0.5*(sum(Alpha.Q*(d2-Aii))-sum(al*Aii)))
}

diag.linear <- function(hyper, input){
  Qstar <- exp(hyper$linear.i) + input^2%*%matrix(exp(hyper$linear.a))
  return(Qstar)
}

diag.pow.ex <- function(hyper, input){
  Qstar <- rep(exp(hyper$pow.ex.v),dim(input)[1])
  return(Qstar)
}

diag.matern <- function(hyper, input){
  Qstar <- rep(exp(hyper$matern.v),dim(input)[1])
  return(Qstar)
}

diag.rat.qu <- function(hyper, input){
  Qstar <- rep(exp(hyper$rat.qu.v),dim(input)[1])
  return(Qstar)
}


#' Plot Gaussian Process regression -- training and prediction
#'
#' Plot Gaussian Process for a given an object of class 'gpr'.
#'
#' @param x The 'gpr' object from either training or predicting of the Gaussian
#'   Process.
#' @param fitted Logical. Plot fitted values or not. Default to FALSE. If FALSE,
#'   plot the predictions.
#' @param col.no Column number of the input matrix. If the input matrix has more
#'   than one columns, than one of them will be used in the plot. Default to be
#'   the first one.
#' @param ylim Range value for y-axis.
#' @param cex.points Graphical parameter
#' @param lwd.points Graphical parameter
#' @param pch Graphical parameter
#' @param lwd Graphical parameter
#' @param realisation Integer identifying which realisation should be plotted
#'   (if there are multiple).
#' @param main Title for the plot
#' @param ... Graphical parameters passed to plot().
#' @importFrom  graphics polygon
#' @importFrom  graphics points
#' @importFrom  graphics matpoints
#' @importFrom  graphics matlines
#' @importFrom  graphics lines
#' @importFrom  grDevices rgb
#' @return A plot
#' @export
#' @examples
#' ## See examples in vignette:
#' # vignette("gpr_ex1", package = "GPFDA")
plot.gpr <- function(x, fitted=F, col.no=1, ylim=NULL, realisation=NULL, main=NULL, 
                     cex.points=NULL, lwd.points=NULL, pch=NULL, lwd=NULL, ...){
  obj <- x
  if(fitted==T){
    if(is.null(obj$fitted.mean)){
      warning('fitted values not found, ploting predicted values')
      type <- 'Prediction'
      mu <- obj$pred.mean
      sd <- obj$pred.sd
      x <- obj$newdata
      X <- obj$train.x
      Y <- obj$train.yOri
    }
    if(!is.null(obj$fitted.mean)){
      type <- 'Fitted values'
      mu <- obj$fitted.mean
      sd <- obj$fitted.sd
      X <- obj$train.x
      Y <- obj$train.yOri
      x <- X
    }
  }else{
    if(is.null(obj$pred.mean)){
      warning('predicted values not found, ploting fitted values')
      type <- 'Fitted values'
      mu <- obj$fitted.mean
      sd <- obj$fitted.sd
      X <- obj$train.x
      Y <- obj$train.yOri
      x <- X
    }
    if(!is.null(obj$pred.mean)){
      type <- 'Prediction'
      mu <- obj$pred.mean
      sd <- obj$pred.sd
      x <- obj$newdata
      X <- obj$train.x
      Y <- obj$train.yOri
    }
  }
  if(dim(X)[1]<=150|length(X)<=150){
    pchType <- 4
    PcexNo <- 1
    LcexNo <- 1.5
    PLcexNo <- 2
  }else{
    pchType <- 20
    PcexNo <- 0.1
    LcexNo <- 0.8
    PLcexNo <- 1
  }
  
  if(!is.null(cex.points)){
    PcexNo <- cex.points
  }
  if(!is.null(lwd.points)){
    PLcexNo <- lwd.points
  }
  if(!is.null(pch)){
    pchType <- pch
  }
  if(!is.null(lwd)){
    LcexNo <- lwd
  }
  
  
  
  
  noiseFreePred <- obj$noiseFreePred
  if(type=='Prediction'){
    if(noiseFreePred){
      type <- "Noise-free prediction"
    }
  }
  
  if(!is.null(main)){
    type <- main
  }
  if(!is.null(realisation)){
    mu <- mu[,realisation]
    Y <- Y[,realisation]
  }
  upper <- mu+1.96*(sd);
  lower <- mu-1.96*(sd);
  if(is.null(ylim)){
    ylim <- range(upper,lower,Y)
  }
  plot(-100,-100,col=0,xlim=range(X[,col.no],x[,col.no]), ylim=ylim, main=type, 
       xlab="input ", ylab="response",...)
  polygon(c(x[,col.no], rev(x[,col.no])), c(upper, rev(lower)), 
          col=rgb(127,127,127,120, maxColorValue=255), border=NA)
  points(X[,col.no], Y, pch=pchType, col=2, cex=PcexNo, lwd=PLcexNo)
  lines(x[,col.no], mu, col=4, lwd=LcexNo)
}


#' Draw an image plot for a given two-dimensional input
#'
#' @param response Data to be plotted (e.g. matrix of predictions)
#' @param input Matrix of two columns representing the input coordinates.
#' @param realisation Integer identifying which realisation should be plotted
#'   (if there are multiple).
#' @param n1 Number of datapoints in the first coordinate direction
#' @param n2 Number of datapoints in the second coordinate direction
#' @param main Title for the plot
#' @param zlim Range of z-axis
#' @param cex.axis Graphical parameter
#' @param cex.lab Graphical parameter
#' @param font.main Graphical parameter
#' @param cex.main Graphical parameter
#' @param legend.cex.axis Graphical parameter
#' @param legend.width Graphical parameter
#' @param mar Graphical parameter
#' @param oma Graphical parameter
#' @param nGrid Dimension of output grid in each coordinate direction
#' @param enlarge_zlim Additional quantity to increase the range of zlim
#' @importFrom  graphics par
#' @importFrom  graphics image
#' @importFrom  graphics title
#' @importFrom  interp interp
#' @importFrom  fields image.plot
#' @importFrom  fields tim.colors
#' @return A plot
#' @export
#' @examples
#' ## See examples in vignette:
#' # vignette("gpr_ex2", package = "GPFDA")
plotImage <- function(response, input, realisation=1, 
                        n1, n2, 
                        main=" ", zlim=NULL,
                        cex.axis=1, cex.lab=2.5, 
                        legend.cex.axis=1, font.main=2, cex.main=2, legend.width=2,
                        mar=c(2.1,2.1,3.1,6.1), 
                        oma=c( 0,1,0,0),
                        nGrid=200,
                        enlarge_zlim=NULL){
  
  if(!is.matrix(input)){
    stop("The argument 'input' must be a matrix")
  }
  if(ncol(input)!=2){
    stop("The argument 'input' must be a matrix of 2 columns")
  }
  
  n <- nrow(input)
  response <- as.matrix(response)
  if(nrow(response)!=n){
    stop("The arguments 'response' and 'input' must have the same sample size.")
  }
  
  opar <- par(no.readonly = TRUE)
  par(mar=mar, oma=oma)
  
  if(is.null(enlarge_zlim)){
    enlarge_zlim <- 0.2*c(-1,1)
  }
  if(is.null(zlim)){
    zlim <- range(response) + enlarge_zlim
  }
  responseMat <- matrix(response[,realisation], nrow=n1, ncol=n2, byrow=F)
  
  akima_li <- interp(x=input[,1], y=input[,2], z=as.numeric(responseMat),
                     nx=nGrid, ny=nGrid, linear = T)
  
  image(x = akima_li$x, y=akima_li$y, z=akima_li$z, col=tim.colors(),
        xlab=NA, ylab=NA, cex.lab=cex.lab, cex.axis=cex.axis, zlim=zlim)
  
  title(main = main, font.main=font.main, cex.main=cex.main, line=1)
  
  image.plot( legend.only=TRUE, zlim=zlim, 
              legend.cex=0.5, legend.width=legend.width,
                      axis.args=list(cex.axis=legend.cex.axis))
  
  par(opar)
}

# ########################### likelihood ######################################


gp.loglikelihood2 <- function(hyper.p,input, response,Cov,gamma,nu){
  
  input <- as.matrix(input)
  datadim <- dim(input)
  
  hp.class <- substr(names(hyper.p),1,8)
  kernel.class <- unique(substr(names(hyper.p),1,6))
  hp.class <- data.frame(class=hp.class,hp=hyper.p)
  names(hp.class) <- c('class','hp')
  hp.list <- split(hp.class$hp,hp.class$class)
  hyper.p <- hp.list
  
  nkernels <- length(Cov)
  CovList <- vector('list',nkernels)
  for(i in 1:nkernels){
    CovList[i]=list(paste0('cov.',Cov[i]))
  }
  CovL <- lapply(CovList,function(j){
    f <- get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper.p, input=input, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper.p, input=input, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){
      return(f(hyper=hyper.p, input=input))}
  })
  Q <- Reduce('+',CovL)
  diag(Q) <- diag(Q)+exp(hyper.p$vv) + 1e-8
 
  n <- nrow(response)
  nrep <- ncol(response)
  
  G <- chol(Q)
  logdetQ <- 2*sum(log(diag(G)))

  tresp_invQ_resp <- t(response)%*%chol2inv(G)%*%response
  if(nrep==1){
    fX <- 0.5*logdetQ + 0.5*tresp_invQ_resp + 0.5*n*log(2*pi)
  }else{
    fX <- nrep*0.5*logdetQ + 0.5*sum(diag( tresp_invQ_resp )) + 
      nrep*0.5*n*log(2*pi)
  }
  fX <- as.numeric(fX)

  return(fX)
}

# ########################### gradient ######################################

gp.Dlikelihood2 <- function(hyper.p,  input, response, Cov, gamma, nu){
  
  input <- as.matrix(input)
  datadim <- dim(input)
  
  hp.class <- substr(names(hyper.p),1,8)
  kernel.class <- unique(substr(names(hyper.p),1,6))
  hp.class <- data.frame(class=hp.class,hp=hyper.p)
  names(hp.class) <- c('class','hp')
  hyper.p <- split(hp.class$hp,hp.class$class)
  
  nkernels <- length(Cov)
  CovList <- vector('list',nkernels)
  for(i in 1:nkernels){
    CovList[i]=list(paste0('cov.',Cov[i]))
  }
  CovL <- lapply(CovList,function(j){
    f <- get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper.p, input=input, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper.p, input=input, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){
      return(f(hyper=hyper.p, input=input))
    }
  })
  Q <- Reduce('+',CovL)
  diag(Q) <- diag(Q)+exp(hyper.p$vv) + 1e-8

  invQ <- chol2inv(chol(Q))
  
  nrep <- ncol(response)

DfxList <- vector('list',nrep)
for(irep in 1:nrep){
  Alpha <- invQ%*%as.matrix(response[,irep])
  AlphaQ <- Alpha%*%t(Alpha)-invQ
  
  Dfx <- lapply(seq_along(hyper.p),function(i){
    Dp <- hyper.p[i];
    name.Dp <- names(Dp)
    f <- get(paste0('Dloglik.',name.Dp))
    if(name.Dp%in%c('pow.ex.w','pow.ex.v') )
      Dpara <- f(hyper=hyper.p, input=input, AlphaQ=AlphaQ, gamma=gamma)
    if(name.Dp%in%c('matern.w','matern.v') )
      Dpara <- f(hyper=hyper.p, input=input, AlphaQ=AlphaQ, nu=nu)
    if(name.Dp%in%c('rat.qu.w','rat.qu.v','rat.qu.a') )
      Dpara <-f(hyper=hyper.p, input=input, AlphaQ=AlphaQ)
    if(name.Dp%in%c('linear.a') )
      Dpara <- f(hyper=hyper.p, input=input, AlphaQ=AlphaQ)
    if(name.Dp%in%c('linear.i') )
      Dpara <- f(hyper=hyper.p, Alpha=Alpha, invQ=invQ)
    if(name.Dp=='vv')
      Dpara <- f(hyper=hyper.p, Alpha=Alpha, invQ=invQ)
    if(substr(name.Dp, 1, 6)=='custom')
      Dpara <- f(hyper=hyper.p, input=input, AlphaQ=AlphaQ)
    return(Dpara)
  })
  
  names(Dfx) <- names(hyper.p)
  Dfx <- - unlist(Dfx) # Minus loglik
  DfxList[[irep]] <- Dfx
}
  
Dfx <- Reduce('+', DfxList)
  
  return(Dfx)
}



######################## Covariance functions ############################ 

#' Stationary powered exponential covariance function.
#'
#' @param hyper The hyperparameters. Must be a list with certain names.
#' @param Data The input data. Must be a vector or a matrix.
#' @param Data.new The data for prediction. Must be a vector or a matrix.
#'   Default to NULL.
#' @param gamma Power parameter that cannot be estimated by simple non-linear
#'   optimization.
#'
#' @details The names for the hyper parameters should be:"linear.a" for linear
#'   covariance function, "pow.ex.w", "pow.ex.v" for power exponential,
#'   "rat.qu.s", "rat.qu.a" for rational quadratic, "vv" for white noise. All
#'   hyper parameters should be in one list.
#' @references Shi, J. Q., and Choi, T. (2011), ``Gaussian Process Regression
#'   Analysis for Functional Data'', CRC Press.
#' @return Covariance matrix
#' @export
#'
cov.pow.ex <- function(hyper,Data,Data.new=NULL,gamma=2){
  
  hyper <- lapply(hyper,exp)
  Data <- as.matrix(Data)
  Q <- ncol(Data)
  if(Q==1){
    A <- as.matrix(hyper$pow.ex.w)
  }else{
    A <- diag(hyper$pow.ex.w)
  }
  
  if (is.null(Data.new)) {
    distMat <- DistMat_sq(X = Data, A=A, power=gamma)
    exp.v.power <- hyper$pow.ex.v*exp(-distMat)
  }else{
    Data.new <- as.matrix(Data.new)
    distMat <- DistMat(X = Data, Xnew = Data.new,A=A, power=gamma)
    exp.v.power <- hyper$pow.ex.v*exp(-distMat)
  }

  return(exp.v.power)
}


#' Stationary rational quadratic covariance function
#'
#' @inheritParams cov.pow.ex
#'
#' @return Covariance matrix
#' @export
#'
cov.rat.qu <- function(hyper,Data,Data.new=NULL){
  
  hyper <- lapply(hyper,exp)
  Data <- as.matrix(Data)
  Q <- ncol(Data)
  if(Q==1){
    A <- as.matrix(hyper$rat.qu.w)
  }else{
    A <- diag(hyper$rat.qu.w)
  }
  
  if (is.null(Data.new)) {
    distMat <- DistMat_sq(X=Data, A=A, power=2)
    covratqu <- hyper$rat.qu.v*(1 + distMat)^(-hyper$rat.qu.a)
  }else{
    Data.new <- as.matrix(Data.new)
    distMat <- DistMat(X=Data, Xnew=Data.new, A=A, power=2)
    covratqu <- hyper$rat.qu.v*(1 + distMat)^(-hyper$rat.qu.a)
  }
  
  return(covratqu)
}

#' Stationary Matern covariance function
#'
#' @inheritParams cov.pow.ex
#' @param nu Smoothness parameter of the Matern class. Must be a positive value.
#' @return Covariance matrix
#' @export
#'
cov.matern <- function(hyper,Data,Data.new=NULL, nu){
  
  Data <- as.matrix(Data)
  cc <- exp(hyper$matern.v)
  Q <- ncol(Data)
  if(Q==1){
    A <- as.matrix(exp(hyper$matern.w))
  }else{
    A <- diag((exp(hyper$matern.w)))
  }
  
  if (is.null(Data.new)) {
    
    if(nu==3/2){
      distmat <- sqrt(DistMat_sq(X=Data, A=A, power=2))
      covmatern <- cc*(1+sqrt(3)*distmat)*exp(-sqrt(3)*distmat)
    }else{
      if(nu==5/2){
        distmat <- sqrt(DistMat_sq(X=Data, A=A, power=2))
        covmatern <- cc*(1+sqrt(5)*distmat+(5/3)*distmat^2)*exp(-sqrt(5)*distmat)
      }else{
        covmatern <- CovMaternCpp_sq(X=Data, cc=cc, A=A, nu=nu)
      }
    }

  }else{
    Data.new <- as.matrix(Data.new)
    if(nu==3/2){
      distmat <- sqrt(DistMat(X=Data, Xnew=Data.new, A=A, power=2))
      covmatern <- cc*(1+sqrt(3)*distmat)*exp(-sqrt(3)*distmat)
    }else{
      if(nu==5/2){
        distmat <- sqrt(DistMat(X=Data, Xnew=Data.new, A=A, power=2))
        covmatern <- cc*(1+sqrt(5)*distmat+(5/3)*distmat^2)*exp(-sqrt(5)*distmat)
      }else{
        covmatern <- CovMaternCpp(X=Data, Xnew=Data.new, cc=cc, A=A, nu=nu)   
      }
    }
  }
  
  return(covmatern)
}

#' Linear covariance function
#'
#' @inheritParams cov.pow.ex
#'
#' @return Covariance matrix
#' @export
#'
cov.linear <- function(hyper,Data,Data.new=NULL){
  
  hyper <- lapply(hyper,exp)
  Data <- as.matrix(Data)
  Q <- ncol(Data)
  if(Q==1){
    A <- as.matrix(hyper$linear.a)
  }else{
    A <- diag(hyper$linear.a)
  }
  
  if (is.null(Data.new)) {
    cov.lin <- DistMatLinear_sq(X=Data, A=A)
  }else{
    Data.new <- as.matrix(Data.new)
    cov.lin <- DistMatLinear(X=Data, Xnew=Data.new, A=A)
  }
  
  cov.lin <- hyper$linear.i + cov.lin
  
  return(cov.lin)
}


#' Gaussian Process regression
#'
#' Gaussian Process regression for a single or multiple independent
#' realisations.
#'
#' @param Data The input data from train data. Matrix or vectors are both
#'   acceptable. Some data.frames are not acceptable.
#' @param response The response data from train data. Matrix or vectors are both
#'   acceptable. Some data.frames are not acceptable.
#' @param Cov Covariance function(s) to use. Default to 'power.ex'.
#' @param m If Subset of Data is to be used, m denotes the subset size and
#'   cannot be larger than the total sample size. Default set to NULL.
#' @param hyper The hyperparameters. Default to NULL. If not NULL, then must be
#'   a list with certain names.
#' @param NewHyper Vector of the names of the new hyper parameters from
#'   customized kernel function. The names of the hyper-parameters must have the
#'   format: xxxxxx.x, i.e. '6 digit' plus 'a dot' plus '1 digit'. This is
#'   required for both 'hyper' and 'NewHyper'
#' @param meanModel Type of mean.
#' @param mean Is the mean taken out when analysis? Default to be 0, which
#'   assumes the mean is zero. if assume mean is a constant, mean=1; if assume
#'   mean is a linear trend, mean='t'.
#' @param gamma Power parameter used in powered exponential kernel function.
#' It must be 0<gamma<=2.
#' @param nu Smoothness parameter of the Matern class. It must be a positive value.
#' @param useGradient Logical. If TRUE, first derivatives will be used in the
#'   optimization.
#' @param itermax Number of maximum iteration in optimization function. Default
#'   to be 100. Normally the number of optimization steps is around 20. If
#'   reduce 'reltol', the iterations needed will be less.
#' @param reltol Relative convergence tolerance. Smaller reltol means more
#'   accurate and slower to converge.
#' @param trace The value of the objective function and the parameters is
#'   printed every trace'th iteration. Defaults to 0 which indicates no trace
#'   information is to be printed.
#' @param nInitCandidates Number of initial hyperparameter vectors. The
#'   optimization starts with the best.
#'
#' @details The most important function in the package, for fitting the GP model
#'   and store everything necessary for prediction. The optimization used in the
#'   function is 'nlminb'. Optimization might break down if the noise for the
#'   curve are too far away from normal. Jitter, LU decomposition and sparse
#'   matrix inverse are used to ensure the matrix inverse can always get an
#'   answer. The names for the hyper parameters should be:"linear.a" for linear
#'   covariance function, "pow.ex.w", "pow.ex.v" for power exponential,
#'   "rat.qu.s", "rat.qu.a" for rational quadratic, "matern" for Matern, "vv"
#'   for Gaussian white noise. All hyper parameters should be in one list.
#'
#' @return A list containing: \describe{ 
#' \item{hyper}{Hyper-parameter estimated from training data}
#' \item{var.hyper}{ Variance of the estimated hyper-parameters} 
#' \item{fitted.mean }{Fitted value of training data } 
#' \item{fitted.sd }{Standard deviation of the fitted value of training data} 
#' \item{train.x }{ Training covariates} 
#' \item{train.y }{ Training response} 
#' \item{ train.yOri}{Original training response } 
#' \item{train.DataOri }{ } 
#' \item{idxSubset }{Index vector identifying which observations were selected 
#' if Subset of Data was used.} 
#' \item{ CovFun}{ Covariance function type} 
#' \item{ gamma}{Parameter used in powered exponential covariance function } 
#' \item{nu }{Parameter used in Matern covariance function } 
#' \item{Q}{Covariance matrix } 
#' \item{mean}{Mean function } 
#' \item{meanModel}{ CHECK: 'lm' object if mean is a linear regression. NULL otherwise.} 
#' \item{meanLinearModel}{ } 
#' \item{conv}{0 means converge; 1 otherwise. } 
#' \item{hyper0}{ starting point of the hyper-parameters}
#'  }
#'
#' @references Shi, J. Q., and Choi, T. (2011), ``Gaussian Process Regression
#'   Analysis for Functional Data'', CRC Press.
#'
#' @export
#' 
gpr <- function(Data, response, Cov='pow.ex', 
                m = NULL, hyper=NULL, NewHyper=NULL, meanModel=0, mean=NULL, 
                gamma=NULL, nu=NULL,
                useGradient=T, itermax=100, reltol=8e-10, trace=0,
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
      warning("Gradient was not used.
      For Matern kernel, the gradient is only available if either nu=3/2 or nu=5/2.
      For other values of 'nu', useGradient is automatically set to FALSE.")
    }}
  Data <- as.matrix(Data)
  
  dimData <- ncol(Data)
  response <- as.matrix(response)
  n <- nrow(response)
  nrep <- ncol(response)
  
  if(!is.null(mean)){
    if(!length(mean)==n){
      stop("'mean' defined by the user must have the same length as the response variable.")
    }
    mean <- matrix(rep(mean, nrep), ncol=nrep, byrow=F)
    meanModel <- 'userDefined'
  }
  
  y.original <- response
  Data.original <- Data
  
  idxSubset <- NULL
  if(!is.null(m)){
    if(m>n){stop("m cannot be bigger than n.")}
    idxSubset <- sort(sample(x=1:n, size=m, replace=F))
    response <- response[idxSubset,,drop=F]
    Data <- Data[idxSubset,,drop=F]
    if(!is.null(mean)){
      mean <- mean[idxSubset,,drop=F]
    }
  }
  
  
  
  if(is.null(hyper)){
    
    ## lower bounds for candidates
    hyper=list()
    if(any(Cov=='linear')){
      hyper$linear.a=rep(log(1e-4), dimData)
      hyper$linear.i=log(1e-4)
    }
    if(any(Cov=='pow.ex')){
      hyper$pow.ex.v=log(1e-4)
      hyper$pow.ex.w=rep(log(1e-4), dimData)
    }
    if(any(Cov=='matern')){
      hyper$matern.v=log(1e-4)
      hyper$matern.w=rep(log(1e-4), dimData)
    }
    if(any(Cov=='rat.qu')){
      hyper$rat.qu.a=log(1e-4)
      hyper$rat.qu.v=log(1e-4)
      hyper$rat.qu.w=rep(log(1e-4), dimData)
    }
    hyper$vv=log(1e-4)
    hyper_low <- unlist(hyper)
    
    if(!is.null(NewHyper)){
      hyper.nam <- c(names(hyper_low),NewHyper)
      for(i in 1:length(NewHyper)){
        hyper_low <- c(hyper_low, -1)
      }
      names(hyper_low) <- hyper.nam
    }
    
    ## upper bounds for candidates
    hyper=list()
    if(any(Cov=='linear')){
      hyper$linear.a=rep(log(1e4), dimData)
      hyper$linear.i=log(1e4)
    }
    if(any(Cov=='pow.ex')){
      hyper$pow.ex.v=log(1e4)
      hyper$pow.ex.w=rep(log(1e4), dimData)
    }
    if(any(Cov=='matern')){
      hyper$matern.v=log(1e4)
      hyper$matern.w=rep(log(1e4), dimData)
    }
    if(any(Cov=='rat.qu')){
      hyper$rat.qu.a=log(1e4)
      hyper$rat.qu.v=log(1e4)
      hyper$rat.qu.w=rep(log(1e4), dimData)
    }
    hyper$vv=log(1e4)
    hyper_upp <- unlist(hyper)
    
    if(!is.null(NewHyper)){
      hyper.nam <- c(names(hyper_upp),NewHyper)
      for(i in 1:length(NewHyper)){
        hyper_upp <- c(hyper_upp, 1)
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
    hyper=hyper[substr(names(hyper),1,6)%in%c(Cov,'vv')]
  }  
  hp.name=names(unlist(hyper))
  
  if(meanModel==0) {response <- response; mean <- 0}
  if(meanModel==1) {
    mean <- mean(response)
    response <- as.matrix(response-mean)
  }
  meanLinearModel <- NULL
  if(meanModel=='t') {
    trend <- data.frame(yyy=c(response), xxx=rep(c(Data), nrep))
    meanLinearModel <- lm(yyy~xxx, data=trend)
    response <- matrix(resid(meanLinearModel), nrow=nrow(response), byrow=F)
    mean <- matrix(fitted(meanLinearModel), nrow=nrow(response), byrow=F)
  }
  
  if(meanModel=='avg') {
    if(nrep<3){
      stop('Mean function can only be the average across replications when
           there are more than two replications.')
    }
    mean <- apply(response, 1, mean)
    mean <- matrix(rep(mean, nrep), ncol=nrep, byrow=F)
    response <- response - mean
  }
  
  #### Try a number of hp vector and start with the best
  candidates <- matrix(0, nInitCandidates, length(hyper_upp))
  for(iCand in 1:nInitCandidates){
    candidates[iCand,] <- runif(n = length(hyper_upp), min = hyper_low, max = hyper_upp)
  }
  colnames(candidates) <- hp.name
  cat(c('\n','--------- Initialising ---------- \n'))
  resCand <- apply(candidates, 1, function(x){
    gp.loglikelihood2(hyper.p=x, Data=Data,response=response,
                      Cov=Cov,gamma=gamma, nu=nu)})
  
  best_init <- candidates[which.min(resCand),]
  
  trace=round(trace)
  if(trace>0)
    # cat(c('\n','title: -likelihood:',hp.name,'\n'),sep='     ')
    cat(c('iter:  -loglik:',hp.name,'\n'),sep='     ')
  
  if(!useGradient){gp.Dlikelihood2 <- NULL}
  CG0 <- nlminb(start=best_init, objective=gp.loglikelihood2, 
                gradient=gp.Dlikelihood2,
                Data=Data,response=response,Cov=Cov,gamma=gamma,nu=nu,
                control=list(iter.max=itermax,rel.tol=reltol,trace=trace))
  
  # if(trace!=F&CG0$convergence==0)
  #   cat('\n','    optimization finished. Converged.','\n')
  # if(trace!=F&CG0$convergence==1)
  #   cat('\n','    optimization finished. Failed Converge.','\n')
  if(trace>0)
    cat('\n','    optimization finished.','\n')
  CG=CG0[[1]]
  names(CG)=hp.name
  CG.df=data.frame(CG=CG,CG.N=substr(hp.name,1,8))
  names(CG.df)=c('CG','CG.N')
  hyper.cg=split(CG.df$CG,CG.df$CG.N)
  
  nkernels=length(Cov)
  CovList=vector('list',nkernels)
  for(i in 1:nkernels) CovList[i]=list(paste0('cov.',Cov[i]))
  CovL=lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper.cg, Data=Data, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper.cg, Data=Data, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){return(f(hyper=hyper.cg, Data=Data))}
  })
  if(length(CovL)==1)
    Q=CovL[[1]]
  if(length(CovL)>1)
    Q=Reduce('+',CovL)
  
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
      
      D2fx <- lapply(seq_along(hyper.cg),function(i){
        Dp=hyper.cg[i]
        name.Dp=names(Dp)
        f=get(paste0('D2',name.Dp))
        if(name.Dp%in%c('pow.ex.w','pow.ex.v') )
          D2para=f(hyper=hyper.cg, data=Data, gamma=gamma, inv.Q=invQ, Alpha.Q=AlphaQ)
        if(name.Dp%in%c('matern.w','matern.v') )
          D2para=f(hyper=hyper.cg, data=Data, nu=nu, inv.Q=invQ, Alpha.Q=AlphaQ)
        if(!name.Dp%in%c('pow.ex.w','pow.ex.v','matern.w','matern.v') & !name.Dp%in%c('linear.a','linear.i'))
          D2para=f(hyper=hyper.cg, data=Data, inv.Q=invQ, Alpha.Q=AlphaQ)
        if(name.Dp%in%c('linear.a'))
          D2para=f(hyper=hyper.cg, data=Data, Alpha.Q=AlphaQ)
        if(name.Dp%in%c('linear.i'))
          D2para=f(hyper=hyper.cg, inv.Q=invQ, Alpha.Q=AlphaQ)
        return(D2para)
      })
      names(D2fx) <- names(hyper.cg)
      D2fx <- unlist(D2fx)
      D2fxList[[irep]] <- D2fx
    }
    
    D2fx <- Reduce('+', D2fxList)
    
    var_hyper <- (-1/(unlist(D2fx)*dim(Data)[1]))
  }else{
    var_hyper <- NULL
  }
  
  
  
  fitted <- (Q-diag(exp(hyper.cg$vv),dim(Q)[1]))%*%invQ%*%(response)+mean
  fitted.var <- exp(hyper.cg$vv)*rowSums((Q-diag(exp(hyper.cg$vv),dim(Q)[1]))*t(invQ))
  result <- list('hyper'=hyper.cg,'var.hyper'=var_hyper,
                 'fitted.mean'=fitted,
                 fitted.sd=sqrt(fitted.var),
                 'train.x'=Data,'train.y'=response,
                 'train.yOri'=y.original, 
                 'train.DataOri'=Data.original, 
                 'idxSubset'=idxSubset,
                 'CovFun'=Cov,'gamma'=gamma,'nu'=nu,'Q'=Q,'inv'=invQ,'mean'=mean,
                 'meanModel'=meanModel, 'meanLinearModel'=meanLinearModel,
                 conv=CG0$convergence,'hyper0'=hyper)
  class(result) <- 'gpr'
  return(result)
}





#' Prediction using Gaussian Process
#'
#' @inheritParams gpr
#' @param train List resulting from training which is a 'gpr' object. Default to be
#'   NULL. If NULL, do training based on the other given arguments; if TRUE,
#'   other arguments (except for Data.new) will replaced by NULL; if FALSE, only
#'   do prediction based on the other given arguments.
#' @param Data.new The test data. Must be a vector or a matrix.
#' @param noiseFreePred Logical. If TRUE, predictions will be noise-free.
#' @param Y  --------- CHECK   Same as Data????  The input data from train data.
#'   Matrix or vectors are both acceptable. Some data.frames are not acceptable.
#'   --------------
#' @param mSR Subset size m if Subset of Regressors method is used. It must be
#'   smaller than the total sample size.
#'
#' @return A list containing  \describe{ 
#' \item{pred.mean}{Mean of predictions}
#' \item{pred.sd}{Standard deviation of predictions}
#' \item{newdata}{New data}
#' \item{noiseFreePred}{Logical. If TRUE, predictions are noise-free.}
#' \item{...}{Objects of 'gpr' class. }
#'   }
#' @export
#'
gppredict <- function(train=NULL,Data.new=NULL,noiseFreePred=F,hyper=NULL, 
                         Data=NULL, Y=NULL, mSR=NULL,
                         Cov=NULL,gamma=NULL,nu=NULL,meanModel=0,mean=0){
  Data.new=as.matrix(Data.new)
  if(mean==1){
    mean=mean(Y)
  }
  if(!is.null(Data)) Data=as.matrix(Data)
  if(!is.null(Y)) Y=as.matrix(Y-mean)
  if("pow.ex"%in%Cov & is.null(gamma)){
    stop("Argument 'gamma' must be informed for pow.ex kernel")
  }
  if("matern"%in%Cov & is.null(nu)){
    stop("Argument 'nu' must be informed for matern kernel")
  }
  if(class(train)=='gpr'){
    hyper=train$hyper
    Data=train$train.x
    Y=train$train.y
    Cov=train$CovFun
    gamma=train$gamma
    nu=train$nu
    mean=train$mean
    meanModel=train$meanModel
    meanLinearModel=train$meanLinearModel
  }
  
  nrep <- ncol(Y)
  
  if(is.null(train)){
    train=gpr(Data=Data,response=Y,Cov=Cov,hyper=hyper,gamma=gamma,nu=nu)
  }
  if(is.null(Data.new)) Data.new=Data

  nkernels <- length(Cov)


if(is.null(mSR)){
  
  CovList <- vector('list',nkernels)
  for(i in 1:nkernels) CovList[i]=list(paste0('cov.',Cov[i]))
  
  CovL1 <- lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper, Data=Data, Data.new=Data.new, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper, Data=Data, Data.new=Data.new, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){return(f(hyper=hyper, Data=Data, Data.new=Data.new))}
  })
  Q1 <- Reduce('+',CovL1)
  
  Q1 <- t(Q1)  #  #  #
  
  CovL <- lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper, Data=Data, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper, Data=Data, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){return(f(hyper=hyper, Data=Data))}
  })
  Q <- Reduce('+',CovL)
  diag(Q) <- diag(Q)+exp(hyper$vv)
  
  for(i in 1:nkernels) CovList[i] <- list(paste0('diag.',Cov[i]))
  CovLn <- lapply(CovList,function(j){
    f=get(j)
    f(hyper=hyper, data=Data.new)
  })
  Qstar <- Reduce('+',CovLn)

  invQ <- chol2inv(chol(Q))
  
  QQ1 <- invQ%*%t(Q1)
  QR <- invQ%*%Y

  if(meanModel=='t'){
    newtrend <- data.frame(xxx=Data.new[,1])
    mean <- predict(meanLinearModel,newdata=newtrend)
    mu <- Q1%*%QR + matrix(rep(mean, nrep), ncol=nrep, byrow=F)
  }else{
    if(meanModel=='avg'){
      if(!all(Data.new[,1]%in%Data[,1])){
        stop("For predicting response values at input locations Data.new, 
             the mean function can only be 'avg' across replications if
             all Data.new are included in Data.")
      }
      mu <- Q1%*%QR + mean
    }else{
      mu <- Q1%*%QR + mean
    }
  }

  if(noiseFreePred){
    sigma2 <- Qstar-as.matrix(diag(Q1%*%invQ%*%t(Q1)))
  }else{
    sigma2 <- Qstar-as.matrix(diag(Q1%*%invQ%*%t(Q1)))+exp(hyper$vv)
  }
  pred.sd <- sqrt(sigma2[,1])
  
}else{ # mSR set
  
  
  if(mSR>=n){stop("mSR must be smaller than n.")}
  n <- nrow(Data)
  idx <- sort(sample(x=1:n, size=mSR, replace=F))
    
  CovList <- vector('list',nkernels)
  for(i in 1:nkernels) CovList[i]=list(paste0('cov.',Cov[i]))

  Cov_m_ns <- lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper, Data=Data[idx,,drop=F], Data.new=Data.new, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper, Data=Data[idx,,drop=F], Data.new=Data.new, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){return(f(hyper=hyper, Data=Data[idx,,drop=F], Data.new=Data.new))}
  })
  K_m_nstar <- Reduce('+',Cov_m_ns)

  Cov_mm <- lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper, Data=Data[idx,,drop=F], gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper, Data=Data[idx,,drop=F], nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){return(f(hyper=hyper, Data=Data[idx,,drop=F]))}
  })
  K_mm <- Reduce('+',Cov_mm)

  Cov_m_n <- lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper, Data=Data[idx,,drop=F], Data.new=Data, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper, Data=Data[idx,,drop=F], Data.new=Data, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){return(f(hyper=hyper, Data=Data[idx,,drop=F], Data.new=Data))}
  })
  K_m_n <- Reduce('+',Cov_m_n)

  toBeInverted <- K_m_n%*%t(K_m_n) + exp(hyper$vv)*K_mm
  diag(toBeInverted) <- diag(toBeInverted) + 1e-8
  invTerm <- chol2inv(chol(toBeInverted))
  
  centYpred <- t(K_m_nstar)%*%invTerm%*%K_m_n%*%Y
  varYpred <- exp(hyper$vv)*t(K_m_nstar)%*%invTerm%*%K_m_nstar
  range(varYpred)
  
  
  if(meanModel=='t'){
    newtrend <- data.frame(xxx=Data.new[,1])
    mean <- predict(meanLinearModel,newdata=newtrend)
    mu <- centYpred + matrix(rep(mean, nrep), ncol=nrep, byrow=F)
  }else{
    if(meanModel=='avg'){
      if(!all(Data.new[,1]%in%Data[,1])){
        stop("For predicting response values at input locations Data.new, 
             the mean function can only be 'avg' across replications if
             all Data.new are included in Data.")
      }
      mu <- centYpred + mean
    }else{
      mu <- centYpred + mean
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

  result=c(list('pred.mean'=mu,
                'pred.sd'=pred.sd,
                'newdata'=Data.new),
           'noiseFreePred'=noiseFreePred,
           unclass(train))
  class(result)='gpr'
  return(result)
}






################################# tools #####################################
rmse <- function(t,a){ 
  y <- sqrt(sum((a-t)^2)/length(t))
  return(y)
}

########################### derivatives #####################################

Dloglik.linear.a <- function(hyper,data,AlphaQ){
  Dlinear.aj=sapply(1:ncol(data),function(i){
    Xi <- as.matrix(data[,i])
    A1 <- as.matrix(1)
    DPsi <- exp(hyper$linear.a[i])*DistMatLinear_sq(X=Xi, A=A1)
    res <- sum(diag(AlphaQ%*%DPsi))
    return(res)
  })
  return(Dlinear.aj)
}

Dloglik.linear.i <- function(hyper,Alpha,invQ){
  constMat <- matrix(exp(hyper$linear.i), nrow=nrow(invQ), ncol=ncol(invQ))
  Dfa0 <- 0.5*sum(diag(  (Alpha%*%t(Alpha) - invQ)%*%constMat ))
  return(Dfa0)
}

Dloglik.pow.ex.w <- function(hyper,data,AlphaQ,gamma){
  Dpow.ex.wj=sapply(1:ncol(data),function(i){
    Xi <- as.matrix(data[,i])
    A1 <- as.matrix(1)
    DPsi <- -cov.pow.ex(hyper=hyper, Data=data, gamma=gamma)*DistMat_sq(X=Xi, A=A1, power=gamma)*exp(hyper$pow.ex.w[i])
    res <- 0.5*sum(diag(AlphaQ%*%DPsi))
    return(res)
  })
  return(Dpow.ex.wj)
}

Dloglik.pow.ex.v <- function(hyper,data,AlphaQ,gamma){
  DDpow.ex.v=cov.pow.ex(hyper,data, gamma=gamma)
  Dpow.ex.v=0.5*sum(diag(AlphaQ%*%DDpow.ex.v))
  return(Dpow.ex.v)
}

#######################

Dloglik.matern.v <- function(hyper,data,AlphaQ,nu){
  DDmatern.v <- cov.matern(hyper=hyper, Data=data, Data.new=NULL, nu=nu)
  Dmatern.v <- 0.5*sum(diag(AlphaQ%*%DDmatern.v))
  return(Dmatern.v)
}

Dloglik.matern.w <- function(hyper,data,AlphaQ,nu){
  
  data <- as.matrix(data)
  cc <- exp(hyper$matern.v)
  dimData <- ncol(data)
  if(dimData==1){
    A <- as.matrix(exp(hyper$matern.w))
  }else{
    A <- diag((exp(hyper$matern.w)))
  }
  
  dist_C <- DistMat_sq(X=data, A=A, power=2)
  dist_D <- sqrt(dist_C)
  
  Dmatern.wj <- sapply(1:ncol(data),function(i){
    Xi <- as.matrix(data[,i])
    A1 <- as.matrix(1)
    dist_i <- DistMat_sq(X=Xi, A=A1, power=2)
    if(nu==3/2){
      DPsi <- -1.5*cc*exp(hyper$matern.w[i])*dist_i*exp(-sqrt(3)*dist_D)
    }
    if(nu==5/2){
      DPsi <- cc*exp(hyper$matern.w[i])*dist_i*exp(-sqrt(5)*dist_D)*(5/6)*(-1-sqrt(5)*dist_D)
    }
    res <- 0.5*sum(diag(AlphaQ%*%DPsi))
    return(res)
  })
  return(Dmatern.wj)
}


Dloglik.rat.qu.w <- function(hyper,data,AlphaQ){

  dimData <- ncol(data)
  d1list <- vector("list", dimData)
  for(i in 1:dimData){
    covmatrix <- cov.rat.qu(hyper=hyper, Data=data)
    
    Xi <- as.matrix(data[,i])
    A1 <- as.matrix(1)
    hyperAdj <- hyper
    hyperAdj$rat.qu.v <- 0
    covmatrixAdj_a <- cov.rat.qu(hyper=hyperAdj, Data=data)
    hyperAdj$rat.qu.a <- 0
    covmatrixAdj_0 <- cov.rat.qu(hyper=hyperAdj, Data=data)
    
    d1 <- -exp(hyper$rat.qu.a)*exp(hyper$rat.qu.v)*covmatrixAdj_a*covmatrixAdj_0*DistMat_sq(X=Xi, A=A1, power=2)*exp(hyper$rat.qu.w[i])
    d1list[[i]] <- d1
  }

  Drat.qu.wj <- sapply(1:dimData, function(i){sum(0.5*AlphaQ*d1list[[i]])} )
  return(Drat.qu.wj)
}

Dloglik.rat.qu.a <- function(hyper,data,AlphaQ){
  
  covmatrix=cov.rat.qu(hyper,data)
  
  hyper <- lapply(hyper,exp)
  Data <- as.matrix(data)
  Q <- ncol(Data)
  if(Q==1){
    A <- as.matrix(hyper$rat.qu.w)
  }else{
    A <- diag(hyper$rat.qu.w)
  }
  v.power <- DistMat_sq(X=Data,A=A,power=2)
  log_term <- log( 1 + v.power )
  
  DDrat.qu.a <- log_term*covmatrix*(-hyper$rat.qu.a)
  Drat.qu.a <- 0.5*sum(diag(AlphaQ%*%DDrat.qu.a))

  return(Drat.qu.a)
}


Dloglik.rat.qu.v <- function(hyper,data,AlphaQ){
  DDrat.qu.v <- cov.rat.qu(hyper,data)
  Drat.qu.v <- 0.5*sum(diag(AlphaQ%*%DDrat.qu.v))
  return(Drat.qu.v)
}



Dloglik.vv <- function(hyper,Alpha,invQ){
  Dfvv <- 0.5*sum(diag(Alpha%*%t(Alpha) - invQ))*exp(hyper$vv)
  return(Dfvv)
}


D2linear.a <- function(hyper,data,Alpha.Q){
  D2linear.aj=sapply(1:ncol(data),function(i){
    Xi <- as.matrix(data[,i])
    A1 <- as.matrix(1)
    DPsi <- exp(hyper$linear.a[i])*DistMatLinear_sq(X=Xi, A=A1)
    res <- sum(diag(Alpha.Q%*%DPsi))
    return(res)
  })
  return(D2linear.aj)
}


D2linear.i <- function(hyper,inv.Q, Alpha.Q){
  constMat <- matrix(exp(hyper$linear.i), nrow=nrow(inv.Q), ncol=ncol(inv.Q))
  D2flinear.i <- D2(constMat, constMat,inv.Q,Alpha.Q)
  return(D2flinear.i)
}



D2pow.ex.w <- function(hyper,data,gamma,inv.Q,Alpha.Q){
  
  dimData <- ncol(data)
  d1list <- vector("list", dimData)
  d2list <- vector("list", dimData)
  for(i in 1:dimData){
    covmatrix <- cov.pow.ex(hyper=hyper, Data=data, gamma=gamma)
    Xi <- as.matrix(data[,i])
    A1 <- as.matrix(1)
    d1 <- -covmatrix*DistMat_sq(X=Xi, A=A1, power=gamma)*exp(hyper$pow.ex.w[i])
    d1list[[i]] <- d1
    d2 <- covmatrix*(DistMat_sq(X=Xi, A=A1, power=2*gamma)*exp(2*hyper$pow.ex.w[i]) - 
                       DistMat_sq(X=Xi, A=A1, power=gamma)*exp(hyper$pow.ex.w[i]))
    d2list[[i]] <- d2
  }

  D2pow.ex.wj=sapply(1:dimData, function(i){
    D2(d1=d1list[[i]], 
       d2=d2list[[i]], 
       inv.Q = inv.Q, Alpha.Q = Alpha.Q)} )

  return(D2pow.ex.wj)
}


D2pow.ex.v <- function(hyper,data,gamma,inv.Q,Alpha.Q){
  DDpow.ex.v <- cov.pow.ex(hyper=hyper, Data=data,Data.new=NULL, gamma=gamma)
  D2pow.ex.v <- D2(DDpow.ex.v,DDpow.ex.v,inv.Q,Alpha.Q)  
  return(D2pow.ex.v)
}

D2matern.v <- function(hyper,data,nu,inv.Q,Alpha.Q){
  DDmatern.v <- cov.matern(hyper=hyper, Data=data, Data.new=NULL, nu=nu)
  D2matern.v <- D2(DDmatern.v,DDmatern.v,inv.Q,Alpha.Q)  
  return(D2matern.v)
}



D2matern.w <- function(hyper,data,nu,inv.Q,Alpha.Q){
  
  data <- as.matrix(data)
  cc <- exp(hyper$matern.v)
  dimData <- ncol(data)
  if(dimData==1){
    A <- as.matrix(exp(hyper$matern.w))
  }else{
    A <- diag((exp(hyper$matern.w)))
  }

  dist_C <- DistMat_sq(X=data, A=A, power=2)
  dist_D <- sqrt(dist_C)

  d1list <- vector("list", dimData)
  d2list <- vector("list", dimData)
  
  for(i in 1:dimData){
    Xi <- as.matrix(data[,i])
    A1 <- as.matrix(1)
    dist_i <- DistMat_sq(X=Xi, A=A1, power=2)
    
    if(nu==3/2){
      d1list[[i]] <- -1.5*cc*exp(hyper$matern.w[i])*dist_i*exp(-sqrt(3)*dist_D)
      d2list[[i]] <- d1list[[i]]*(1 - 0.5*sqrt(3)*dist_C^(-0.5))*exp(hyper$matern.w[i])*dist_i
      diag(d2list[[i]]) <- 0
    }
    
    if(nu==5/2){
      d1list[[i]] <- cc*exp(hyper$matern.w[i])*dist_i*exp(-sqrt(5)*dist_D)*(5/6)*(-1-sqrt(5)*dist_D)
      d2list[[i]] <- (-5/6)*cc*dist_i*exp(hyper$matern.w[i])*exp(-sqrt(5)*dist_D)*(
        1+sqrt(5)*dist_D - 2.5*exp(hyper$matern.w[i])*dist_i)
    }
    
  }
  D2matern.wj=sapply(1:dimData, function(i){
    D2(d1=d1list[[i]], 
       d2=d2list[[i]], 
       inv.Q = inv.Q, Alpha.Q = Alpha.Q)} )
  
  return(D2matern.wj)
  
}


######################################################################

# covmatrixAdj_a will use unit signal variance v and usual alpha parameter
# covmatrixAdj_0 will use unit signal variance v and alpha=1 (ie, log(alpha)=0, exp(log(1)) = exp(0))
D2rat.qu.w <- function(hyper,data,inv.Q,Alpha.Q){
  
  dimData <- ncol(data)
  d1list <- vector("list", dimData)
  d2list <- vector("list", dimData)
  for(i in 1:dimData){
    covmatrix <- cov.rat.qu(hyper=hyper, Data=data)
    
    Xi <- as.matrix(data[,i])
    A1 <- as.matrix(1)
    # calc Drat.qu
    hyperAdj <- hyper
    hyperAdj$rat.qu.v <- 0
    covmatrixAdj_a <- cov.rat.qu(hyper=hyperAdj, Data=data)
    hyperAdj$rat.qu.a <- 0
    covmatrixAdj_0 <- cov.rat.qu(hyper=hyperAdj, Data=data)
    
    d1 <- -exp(hyper$rat.qu.a)*exp(hyper$rat.qu.v)*covmatrixAdj_a*covmatrixAdj_0*DistMat_sq(X=Xi, A=A1, power=2)*exp(hyper$rat.qu.w[i])
    d1list[[i]] <- d1
    
    # calc D2rat.qu
    d2 <- -exp(hyper$rat.qu.a)*exp(hyper$rat.qu.v)*covmatrixAdj_a*covmatrixAdj_0*(
      (-exp(hyper$rat.qu.a)-1)*covmatrixAdj_0*DistMat_sq(X=Xi, A=A1, power=4)*exp(2*hyper$rat.qu.w[i]) + 
        DistMat_sq(X=Xi, A=A1, power=2)*exp(hyper$rat.qu.w[i])
    )
    
    d2list[[i]] <- d2
  }
  
  D2rat.qu.wj=sapply(1:dimData, function(i){
    D2(d1=d1list[[i]], 
       d2=d2list[[i]], 
       inv.Q = inv.Q, Alpha.Q = Alpha.Q)} )
  
  return(D2rat.qu.wj)
}


D2rat.qu.a <- function(hyper,data,inv.Q,Alpha.Q){
  
  covmatrix <- cov.rat.qu(hyper,data)

  hyper <- lapply(hyper,exp)
  Data <- as.matrix(data)
  Q <- ncol(Data)
  if(Q==1){
    A <- as.matrix(hyper$rat.qu.w)
  }else{
    A <- diag(hyper$rat.qu.w)
  }
  v.power <- DistMat_sq(X=Data, A=A, power=2)
  log_term <- log( 1 + v.power )
  
  d1 <- log_term*covmatrix*(-hyper$rat.qu.a)
  d2 <- (d1 + covmatrix)*(-hyper$rat.qu.a)*log_term

  D2rat.qu.a <- D2(d1,d2,inv.Q,Alpha.Q)
  
  return(D2rat.qu.a)
}

D2rat.qu.v <- function(hyper,data,inv.Q,Alpha.Q){
  covmatrix <- cov.rat.qu(hyper,data)
  D2rat.qu.v <- D2(covmatrix,covmatrix,inv.Q,Alpha.Q)  
  return(D2rat.qu.v)
}


D2vv <- function(hyper,data,inv.Q,Alpha.Q){
  D2fvv=D2(diag(exp(hyper$vv),dim(data)[1]),diag(exp(hyper$vv),dim(data)[1]),inv.Q,Alpha.Q)
  return(D2fvv)
}


#'Second derivative of the likelihood
#'
#'Calculates the second derivative of the likelihood function with respect to
#'one of the hyperparameters, given the first and second derivative of the
#'kernel with respect to that hyperparameter.
#'
#'@param d1 First derivative of the kernel function with respect to the required
#'  hyper-parameter.
#'@param d2 Second derivative of the kernel function with respect to the
#'  required hyper-parameter.
#'@param inv.Q Inverse matrix of the covariance matrix
#'@param Alpha.Q  This is iQY %*% t(iQY)-iQ, where iQ is the inverse of the
#'  covariance matrix, Y is the response.
#'
#'@details The function is to calculate the second derivative of the normal
#'  likelihood, using the first and second derivative of the kernel functions.
#'  The first and second derivative need to be pre-defined, for example of
#'  customized covariance function, see "demo('co2')".
#'@return A number
#'
#'@references Shi, J. Q., and Choi, T. (2011), ``Gaussian Process Regression
#'  Analysis for Functional Data'', CRC Press.
#'
#'@export
#'
D2 <- function(d1,d2,inv.Q,Alpha.Q){
  Aii=t(d1)%*%inv.Q%*%d1
  al=Alpha.Q+inv.Q
  return(0.5*(sum(Alpha.Q*(d2-Aii))-sum(al*Aii)))
}

diag.linear <- function(hyper,data){
  Qstar=exp(hyper$linear.i) + data^2%*%matrix(exp(hyper$linear.a))
  return(Qstar)
}

diag.pow.ex <- function(hyper,data){
  Qstar=rep(exp(hyper$pow.ex.v),dim(data)[1])
  return(Qstar)
}

diag.matern <- function(hyper,data){
  Qstar=rep(exp(hyper$matern.v),dim(data)[1])
  return(Qstar)
}

diag.rat.qu <- function(hyper,data){
  Qstar=rep(exp(hyper$rat.qu.v),dim(data)[1])
  return(Qstar)
}


#' Plot Gaussian Process regression -- training and prediction
#'
#' Plot Gaussian Process for a given an object of class 'gpr'.
#'
#' @param x The 'gpr' object from either training or predicting of the Gaussian
#'   Process.
#' @param fitted Logical. Plot fitted value or not. Default to FALSE, which is
#'   to plot the predictions.
#' @param col.no Column number of the input matrix. If the input matrix has more
#'   than one columns, than one of them will be used in the plot. Default to be
#'   the first one.
#' @param ylim Range value for y-axis.
#' @param realisation Which realisation should be plotted (if there are multiple).
#' @param ... Graphical parameters passed to plot().
#' @importFrom  graphics polygon
#' @importFrom  graphics points
#' @importFrom  graphics matpoints
#' @importFrom  graphics matlines
#' @importFrom  graphics lines
#' @importFrom  grDevices rgb
#' @return A plot
#' @export
#'
plot.gpr <- function(x,fitted=F,col.no=1, ylim=NULL, realisation=NULL, ...){
  obj=x
  if(fitted==T){
    if(is.null(obj$fitted.mean)){
      warning('fitted values not found, ploting predicted values')
      type='Prediction'
      mu=obj$pred.mean
      sd=obj$pred.sd
      x=obj$newdata
      X=obj$train.x
      Y=obj$train.yOri
    }
    if(!is.null(obj$fitted.mean)){
      type='Fitted values'
      mu=obj$fitted.mean
      sd=obj$fitted.sd
      X=obj$train.x
      Y=obj$train.yOri
      x=X
    }
  }else{
    if(is.null(obj$pred.mean)){
      warning('predicted values not found, ploting fitted values')
      type='Fitted values'
      mu=obj$fitted.mean
      sd=obj$fitted.sd
      X=obj$train.x
      Y=obj$train.yOri
      x=X
    }
    if(!is.null(obj$pred.mean)){
      type='Prediction'
      mu=obj$pred.mean
      sd=obj$pred.sd
      x=obj$newdata
      X=obj$train.x
      Y=obj$train.yOri
    }
  }
  if(dim(X)[1]<=150|length(X)<=150){
    pchType=4
    PcexNo=0.8
    LcexNo=1.5
  }else{
    pchType=20
    PcexNo=0.1
    LcexNo=0.8
  }
  noiseFreePred <- obj$noiseFreePred
  if(type=='Prediction'){
    if(noiseFreePred){
      type <- "Noise-free prediction"
    }
  }
  
  if(!is.null(realisation)){
    mu <- mu[,realisation]
    Y <- Y[,realisation]
  }
  upper=mu+1.96*(sd);
  lower=mu-1.96*(sd);
  if(is.null(ylim)){
    ylim <- range(upper,lower,Y)
  }
  plot(-100,-100,col=0,xlim=range(X[,col.no],x[,col.no]),ylim=ylim,main=type, xlab="input ",ylab="response",...)
  #
  polygon(c(x[,col.no], rev(x[,col.no])), c(upper, rev(lower)),col = rgb(127,127,127,120, maxColorValue = 255), border = NA)
  #
  points(X[,col.no],Y,pch=pchType,col=2,cex=PcexNo)
  # lines(X[,1],Y)
  lines(x[,col.no],mu,col=4,lwd=LcexNo)  
}


# ########################### likelihood ######################################


gp.loglikelihood2 <- function(hyper.p,Data, response,Cov,gamma,nu){
  
  Data=as.matrix(Data)
  datadim=dim(Data)
  
  hp.class=substr(names(hyper.p),1,8)
  kernel.class=unique(substr(names(hyper.p),1,6))
  hp.class=data.frame(class=hp.class,hp=hyper.p)
  names(hp.class)=c('class','hp')
  hp.list=split(hp.class$hp,hp.class$class)
  hyper.p=hp.list
  
  nkernels=length(Cov)
  CovList=vector('list',nkernels)
  for(i in 1:nkernels) CovList[i]=list(paste0('cov.',Cov[i]))
  CovL=lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper.p, Data=Data, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper.p, Data=Data, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){return(f(hyper=hyper.p, Data=Data))}
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
    fX <- nrep*0.5*logdetQ + 0.5*sum(diag( tresp_invQ_resp )) + nrep*0.5*n*log(2*pi)
  }
  fX <- as.numeric(fX)

  return(fX)
}

# ########################### gradient ######################################

gp.Dlikelihood2 <- function(hyper.p,  Data, response, Cov, gamma, nu){
  
  Data=as.matrix(Data)
  datadim=dim(Data);
  
  hp.class=substr(names(hyper.p),1,8)
  kernel.class=unique(substr(names(hyper.p),1,6))
  hp.class=data.frame(class=hp.class,hp=hyper.p)
  names(hp.class)=c('class','hp')
  hyper.p=split(hp.class$hp,hp.class$class)
  
  nkernels=length(Cov)
  CovList=vector('list',nkernels)
  for(i in 1:nkernels) CovList[i]=list(paste0('cov.',Cov[i]))
  CovL=lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper.p, Data=Data, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper.p, Data=Data, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){return(f(hyper=hyper.p, Data=Data))}
  })
  Q <- Reduce('+',CovL)
  diag(Q) <- diag(Q)+exp(hyper.p$vv) + 1e-8

  invQ <- chol2inv(chol(Q))
  
  nrep <- ncol(response)

DfxList <- vector('list',nrep)
for(irep in 1:nrep){
  Alpha <- invQ%*%as.matrix(response[,irep])
  AlphaQ <- Alpha%*%t(Alpha)-invQ
  
  Dfx=lapply(seq_along(hyper.p),function(i){
    Dp=hyper.p[i];
    name.Dp=names(Dp)
    f=get(paste0('Dloglik.',name.Dp))
    if(name.Dp%in%c('pow.ex.w','pow.ex.v') )
      Dpara=f(hyper=hyper.p, data=Data, AlphaQ=AlphaQ, gamma=gamma)
    if(name.Dp%in%c('matern.w','matern.v') )
      Dpara=f(hyper=hyper.p, data=Data, AlphaQ=AlphaQ, nu=nu)
    if(name.Dp%in%c('rat.qu.w','rat.qu.v','rat.qu.a') )
      Dpara=f(hyper=hyper.p, data=Data, AlphaQ=AlphaQ)
    if(name.Dp%in%c('linear.a') )
      Dpara=f(hyper=hyper.p, data=Data, AlphaQ=AlphaQ)
    if(name.Dp%in%c('linear.i') )
      Dpara=f(hyper=hyper.p, Alpha=Alpha, invQ=invQ)
    if(name.Dp=='vv')
      Dpara=f(hyper=hyper.p, Alpha=Alpha, invQ=invQ)
    if(substr(name.Dp, 1, 6)=='custom')
      Dpara=f(hyper=hyper.p, data=Data, AlphaQ=AlphaQ)
    return(Dpara)
  })
  
  names(Dfx)=names(hyper.p)
  Dfx <- - unlist(Dfx) # Minus loglik
  DfxList[[irep]] <- Dfx
}
  
Dfx <- Reduce('+', DfxList)
  
  return(Dfx)
}




listFunctions <- function(filename) {
  temp.env <- new.env()
  sys.source(filename, envir = temp.env)
  functions <- lsf.str(envir=temp.env)
  rm(temp.env)
  return(functions)
}
# listFunctions("gp.functions5NEW_1.R") # DO NOT UNCOMMENT, run in console instead


######################## Covariance functions ############################

cov.pow.ex <- function(hyper,Data,Data.new=NULL,gamma){
  
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
# rat.qu.s will no longer be used

# source("CovMaternCpp_sq.R")
# source("CovMaternCpp.R")

# distmat <- sqrt(2*nu)*sqrt(DistMat_sq(X=Data, A=A, power=2))
# covmaternNew <- cc*2^(1-nu)/(gamma(nu))*(distmat^nu)*besselK(x=distmat, nu=nu)
# # besselK at zero produces inf. At zero covmatern=cc
# range(covmaternNew)
# corner(covmaternNew)
# corner(covmatern0)
# max(abs(covmatern0-covmaternNew))

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




# 

#k(x^1,x^2)=sum_{q=1}^{Q} a_q * x^1_q * x^2_q
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

######################## Predict ############################
gppredictNEW <- function(train=NULL,Data.new=NULL,noiseFreePred=F,hyper=NULL, 
                         Data=NULL, Y=NULL, mSR=NULL,
                         Cov=NULL,gamma=NULL,nu=NULL,meanModel=NULL,mean=0){
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
    train=gprNEW(Data=Data,response=Y,Cov=Cov,hyper=hyper,gamma=gamma,nu=nu)
  }
  if(is.null(Data.new)) Data.new=Data
  # n=dim(Data)[1];
  # nn=dim(Data.new)[1];
  # n.var=dim(Data)[2]
  
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
  
  ### K_m_nstar
  Cov_m_ns <- lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper, Data=Data[idx,,drop=F], Data.new=Data.new, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper, Data=Data[idx,,drop=F], Data.new=Data.new, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){return(f(hyper=hyper, Data=Data[idx,,drop=F], Data.new=Data.new))}
  })
  K_m_nstar <- Reduce('+',Cov_m_ns)
  
  ### K_m_m
  Cov_mm <- lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper, Data=Data[idx,,drop=F], gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper, Data=Data[idx,,drop=F], nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){return(f(hyper=hyper, Data=Data[idx,,drop=F]))}
  })
  K_mm <- Reduce('+',Cov_mm)
  # diag(Q_mm) <- diag(Q_mm) + exp(hyper$vv)
  # diag(K_mm) <- diag(K_mm) + 1e-8
  
  ### K_m_n
  Cov_m_n <- lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex'){return(f(hyper=hyper, Data=Data[idx,,drop=F], Data.new=Data, gamma=gamma))}
    if(j=='cov.matern'){return(f(hyper=hyper, Data=Data[idx,,drop=F], Data.new=Data, nu=nu))}
    if(!(j%in%c('cov.pow.ex', 'cov.matern'))){return(f(hyper=hyper, Data=Data[idx,,drop=F], Data.new=Data))}
  })
  K_m_n <- Reduce('+',Cov_m_n)
  # whichAdj <- cbind(1:mSR, idx)
  # K_m_n[whichAdj] <- K_m_n[whichAdj] + 1e-6
  
  
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
  

  
  result=c(list('noiseFreePred'=noiseFreePred, 
                # 'pred.mean'=mu[,1],
                'pred.mean'=mu,
                'pred.sd'=pred.sd)
                # ,'newdata'=Data.new,unclass(train),nsigma=any(Qstar<as.matrix(colSums(t(Q1)*QQ1))))
           )
  class(result)='gpr'
  return(result)
}

gprNEW <- function(Data, response, Cov='pow.ex', 
                m = NULL, hyper=NULL, NewHyper=NULL, meanModel=0, mean=NULL, 
                gamma=NULL, nu=NULL,
                useGradient=T, itermax=100,reltol=8e-10,trace=0,
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
      hyper$matern.w=rep(log(1e-4), dimData)   # cannot be too small => causes NegDef
    }
    if(any(Cov=='rat.qu')){
      hyper$rat.qu.a=log(1e-4)
      hyper$rat.qu.v=log(1e-4)
      hyper$rat.qu.w=rep(log(1e-4), dimData)
    }
    hyper$vv=log(1e-4)
    hyper_low <- unlist(hyper)

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
    
    if(length(hyper_upp)!=length(hyper_low)){
      stop("hyper_upp and hyper_low must have the same dimension")
    }
    
    
    hyper.nam <- names(hyper_low)
    
    # # revise NewHyper later ----
    # if(!is.null(NewHyper)){
    #   hyper.nam=c(hyper.nam,NewHyper)
    #   nh.length=length(NewHyper)
    #   for(i in 1:nh.length){
    #     hyper=c(hyper,runif(1,-1,1))
    #   }
    #   names(hyper)=hyper.nam
    # }
    # #----------------------------
  }
  
  if(!is.null(hyper)){
    hyper=hyper[substr(names(hyper),1,6)%in%c(Cov,'vv')]
  }  
  hp.name=names(unlist(hyper))
  
  if(meanModel==0) {response <- response; mean <- 0; meanModel <- 0}
  if(meanModel==1) {
    mean <- mean(response)
    response <- as.matrix(response-mean)
  }
  meanLinearModel <- NULL
  if(meanModel=='t') {
    # trend <- data.frame(yyy=response, xxx=Data[,1])
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
    # trend <- data.frame(yyy=response, xxx=Data[,1])
    mean <- apply(response, 1, base::mean)
    mean <- matrix(rep(mean, nrep), ncol=nrep, byrow=F)
    response <- response - mean
  }
  
  #### Try a number of hp vector and start with the best
  # each row is a hp vector
  
  candidates <- matrix(0, nInitCandidates, length(hyper_upp))
  for(iCand in 1:nInitCandidates){
    candidates[iCand,] <- runif(n = length(hyper_upp), min = hyper_low, max = hyper_upp)
  }
  colnames(candidates) <- hyper.nam
  cat(c('\n','--------- Initialising ---------- \n'))
  resCand <- apply(candidates, 1, function(x){
    # cat(paste0(round(unlist(x),4)), paste0("\n"))
    gp.loglikelihood2NEW(hyper.p=x, Data=Data,response=response,
                         Cov=Cov,gamma=gamma, nu=nu)})
  
  # print(resCand)
  
  best_init <- candidates[which.min(resCand),]
  ###
  
  trace=round(trace)
  if(trace>0)
    # cat(c('\n','title: -likelihood:',hp.name,'\n'),sep='     ')
    cat(c('iter:  -loglik:',hp.name,'\n'),sep='     ')

  if(!useGradient){gp.Dlikelihood2NEW <- NULL}
  CG0 <- nlminb(start=best_init, objective=gp.loglikelihood2NEW, 
                gradient=gp.Dlikelihood2NEW,
                Data=Data,response=response,Cov=Cov,gamma=gamma,nu=nu,
                control=list(iter.max=itermax,rel.tol=reltol,trace=trace))
  
  # CG0 <- nlminb(unlist(hyper), gp.loglikelihood2, gp.Dlikelihood2,Data=Data,response=response,Cov=Cov,gamma=gamma,control=list(iter.max=itermax,rel.tol=reltol,trace=trace))
  # CG0 <- optim(unlist(hyper), gp.loglikelihood2, gp.Dlikelihood2,Data=Data,response=response,Cov=Cov,gamma=gamma,method='CG',control=list(maxiter=itermax,reltol=reltol,trace=trace))
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

  # response <- as.matrix(response)
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
      D2fx <- unlist(D2fx)  # minus?
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
                 # 'fitted.mean'=fitted[,1],
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

# ########################### likelihood ######################################
# gp.loglikelihood2 <- function(hyper.p,Data, response,Cov,gamma=1){
#   #this function doesn't return anything, it's for the conjugate gradian
#   #hyper is a list of hyper-parameters
#   #Data should have the form that, each column is a variable
#   #response is the given response vector
#   #Cov is a function contains all the covariance matrix, defult is:
#   ###cov.linear(hyper,Data)+cov.pow.ex(hyper,Data), but 
#   #####it could be other forms.
#   
#   Data=as.matrix(Data)
#   datadim=dim(Data)
# 
#   hp.class=substr(names(hyper.p),1,8)
#   kernel.class=unique(substr(names(hyper.p),1,6))
#   hp.class=data.frame(class=hp.class,hp=hyper.p)
#   names(hp.class)=c('class','hp')
#   hp.list=split(hp.class$hp,hp.class$class)
#   hyper.p=hp.list
# 
#   n=length(Cov)
#   CovList=vector('list',n)
#   for(i in 1:n) CovList[i]=list(paste0('cov.',Cov[i]))
#   CovL=lapply(CovList,function(j){
#     f=get(j)
#     if(j=='cov.pow.ex')
#       return(f(hyper.p,Data,Data,gamma=gamma))
#     if(j!='cov.pow.ex')
#       return(f(hyper.p,Data,Data))
#   }  )
#   Q=Reduce('+',CovL)
#   diag(Q)=diag(Q)+exp(hyper.p$vv)
# 
#   response=as.matrix(response)
#   # invQ=pseudoinverse(Q+diag(1e-9,ncol=ncol(Q),nrow=nrow(Q)))
#   # invQ=pseudoinverse(Q)
#   # invQ=mymatrix2(Q)$res
#   invQ <- chol2inv(chol(Q))
#   invQ.response=invQ%*%response
#   logdetQ=sum(determinant(Q,logarithm=T)$modulus)
#   
#   # fX=c(0.5*logdetQ + 0.5*t(response)%*%invQ%*%response + 0.5*dim(Data)[1]*log(2*pi))
#   fX=0.5*logdetQ + 0.5*t(response)%*%invQ.response + 0.5*nrow(Data)*log(2*pi)
#   
#   # temp=0
#   # if(any(is.na(Xprior2))==F){
#   # if(any(Cov=='linear')){
#   #   for (d in 1:n.hyper){
#   #     temp=temp+hyper$linear.a[d]+0.5*((hyper$linear.a[d]-Xprior2$mua[d])/Xprior2$sigma[d])^2}}
#   #     # updating (mu, sigma, log(a))
#   # if(any(Cov=='pow.ex')){
#   # 	for (d in 1:n.hyper){
#   # 		temp=temp+(Xprior2$linear.alpha[d]+1)*hyper$w[d]+Xprior2$mu[d]/(Xprior2$linear.alpha[d]*exp(hyper$w[d]))}
#   # 		# updating (alpha, mu, log(w))
#   #     temp=temp+hyper$v1[1]+0.5*((hyper$v1[1]-Xprior2$muv1[1])/Xprior2$sigmav1[1])^2
#   # 	# updating (muv1[1],sigmav1[1], v1[1])
#   #     temp=temp+hyper$v0+0.5*((hyper$v0-Xprior2$muv0[1])/Xprior2$sigmav0[1])^2
#   # 	# updating (muv0[1],sigmav0[1], v0[1])
#   # }}
#   # 
#   # fX=fX+temp
#   return(fX)
# }
# 
# gp.Dlikelihood2 <- function(hyper.p,  Data, response,Cov,gamma){
#   #this function doesn't return anything, it's for the conjugate gradian
#   #hyper is a list of hyper-parameters
#   #Data should have the form that, each column is a variable
#   #response is the given response vector
#   #Cov is a function contains all the covariance matrix, defult is:
#   ###cov.linear(hyper,Data)+cov.pow.ex(hyper,Data), but 
#   #####it could be other forms.
#   
#   Data=as.matrix(Data)
#   datadim=dim(Data);
#   
#   hp.class=substr(names(hyper.p),1,8)
#   kernel.class=unique(substr(names(hyper.p),1,6))
#   hp.class=data.frame(class=hp.class,hp=hyper.p)
#   names(hp.class)=c('class','hp')
#   hyper.p=split(hp.class$hp,hp.class$class)
#   
#   n=length(Cov)
#   CovList=vector('list',n)
#   for(i in 1:n) CovList[i]=list(paste0('cov.',Cov[i]))
#   CovL=lapply(CovList,function(j){
#     f=get(j)
#     f(hyper.p,Data,Data)
#   }  )
#   Q=Reduce('+',CovL)
#   
#   diag(Q)=diag(Q) + exp(hyper.p$vv)
# 
#   response=as.matrix(response)
#   # invQ=pseudoinverse(Q+diag(1e-9,ncol=ncol(Q),nrow=nrow(Q)))
#   # invQ=pseudoinverse(Q)
#   # invQ=mymatrix2(Q)$res
#   invQ <- chol2inv(chol(Q))
#   Alpha=invQ%*%response
#   #   Alpha2=t(Alpha)%*%Alpha
#   AlphaQ=Alpha%*%t(Alpha)-invQ
# 
#   Dfx=lapply(seq_along(hyper.p),function(i){
#     Dp=hyper.p[i];
#     name.Dp=names(Dp)
#     f=get(paste0('Dloglik.',name.Dp))
#     if(name.Dp%in%c('pow.ex.w','pow.ex.v') )
#       Dpara=f(hyper.p,Data,AlphaQ,gamma=gamma)
#     if(name.Dp=='vv')
#       Dpara=f(hyper.p,Alpha,invQ)
#     if(!name.Dp%in%c('pow.ex.w','pow.ex.v','vv'))
#       Dpara=f(hyper.p,Data,AlphaQ)
#     return(Dpara)
#   })
#   
#   names(Dfx)=names(hyper.p)
#   Dfx=-0.5*unlist(Dfx)
#   Dfx
# }






################################# tools ###########################################
rmse <- function(t,a){ 
#compute the root mean squar error between two vectors
y = sqrt(sum((a-t)^2)/length(t))
return(y)
}

# mymatrix <- function(matrix,log=T){
# #singular decomposition 
# 	m=matrix
# 	a=svd(m)
# 	U=a$u; V=a$v; D=a$d
# 	l=length(D)
# 	idx=which(D<2e-13)
# 	if(length(idx)>0)
# 		D=D+1e-12
# 	else D=D
# 	inv=V%*%diag(1/D)%*%t(U)
# 	if(log==T)
# 	  det=sum(log(D))
# 	else
# 		det=prod(D)
# 	return(list("inv"=inv,"det"=det))
# }


# mymatrix2 <- function(smatrix,sB='sB',det=F,log=T,jitter=1e-10){
#   mat=smatrix+diag(jitter,dim(smatrix)[1])
#   smatrix=as.spam(mat,eps=1e-8)
#   if(is.character(sB)) sB=diag(1,dim(mat)[1])
#   else sB=as.matrix(sB)
#   sB=as.spam(sB,eps=1e-8)
#   x=solve.spam(smatrix,sB)
#   d=NULL
#   if(det==T){
#     L=chol(smatrix)
#     if(log==T)
#       d=2*sum(log(diag(L)))
#     if(log==F)
#       d=prod(diag(L))
#   }
#   return(list('res'=as.matrix(x),'det'=d))
# }

# xixj <- function(mat,mat.new=NULL,a=NULL){
#   mat=as.matrix(mat)
#   mdim=dim(mat)
#   #   err=1
#   
#   if(is.null(mat.new)){
#     #     err=0
#     mat.new=mat
#   }
#   
#   if(is.null(a))  a=rep(1,mdim[2])
#   if(length(a)<mdim[2]) {
#     a1=rep(1,mdim[2])
#     a1[1:length(a)]=a
#     a=a1;rm(a1)
#     warning('number of "a" is less than the number of columns, use 1 as the missing "a"')
#   }
#   if(length(a)>mdim[2]) {
#     a=a[1:mdim[2]]
#     warning('number of "a" is more than the number of columns, omit the extra "a"')
#   }
#   
#   aa=matrix(rep(a,mdim[1]),ncol=mdim[2],byrow=T)
#   out=(aa*mat)%*%t(mat.new)
#   return(out)
# }

# xixj_sta <- function(mat,mat.new=NULL,w=NULL,power=NULL){
#   mat=as.matrix(mat)
#   if(is.null(mat.new)) mat.new=mat
#   mdim=dim(mat);mdim.new=dim(mat.new)
#   cov.=matrix(sapply(1:mdim[1],function(i) matrix(rep(mat[i,],mdim.new[1]),nrow=mdim.new[1],byrow=T)-mat.new),ncol=mdim[1])
#   if(is.null(power)) power=1
#   cov.=((cov.)^2)^power;
#   if(is.null(w)) {
#     w=rep(1,mdim[2])
#     warning('missing "weight", use 1 instead')
#   }
#   if(length(w)==1&mdim[2]>1){
#     w=rep(w,mdim[2])
#     warning('only one "weight" found, applied to all columns')
#   }
# 
#   if(length(w)>1&length(w)<mdim[2]){
#     w1=rep(1,mdim[2])
#     w1[1:length(w)]=w
#     w=w1;rm(w1)
#     warning('number of "weight" is less than the number of columns, use 1 as the missing "weight"')
#   }
#   if(length(w)>mdim[2]){
#     w=w[1:mdim[2]]
#     warning('number of "weight" is more than the number of columns, omit the extra "weight"')
#   }
#   
#   wmat=matrix(rep(w,each=dim(cov.)[1]*dim(cov.)[2]/mdim[2]),ncol=dim(cov.)[2],byrow=T)
# 
#   cov.=wmat*cov.
# 
#   cov..=matrix(0,ncol=mdim[1],nrow=mdim.new[1])
#   if(mdim[2]>1){
#     for(i in 1:(mdim[2]-1)){
#       cov..=cov..+cov.[1:mdim.new[1],];cov.=cov.[-(1:mdim.new[1]),]}
#     cov.=cov..+cov.
#   }
#   return(cov.)  
# }

# Dpow.ex <- function(vec,Data,hyper,Q=NULL,gamma){
#   DQ=cov.pow.ex(hyper,Data,gamma=gamma);
#   DQ=-0.5*DQ*xixj_sta(as.matrix(vec),w=exp(hyper$pow.ex.w[which(apply(Data,2,mean)==mean(vec) & apply(Data,2,max)==max(vec) & apply(Data,2,min)==  min(vec))]),power=gamma)
#   return(DQ)
# }







# Dloglik.linear.a <- function(hyper,data,AlphaQ){
#   Dlinear.aj=apply(data,2,function(i) sum(AlphaQ*exp(hyper$linear.a[which(data[1,]==i[1])])*DistMatLinear_sq(X=as.matrix(i),A=as.matrix(1))) )
#   return(Dlinear.aj)
# }

Dloglik.linear.a <- function(hyper,data,AlphaQ){
  Dlinear.aj=sapply(1:ncol(data),function(i){
    Xi <- as.matrix(data[,i])
    A1 <- as.matrix(1)
    # res <- sum(AlphaQ*exp(hyper$linear.a[which(data[1,]==i[1])])*DistMatLinear_sq(X=mat_i,A=mat_1))
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


# Dloglik.pow.ex.w <- function(hyper,data,gamma=1,AlphaQ){
#   Dpow.ex.wj=apply(data,2,function(i) sum(AlphaQ*Dpow.ex(as.matrix(i),data,hyper,gamma=gamma)) )
#   return(Dpow.ex.wj)
# }
Dloglik.pow.ex.w <- function(hyper,data,AlphaQ,gamma){
  Dpow.ex.wj=sapply(1:ncol(data),function(i){
    Xi <- as.matrix(data[,i])
    A1 <- as.matrix(1)
    DPsi <- -cov.pow.ex(hyper=hyper, Data=data, gamma=gamma)*DistMat_sq(X=Xi, A=A1, power=gamma)*exp(hyper$pow.ex.w[i])
    res <- 0.5*sum(diag(AlphaQ%*%DPsi))
    # res <- sum(AlphaQ*DPsi)
    return(res)
  })
  return(Dpow.ex.wj)
}

Dloglik.pow.ex.v <- function(hyper,data,AlphaQ,gamma){
  DDpow.ex.v=cov.pow.ex(hyper,data, gamma=gamma)
  Dpow.ex.v=0.5*sum(diag(AlphaQ%*%DDpow.ex.v))
  # Dpow.ex.v=sum(AlphaQ*DDpow.ex.v)
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
      # DPsi <- -3*cc*dist_i*dist_all*exp(-sqrt(3)*dist_all)
      DPsi <- -1.5*cc*exp(hyper$matern.w[i])*dist_i*exp(-sqrt(3)*dist_D)
      # diag(DPsi) <- 0
    }
    if(nu==5/2){
      # expTerm52 <- exp(-sqrt(5)*DistMat_sq(X=data, A=A, power=1)) # CHANGE power to 2
      # DPsi <- cc*((sqrt(5)*dist_i + (10/3)*dist_all*dist_i)*expTerm52 + 
      #   (1+sqrt(5)*dist_all+(5/3)*(dist_all^2))*expTerm52*(-sqrt(5)*dist_i))
      # DPsi <- cc*dist_i*(0.5*sqrt(5)*dist_C^(-0.5) + 
      #                      exp(-sqrt(5)*dist_D)*(
      #                        (5/3) - (sqrt(5)/2)*(dist_C^(-0.5) + sqrt(5) + (5/3)*dist_D)
      #                      ))
      # DPsi <- cc*dist_i*(0.5*sqrt(5)*dist_D^(-1) + 
      #                      exp(-sqrt(5)*dist_D)*(
      #                        (5/3) - (sqrt(5)/2)*(dist_D^(-1) + sqrt(5) + (5/3)*dist_D)
      #                      ))
      DPsi <- cc*exp(hyper$matern.w[i])*dist_i*exp(-sqrt(5)*dist_D)*(5/6)*(-1-sqrt(5)*dist_D)
      # diag(DPsi) <- 0
    }
    res <- 0.5*sum(diag(AlphaQ%*%DPsi))
    # res <- sum(AlphaQ*DPsi)
    return(res)
  })
  return(Dmatern.wj)
}


#######################

Dloglik.rat.qu.w <- function(hyper,data,AlphaQ){

  dimData <- ncol(data)
  d1list <- vector("list", dimData)
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
  }

  Drat.qu.wj <- sapply(1:dimData, function(i){sum(0.5*AlphaQ*d1list[[i]])} )
  # Drat.qu.wj=apply(data,2,function(i) sum(AlphaQ*Drat.qu(i,data,hyper)) )
  return(Drat.qu.wj)
}
# Dloglik.rat.qu.w <- function(hyper,data,AlphaQ){ 
#   Drat.qu.wj=apply(data,2,function(i) sum(AlphaQ*Drat.qu(i,data,hyper)) )
#   return(Drat.qu.wj)
# }


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
# Dloglik.rat.qu.a <- function(hyper,data,AlphaQ){
#   DDrat.qu.a=cov.rat.qu(hyper,data)
#   DDrat.qu.a=log(DDrat.qu.a)%*%DDrat.qu.a
#   Drat.qu.a=sum(AlphaQ*DDrat.qu.a)
#   return(Drat.qu.a)
# }

Dloglik.rat.qu.v <- function(hyper,data,AlphaQ){
  DDrat.qu.v <- cov.rat.qu(hyper,data)
  Drat.qu.v <- 0.5*sum(diag(AlphaQ%*%DDrat.qu.v))
  return(Drat.qu.v)
}



Dloglik.vv <- function(hyper,Alpha,invQ){
  # Dfvv=-sum(diag(invQ))*exp(hyper$vv) + t(Alpha)%*%Alpha*exp(hyper$vv)
  Dfvv=0.5*sum(diag(Alpha%*%t(Alpha) - invQ))*exp(hyper$vv)
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
  
  # dist_all <- DistMat_sq(X=data, A=A, power=1)
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
      # diag(d1list[[i]]) <- 0
      d2list[[i]] <- (-5/6)*cc*dist_i*exp(hyper$matern.w[i])*exp(-sqrt(5)*dist_D)*(
        1+sqrt(5)*dist_D - 2.5*exp(hyper$matern.w[i])*dist_i)
      # diag(d2list[[i]]) <- 0
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



#  #  #
#
D2rat.qu.a <- function(hyper,data,inv.Q,Alpha.Q){
  
  covmatrix <- cov.rat.qu(hyper,data)
  
  #
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
  #
  
  d1 <- log_term*covmatrix*(-hyper$rat.qu.a)
  d2 <- (d1 + covmatrix)*(-hyper$rat.qu.a)*log_term

  D2rat.qu.a <- D2(d1,d2,inv.Q,Alpha.Q)
  
  return(D2rat.qu.a)
}
#
#  #  #
# D2rat.qu.a <- function(hyper,data,inv.Q,Alpha.Q){
#   Q=cov.rat.qu(hyper,data)
#   d1=log(Q)%*%Q
#   d2=log(Q)%*%(log(Q)%*%Q-Q*exp(hyper$rat.qu.a))/exp(hyper$rat.qu.a)
#   D2rat.qu.a=D2(d1,d2,inv.Q,Alpha.Q)
#   return(D2rat.qu.a)
# }


D2rat.qu.v <- function(hyper,data,inv.Q,Alpha.Q){
  covmatrix <- cov.rat.qu(hyper,data)
  D2rat.qu.v <- D2(covmatrix,covmatrix,inv.Q,Alpha.Q)  
  return(D2rat.qu.v)
}






D2vv <- function(hyper,data,inv.Q,Alpha.Q){
  D2fvv=D2(diag(exp(hyper$vv),dim(data)[1]),diag(exp(hyper$vv),dim(data)[1]),inv.Q,Alpha.Q)
  return(D2fvv)
}


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

##################### plot ##########################
plot.gpr <- function(x,...,fitted=F,col.no=1, ylim=NULL){
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
  }
  else{
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
  }
  else{
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
  
  upper=mu+1.96*(sd);
  lower=mu-1.96*(sd);
  if(is.null(ylim)){
    range(upper,lower,Y)
  }
  plot(-100,-100,col=0,xlim=range(X[,col.no],x[,col.no]),ylim=ylim,main=type, xlab="input ",ylab="response",...)
  #
  polygon(c(x[,col.no], rev(x[,col.no])), c(upper, rev(lower)),col = rgb(127,127,127,120, maxColorValue = 255), border = NA)
  #
  points(X[,col.no],Y,pch=pchType,col=2,cex=PcexNo)
  # lines(X[,1],Y)
  lines(x[,col.no],mu,col=4,lwd=LcexNo)  
}






# library("Rcpp")
# library("inline")
# source("DistMat_sq.R")
# source("DistMat.R")
# source("DistMatLinear_sq.R")
# source("DistMatLinear.R")




#hyper is a list of hyper-parameters
#Data should have the form that, each column is a variable
#response is the given response vector
#Cov is a function contains all the covariance matrix, defult is:
###cov.linear(hyper,Data)+cov.pow.ex(hyper,Data)
gp.loglikelihood2NEW <- function(hyper.p,Data, response,Cov,gamma,nu){
  
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
  
  # response <- as.matrix(response)
  # G <- chol(Q)
  # logdetQ <- 2*sum(log(diag(G)))
  # tresp_invQ_resp <- t(response)%*%chol2inv(G)%*%response
  # fX <- 0.5*logdetQ + 0.5*tresp_invQ_resp + 0.5*nrow(Data)*log(2*pi)
  
  
  #  #   #   #   #   
  n <- nrow(response)
  nrep <- ncol(response)
  
  G <- chol(Q)
  logdetQ <- 2*sum(log(diag(G)))
  
  # tresp_invQ_resp <- t(response)%*%chol2inv(G)%*%response
  # fX <- 0
  # for(i in 1:nrep){
  #   tresp_invQ_resp <- t(response[,i])%*%chol2inv(G)%*%response[,i]
  #   fX_i <- 0.5*logdetQ + 0.5*tresp_invQ_resp + 0.5*n*log(2*pi)
  #   fX <- fX + fX_i
  # }
  # fX
  
  tresp_invQ_resp <- t(response)%*%chol2inv(G)%*%response
  if(nrep==1){
    fX <- 0.5*logdetQ + 0.5*tresp_invQ_resp + 0.5*n*log(2*pi)
  }else{
    fX <- nrep*0.5*logdetQ + 0.5*sum(diag( tresp_invQ_resp )) + nrep*0.5*n*log(2*pi)
  }
  fX <- as.numeric(fX)
  # fX
  
  #  #   #   #   #   

  return(fX)
}


gp.Dlikelihood2NEW <- function(hyper.p,  Data, response, Cov, gamma, nu){
  
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
  
  # response <- as.matrix(response)
  
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
    return(Dpara)
  })
  
  names(Dfx)=names(hyper.p)
  # Dfx=-0.5*unlist(Dfx)
  Dfx <- - unlist(Dfx) # Minus loglik
  
  DfxList[[irep]] <- Dfx
}
  
Dfx <- Reduce('+', DfxList)
  
  
  return(Dfx)
}






#' Convolved (multivariate) Gaussian process regression
#'
#' @param Data List of two elements: (1) response variable; and (2) input
#'   variables
#' @param m If Subset of Data is to be used, m denotes the subset size and
#'   cannot be larger than the total sample size. Default set to NULL.
#' @param meanModel Model for mean.
#' @param mean Is the mean taken out when analysis? Default to be 0, which
#'   assumes the mean is zero. if assume mean is a constant, mean=1; if assume
#'   mean is a linear trend, mean='t'.
#'
#' @return A list containing: \describe{ 
#' \item{fitted.mean }{Fitted value of training data } 
#' \item{fitted.sd }{Standard deviation of the fitted value of training data} 
#' \item{X}{Original input variables}
#' \item{Y}{Original response}
#' \item{idx}{Original response}
#' 
#' \item{Cov}{Original response}
#' \item{mean}{CHECK}
#' \item{meanModel}{CHECK  Mean function type}
#' \item{meanLinearModel}{CHECK Mean model}
#' }
#' 
#' @references Shi, J. Q., and Choi, T. (2011), ``Gaussian Process Regression
#'  Analysis for Functional Data'', CRC Press.
#' @export
#'
CGPR <- function(Data, m=NULL, meanModel=0, mean=NULL){
  
  N <- length(Data$input)
  X <- as.matrix(unlist(Data$input))
  ns <- sapply(Data$input, length)
  idx <- c(unlist(sapply(1:N, function(i) rep(i, ns[i]))))
  
  response <- Reduce('rbind', Data$response)
  
  Q <- ncol(X)
  nrep <- ncol(response)
  n <- nrow(response)
  
  Y.original <- response
  X.original <- X
  idx.original <- idx
  n.original <- n
  
  if(!is.null(mean)){
    if(!length(mean)==n){
      stop("'mean' defined by the user must have the same length as the response variable.")
    }
    mean <- matrix(rep(mean, nrep), ncol=nrep, byrow=F)
    meanModel <- 'userDefined'
  }
  
  if(meanModel==0) {
    mean <- 0
    mean <- matrix(mean, nrow=n, ncol=nrep, byrow=F)
  }

  if(meanModel==1) {
    responseNew <- NULL
    mean <- NULL
    for(j in 1:N){
      
      mean_j <- mean(response[idx==j,])
      nj <- nrow(response[idx==j,,drop=F])
      mean_j <- matrix( rep(mean_j, nj*nrep), nrow=nj, byrow=F)
      response_j <- response[idx==j,,drop=F] - mean_j
      responseNew <- rbind(responseNew, response_j)
      mean <- rbind(mean, mean_j)
    }
    response <- responseNew
    
  }

  if(meanModel=='t') {
    meanLinearModel <- list()  
    responseNew <- NULL
    mean <- NULL
    for(j in 1:N){
      trend <- data.frame(yyy=c(response[idx==j,]), xxx=rep(c(X[idx==j,]), nrep))
      meanLinearModel_j <- lm(yyy~xxx, data=trend)
      meanLinearModel[[j]] <- meanLinearModel_j
      response_j <- matrix(resid(meanLinearModel_j), nrow=nrow(response[idx==j,,drop=F]), byrow=F)
      mean_j <- matrix(fitted(meanLinearModel_j), nrow=nrow(response[idx==j,,drop=F]), byrow=F)
      
      responseNew <- rbind(responseNew, response_j)
      mean <- rbind(mean, mean_j)
    }
    response <- responseNew
  }else{
    meanLinearModel <- NULL
  }
  
  if(meanModel=='avg') {
    if(nrep<3){
      stop('Mean function can only be the average across replications when
           there are more than two replications.')
    }
    mean <- apply(response, 1, base::mean)
    mean <- matrix(rep(mean, nrep), ncol=nrep, byrow=F)
    response <- response - mean
  }
  
  mean.original <- mean
  
  idxSubset <- NULL
  if(!is.null(m)){
    if(m>n){stop("m cannot be bigger than n.")}
    idxSubset <- sort(sample(x=1:n, size=m, replace=F))
    response <- response[idxSubset,,drop=F]
    X <- X[idxSubset,,drop=F]
    idx <- idx[idxSubset]
    # if(!is.null(mean)){
      mean <- mean[idxSubset,,drop=F]
    # }
    n <- nrow(response)
  }
    
  # c(va0s, va1s, Apars, sig)
  lowerlimits <- c(rep(-50, N), rep(log(1e-3), N+2*N*Q), log(1e-8))
  upperlimits <- c(rep(50, N),  rep(log(3000), N+2*N*Q), log(3))
  
  # Select the starting values
  nCand <- 100
  candidates <- matrix(0, nCand, 2*N + 2*N*Q + 1)
  for(iCand in 1:nCand){
    nu0s <- sample(c(-1,1), size = N, replace = T)*rep(runif(1, min=-2,max=2), N)
    nu1s <- log(abs(nu0s)*5)
    a0s <- runif(N, min=log(1.1), max=log(30))
    a1s <- runif(N, min=log(1.1), max=log(30))
    sigm <- runif(1, min=log(1e-2), max=log(0.1))
    candidates[iCand,] <- c(nu0s, nu1s, rep(c(a0s,a1s), Q), sigm)
  }
  
  resCand <- apply(candidates, 1, function(x) LogLikCGP(x, response, X, idx))
  hp_init_log <- candidates[which.min(resCand),]
  
  res <- nlminb(start=hp_init_log, objective=LogLikCGP, gradient = NULL, hessian = NULL,
                       control = list(eval.max=1000, iter.max=1000,
                                      rel.tol=1e-8, x.tol=1e-8, xf.tol=1e-8),
                       lower = lowerlimits, upper = upperlimits, response=response, X=X, idx=idx)
  
  hp_opt <- res$par
  
  K <- CGPCovMat(Data=Data, hp=hp_opt)
  invK <- chol2inv(chol(K))
  
  varEpsilon <- exp(hp_opt[length(hp_opt)])^2
    
  fitted <- (K-diag(varEpsilon, n.original))%*%invK%*%Y.original + mean.original
  fitted.var <- varEpsilon*rowSums((K-diag(varEpsilon, n.original))*t(invK))

  result <- list('hyper'=hp_opt,
              'fitted.mean'=fitted,
              'fitted.sd'=sqrt(fitted.var),
              'X'=X.original, 'Y'=Y.original, 'idx'=idx.original, 'Cov'=K, 'mean'=mean.original[,1], 
              'meanModel'=meanModel, 'meanLinearModel'=meanLinearModel)
  class(result)='mgpr'
  return(result=result)
}



#' Prediction of Convolved Gaussian process
#'
#' @inheritParams CGPR
#' @param train List resulting from training which is an 'mgpr' object.
#' @param Data.train List of training data
#' @param Data.new List of test input data
#' @param noiseFreePred Logical. If TRUE, predictions will be noise-free.
#'
#' @export
#'
#' @return A list containing  \describe{ 
#' \item{pred.mean}{Mean of predictions}
#' \item{pred.sd}{Standard deviation of predictions}
#' \item{noiseFreePred}{Logical. If TRUE, predictions are noise-free.}
#' }
#' 
CGPprediction <- function(train=NULL, 
                       Data.train=NULL,
                       Data.new=NULL,
                       noiseFreePred=F, 
                       meanModel=NULL, mean=0){
  
  if(class(train)=='mgpr'){
    hyper=train$hyper
    X=train$X
    Y=train$Y
    idx=train$idx
    Cov=train$Cov
    mean=train$mean
    meanModel=train$meanModel
    meanLinearModel=train$meanLinearModel
  }
  
  
  N <- length(Data.train$input)
  
  X.new <- as.matrix(unlist(Data.new$input))
  ns.new <- sapply(Data.new$input, length)
  idx.new <- c(unlist(sapply(1:N, function(i) rep(i, ns.new[i]))))
  
  X <- as.matrix(unlist(Data.train$input))
  ns <- sapply(Data.train$input, length)
  nsTest <- sapply(Data.new$input, length)
  idx <- c(unlist(sapply(1:N, function(i) rep(i, ns[i]))))
  Y <- Reduce('rbind', Data.train$response)
  
  nrep <- ncol(Y)
  
  if(meanModel==0){
    meanList <- list()
    for(j in 1:N){
      meanList[[j]] <- rep(0, ns[j])
    }
    meanY <- do.call(cbind, replicate(nrep, unlist(meanList), simplify=FALSE))
  }
  if(meanModel=='t'){
    meanList <- list()
    for(j in 1:N){
      newtrend <- data.frame(xxx=Data.train$input[[j]])
      meanList[[j]] <- predict(meanLinearModel[[j]], newdata=newtrend)
    }
    meanY <- do.call(cbind, replicate(nrep, unlist(meanList), simplify=FALSE))
  }
  Y <- Y - meanY
  
  
  
  hp <- TransfToNatScaleCGP(hyper, N)
  
  Q <- ncol(X)
  va0s <- hp[1:N]
  va1s <- hp[(N+1):(2*N)]
  AparsMat <- matrix(hp[seq(2*N+1, by=1, length.out=2*N*Q)], ncol=2*Q, byrow=T)
  A0s <- A1s <- list()
  
  if(Q==1){
    for(j in 1:N){
      A0s[[j]] <- as.matrix(AparsMat[j,1:Q])
      A1s[[j]] <- as.matrix(AparsMat[j,(Q+1):(2*Q)])
    }
  }else{
    for(j in 1:N){
      A0s[[j]] <- diag(AparsMat[j,1:Q])
      A1s[[j]] <- diag(AparsMat[j,(Q+1):(2*Q)])
    }
  }
  
  sig <- hp[length(hp)]
  
  Psi <- KCGP(X=X, idx=idx, va0s=va0s, va1s=va1s, A0s=A0s, A1s=A1s, sig=sig)
  Knm <- KCGPnm(X=X, Xp = X.new, idx=idx, idx_new = idx.new, va0s=va0s, va1s=va1s, A0s=A0s, A1s=A1s, sig=0)
  
  Kstar <- KCGP(X=X.new, idx=idx.new, va0s=va0s, va1s=va1s, A0s=A0s, A1s=A1s, sig=sig)

  invPsi <- chol2inv(chol(Psi))
  QR <- invPsi%*%Y
  
  if(meanModel==0){
    meanList <- list()
    for(j in 1:N){
      meanList[[j]] <- rep(0, nsTest[j])
    }
    mean <- do.call(cbind, replicate(nrep, unlist(meanList), simplify=FALSE))
  }
  if(meanModel=='t'){
    meanList <- list()
    for(j in 1:N){
      newtrend <- data.frame(xxx=Data.new$input[[j]])
      meanList[[j]] <- predict(meanLinearModel[[j]], newdata=newtrend)
    }
    mean <- do.call(cbind, replicate(nrep, unlist(meanList), simplify=FALSE))
  }
  
  mu <- t(Knm)%*%QR + mean
  
  if(noiseFreePred){
    sigma2 <- diag(Kstar)-diag(t(Knm)%*%invPsi%*%Knm) - sig^2
  }else{
    sigma2 <- diag(Kstar)-diag(t(Knm)%*%invPsi%*%Knm)
  }
  pred.sd. <- sqrt(sigma2)
    

  
  pred.mean <- list()
  pred.sd <- list()
  for(j in 1:N){
    pred.mean[[j]] <- mu[idx.new==j,,drop=F]
    pred.sd[[j]] <- pred.sd.[idx.new==j]
  }
  
  
  result=c(list('pred.mean'=pred.mean,
                'pred.sd'=pred.sd,
                'noiseFreePred'=noiseFreePred))
  class(result)='mgpr'
  return(result)
  
}



#' Calculation of a Convolved Gaussian processes covariance matrix given a
#' vector of hyperparameters
#'
#' @inheritParams CGPR
#' @param hp Vector of hyperparameters
#' @references Shi, J. Q., and Choi, T. (2011), ``Gaussian Process Regression
#'  Analysis for Functional Data'', CRC Press.
#'  
#' @return Covariance matrix
#' @export
CGPCovMat <- function(Data, hp){

  N <- length(Data$input)
  X <- as.matrix(unlist(Data$input))
  ns <- sapply(Data$input, length)
  
  idx <- c(unlist(sapply(1:N, function(i) rep(i, ns[i]))))
  hp <- TransfToNatScaleCGP(hp, N)
  
  Q <- ncol(X)
  va0s <- hp[1:N]
  va1s <- hp[(N+1):(2*N)]
  AparsMat <- matrix(hp[seq(2*N+1, by=1, length.out=2*N*Q)], ncol=2*Q, byrow=T)
  A0s <- A1s <- list()
  
  if(Q==1){
    for(j in 1:N){
      A0s[[j]] <- as.matrix(AparsMat[j,1:Q])
      A1s[[j]] <- as.matrix(AparsMat[j,(Q+1):(2*Q)])
    }
  }else{
    for(j in 1:N){
      A0s[[j]] <- diag(AparsMat[j,1:Q])
      A1s[[j]] <- diag(AparsMat[j,(Q+1):(2*Q)])
    }
  }

  sig <- hp[length(hp)]
  K <- KCGP(X=X, idx=idx, va0s=va0s, va1s=va1s, A0s=A0s, A1s=A1s, sig=sig)
  return(K)
}


TransfToNatScaleCGP <- function(hp, N){
  hp[(N+1):length(hp)] <- exp(hp[(N+1):length(hp)])
  return(hp)
}



LogLikCGP <- function(hp, response, X, idx){
  
  Q <- ncol(X)
  N <- length(unique(idx))
  hp <- TransfToNatScaleCGP(hp, N)
  
  va0s <- hp[1:N]
  va1s <- hp[(N+1):(2*N)]
  AparsMat <- matrix(hp[seq(2*N+1, by=1, length.out=2*N*Q)], ncol=2*Q, byrow=T)
  A0s <- A1s <- list()
  if(Q==1){
    for(j in 1:N){
      A0s[[j]] <- as.matrix(AparsMat[j,1:Q])
      A1s[[j]] <- as.matrix(AparsMat[j,(Q+1):(2*Q)])
    }
  }else{
    for(j in 1:N){
      A0s[[j]] <- diag(AparsMat[j,1:Q])
      A1s[[j]] <- diag(AparsMat[j,(Q+1):(2*Q)])
    }
  }
  sig <- hp[length(hp)]
  
  K <- KCGP(X=X, idx=idx, va0s=va0s, va1s=va1s, A0s=A0s, A1s=A1s, sig=sig)

  G <- chol(K)
  logdetK <- 2*sum(log(diag(G)))
  yt.invK.y <- t(response)%*%chol2inv(G)%*%response

  n <- nrow(response)
  nrep <- ncol(response)
  
  if(nrep==1){
    fX <- 0.5*logdetK + 0.5*yt.invK.y + 0.5*n*log(2*pi)
  }else{
    fX <- nrep*0.5*logdetK + 0.5*sum(diag( yt.invK.y )) + nrep*0.5*n*log(2*pi)
  }
  fX <- as.numeric(fX)
  
  return(fX)
}




#' Plot Convolved (multivariate) Gaussian Process regression
#'
#' Plot Gaussian Process for a given an object of class 'gpr'.
#'
#' @param train The 'mgpr' object 
#' @param Data.train List of training data
#' @param Data.new List of test data
#' @param i Which realisation should be plotted.
#' @param ylim Range of y-axis
#' @param mfrow Graphical parameter
#' @param cex  Graphical parameter
#' @param cex.lab  Graphical parameter
#' @param cex.axis  Graphical parameter
#'
#' @return A plot showing predictions of each element of the multivariate process.
#' @export
plotCGPprediction <- function(train, Data.train, Data.new, i, ylim=NULL, mfrow=NULL,
                               cex=1, cex.lab=1, cex.axis=1){
  
  op <- par(mar=c(4.5,5.1,0.2,0.8), 
            oma=c(0,0,0,0),
            cex.lab=2, cex.axis=2, cex.main=2)
  
  predCGP <- CGPprediction(train=train, 
                        Data.train=Data.train,
                        Data.new=Data.new)
  
  N <- length(predCGP$pred.mean)
  
  if(is.null(mfrow)){
    if(N<4){
      par(mfrow=c(1,N))
    }
  }else{
    par(mfrow=mfrow)
  }
  
  for(variable in 1:N){
    
    predMean <- predCGP$pred.mean[[variable]][,i]
    upper <- predMean+1.96*predCGP$pred.sd[[variable]]
    lower <- predMean-1.96*predCGP$pred.sd[[variable]]
    
    if(is.null(ylim)){
      ylim_i <- range(c(lower, upper)) #+ c(-0.5,0.5)
    }else{
      ylim_i <- ylim[[variable]]
    }
    
    xlim_i <- range(Data.train$input[[variable]], Data.new$input[[variable]])
    plot(Data.train$input[[variable]], Data.train$response[[variable]][,i], type="p",
         xlab="t", ylab=bquote(X[.(variable)]), ylim=ylim_i, xlim=xlim_i, pch=19, 
         cex=cex, cex.axis=cex.axis, cex.lab=cex.lab)
    lines(Data.new$input[[variable]], predMean, col="blue", lwd=2)
    
    polygon(x=c(Data.new$input[[variable]], rev(Data.new$input[[variable]])), 
            y=c(upper, rev(lower)),
            col = rgb(127,127,127,120, maxColorValue = 255), border = NA)
  }
  
  par(mfrow=c(1,1))
  par(op)
}





#' Plot of auto- or cross-covariance function of a Convolved (multivariate)
#' Gaussian process
#'
#' @param type Logical. It can be either 'Cov' (for covariance function) or
#'   'Cor' (for corresponding correlation function)
#' @param output Integer identifying one element of the multivariate process
#' @param outputp Integer identifying one element of the multivariate process.
#'   If 'output' and 'outputp' are the same, the auto-covariance function will
#'   be plotted. Otherwise, the cross-covariance function will be plotted.
#' @param Data List of Data
#' @param hp Vector of hyperparameters
#' @param ylim Graphical parameter
#' @param xlim Graphical parameter
#'
#' @return A plot
#' @export
#' 
plotCGPCovFun <- function(type="Cov", output, outputp, Data, hp, ylim=NULL, xlim=NULL){
  
  op <- par(mar=c(4.5,5.1,0.2,0.8), 
            oma=c(0,0,0,0),
            cex.lab=2, cex.axis=2, cex.main=2)
  
  Psi <- CGPCovMat(Data=Data, hp=hp)
  
  if(type=="Cor"){
    Psi <- cov2cor(Psi)
  }
  
  toPlot <- Psi[idx==output, idx==outputp][,1]
  tp0 <- Data$input[[outputp]][1]
  
  if(!is.null(ylim)){
    ylim <- range(toPlot)
  }
  plot(Data$input[[output]] - tp0, toPlot, ylim=ylim, xlim=xlim, type="l",
       xlab=bquote("t-("*.(tp0)*")"),
       ylab=bquote(.(type)*"["~X[.(output)]*"(t),"~
                     X[.(outputp)]*"("*.(tp0)*")]"))

  par(op)
}




#' Multivariate Gaussian process regression (MGPR) model
#'
#' Multivariate Gaussian process regression where each of the N
#' outputs is unidimensional. The multivariate output is allowed to have
#' multiple independent realisations.
#'
#' @param Data List of two elements: 'input' and 'response'. The element 'input'
#'   is a list of N vectors, where each vector represents the input covariate
#'   values for a particular output. The element 'response' is the corresponding
#'   list of N matrices (if there are multiple realisations) or vectors (for a
#'   single realisation) representing the response variables.
#' @param m If Subset of Data is to be used in the estimation, m denotes the
#'   subset size. It cannot be larger than the total sample size. Default to
#'   NULL (Subsetting is not used).
#' @param meanModel Type of mean function applied to all outputs. It can be
#'   \describe{ \item{0}{Zero mean function for each output.} \item{1}{Constant
#'   mean function to be estimated for each output.} \item{'t'}{Linear model for
#'   the mean function of each output.} \item{'avg'}{The average across
#'   replications is used as the mean function of each output. This can only be
#'   used if there are more than two realisations observed at the same input
#'   values.} } Default to 0. If argument 'mu' is specified, then 'meanModel'
#'   will be set to 'userDefined'.
#' @param mu Vector of concatenated mean function values defined by the user.
#'   Default to NULL.
#'
#' @return A list containing: \describe{ \item{fitted.mean }{Fitted values for
#'   the training data } \item{fitted.sd }{Standard deviation of the fitted
#'   values for training data} \item{N}{Number of response variables}
#'   \item{X}{Original input variables}
#'   \item{Y}{Original response} \item{idx}{Index vector identifying to which 
#'   output the elements of concatenated vectors correspond to.} 
#'   \item{Cov}{Covariance matrix} \item{mean}{Concatenated mean function } 
#'   \item{meanModel}{Mean model used for each output} 
#'   \item{constantMeanEstimates}{Constant mean estimated for each output if 
#'   the constant mean model is used for the mean functions.}
#'   \item{meanLinearModel}{'lm' object for each output if the linear regression 
#'   model is used for the mean functions.} }
#'
#' @references Shi, J. Q., and Choi, T. (2011), ``Gaussian Process Regression
#'   Analysis for Functional Data'', CRC Press.
#' @export
#' @examples
#' ## See examples in vignette:
#' # vignette("mgpr", package = "GPFDA")
mgpr <- function(Data, m=NULL, meanModel=0, mu=NULL){
  
  N <- length(Data$input)
  X <- as.matrix(unlist(Data$input))
  ns <- sapply(Data$input, length)
  idx <- c(unlist(sapply(1:N, function(i) rep(i, ns[i]))))
  
  response <- Reduce('rbind', Data$response)
  
  Q <- 1
  nrep <- ncol(response)
  n <- nrow(response)
  
  Y.original <- response
  X.original <- X
  idx.original <- idx
  n.original <- n
  
  if(!is.null(mu)){
    if(!length(mu)==n){
      stop("'mu' defined by the user must have the same length as the response 
           variable.")
    }
    mu <- matrix(rep(mu, nrep), ncol=nrep, byrow=F)
    meanModel <- 'userDefined'
  }
  
  if(meanModel==0) {
    mu <- 0
    mu <- matrix(mu, nrow=n, ncol=nrep, byrow=F)
  }

  if(meanModel==1) {
    constantMeanEstimates <- list()
    responseNew <- NULL
    mu <- NULL
    for(j in 1:N){
      
      mean_j <- mean(response[idx==j,])
      constantMeanEstimates[[j]] <- mean_j
      nj <- nrow(response[idx==j,,drop=F])
      mean_j <- matrix( rep(mean_j, nj*nrep), nrow=nj, byrow=F)
      response_j <- response[idx==j,,drop=F] - mean_j
      responseNew <- rbind(responseNew, response_j)
      mu <- rbind(mu, mean_j)
    }
    response <- responseNew
    
  }else{
    constantMeanEstimates <- NULL
  }

  if(meanModel=='t') {
    meanLinearModel <- list()  
    responseNew <- NULL
    mu <- NULL
    for(j in 1:N){
      trend <- data.frame(yyy=c(response[idx==j,]), xxx=rep(c(X[idx==j,]), 
                                                            nrep))
      meanLinearModel_j <- lm(yyy~xxx, data=trend)
      meanLinearModel[[j]] <- meanLinearModel_j
      response_j <- matrix(resid(meanLinearModel_j), 
                           nrow=nrow(response[idx==j,,drop=F]), byrow=F)
      mean_j <- matrix(fitted(meanLinearModel_j), 
                       nrow=nrow(response[idx==j,,drop=F]), byrow=F)
      
      responseNew <- rbind(responseNew, response_j)
      mu <- rbind(mu, mean_j)
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
    mu <- apply(response, 1, mean)
    mu <- matrix(rep(mu, nrep), ncol=nrep, byrow=F)
    response <- response - mu
  }
  
  mean.original <- mu
  
  idxSubset <- NULL
  if(!is.null(m)){
    if(m>n){stop("m cannot be bigger than n.")}
    idxSubset <- sort(sample(x=1:n, size=m, replace=F))
    response <- response[idxSubset,,drop=F]
    X <- X[idxSubset,,drop=F]
    idx <- idx[idxSubset]
    # if(!is.null(mu)){
      mu <- mu[idxSubset,,drop=F]
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
    nu0s <- sample(c(-1,1), size=N, replace=T)*rep(runif(1, min=-2,max=2), N)
    nu1s <- log(abs(nu0s)*5)
    a0s <- runif(N, min=log(1.1), max=log(30))
    a1s <- runif(N, min=log(1.1), max=log(30))
    sigm <- runif(1, min=log(1e-2), max=log(0.1))
    candidates[iCand,] <- c(nu0s, nu1s, rep(c(a0s,a1s), Q), sigm)
  }
  
  resCand <- apply(candidates, 1, function(x) LogLikCGP(x, response, X, idx))
  hp_init_log <- candidates[which.min(resCand),]
  
  res <- nlminb(start=hp_init_log, objective=LogLikCGP, gradient=NULL, 
                hessian=NULL,
                control=list(eval.max=1000, iter.max=1000,
                               rel.tol=1e-8, x.tol=1e-8, xf.tol=1e-8),
                lower=lowerlimits, upper=upperlimits, response=response, X=X, 
                idx=idx)
  
  hp_opt <- res$par

  K <- mgpCovMat(Data=Data, hp=hp_opt)
  invK <- chol2inv(chol(K))
  
  varEpsilon <- exp(hp_opt[length(hp_opt)])^2
    
  fitted <- (K-diag(varEpsilon, n.original))%*%invK%*%Y.original + mean.original
  fitted.var <- varEpsilon*rowSums((K-diag(varEpsilon, n.original))*t(invK))

  nu0names <- paste0("nu0_", 1:N)
  nu1names <- paste0("nu1_", 1:N)
  a0names <- paste0("a0_", 1:N)
  a1names <- paste0("a1_", 1:N)
  signame <- "sigma"
  names(hp_opt) <- c(nu0names, nu1names, a0names, a1names, signame)
  
  result <- list('hyper'=hp_opt,
              'fitted.mean'=fitted,
              'fitted.sd'=sqrt(fitted.var),
              'N'=N,
              'X'=X.original, 'Y'=Y.original, 'idx'=idx.original, 'Cov'=K, 
              'mu'=mean.original[,1], 
              'meanModel'=meanModel, 'constantMeanEstimates'=constantMeanEstimates,
              'meanLinearModel'=meanLinearModel)
  class(result)='mgpr'
  
  return(result)
}


#' Print method for 'mgpr' objects
#'
#' @param x An 'mgpr' object.
#' @param ... Not used here.
#'
#' @return Returns a summary of an object of class 'mgpr'.
#' @export
#'
#' @seealso \link[GPFDA]{mgpr}
print.mgpr <- function(x, ...) {
  
  if(class(x) != "mgpr"){
    stop("'x' must be of class 'mgpr'")
  }
  
  cat(paste0("MGPR model specifications:", "\n"))
  cat("\n")
  
  
  cat(paste0("Number of response variables: ", x$N, "\n"))
  cat("\n")
  
  if(x$meanModel==0){
    cat(paste0("Zero mean function of each response \n"))
  }else if(x$meanModel==1){
    cat(paste0("Constant mean function of each response \n"))
  }else if(x$meanModel=="t"){
    cat(paste0("Linear model for the mean function of each response \n"))
  }else if(x$meanModel=="avg"){
    cat(paste0("Average across replications is used as the mean function of each response \n"))
  }
  
}




#' Summary method for 'mgpr' objects
#'
#' @param object An 'mgpr' object.
#' @param ... Not used here.
#' 
#' @return Some fitting results of the MGPR model.
#' @export
#'
#' @seealso \link[GPFDA]{mgpr}
summary.mgpr <- function(object, ...){
  
  if(class(object) != "mgpr"){
    stop("'object' must be of type mgpr")
  }
  
  N <- object$N
  
  if(object$meanModel==0){
    cat(paste0("Mean function assumed to be zero.", "\n"))
    cat("\n")
  }
  
  if(object$meanModel==1){
    cat(paste0("Estimated constant mean function of each response variable:", "\n"))
    
    for(j in 1:N){
      cat(paste0("For response variable ", j, ":\n"))
      cat(paste0("mean = ", round(object$constantMeanEstimates[[j]], digits = 6), "\n"))  
    }
    
    cat("\n")
  }
  
  if(object$meanModel=="t"){
    cat(paste0("Coefficients of the fitted linear model for the mean function of each response variable:", "\n"))
    cat("\n")
    for(j in 1:N){
      cat(paste0("Linear model for response variable ", j, ":\n"))
      linModel <- object$meanLinearModel[[j]]$coefficients
      linModel <- round(linModel, digits = 6)
      names(linModel) <- c("intercept", paste0("input"))
      print(linModel)
      cat("\n")
    }
    cat("\n")
  }
  
  cat(paste0("Marginal ML estimates of hyperparameters:", "\n"))
  print(round(object$hyper, digits = 6))
}





  

#' Prediction of MGPR model
#'
#' @inheritParams mgpr
#' @param object An 'mgpr' object obtained from 'mgpr' function. 
#'   If NULL, predictions are made based on DataObs informed by the user.
#' @param DataObs List of observed data. Default to NULL. If NULL,
#'   predictions are made based on the trained data 
#'   (included in the object of class 'mgpr') used for learning.
#' @param DataNew List of test input data.
#' @param noiseFreePred Logical. If TRUE, predictions will be noise-free.
#' @param ... Not used here.
#' @export
#'
#' @return A list containing  \describe{ 
#' \item{pred.mean}{Mean of predictions for the test set.}
#' \item{pred.sd}{Standard deviation of predictions for the test set.}
#' \item{noiseFreePred}{Logical. If TRUE, predictions are noise-free.} }
#'   
#' @examples
#' ## See examples in vignette:
#' # vignette("mgpr", package = "GPFDA")
predict.mgpr <- function(object, 
                       DataObs=NULL,
                       DataNew,
                       noiseFreePred=F, 
                       meanModel=NULL, mu=0, ...){
  
  
  if( class(object)!="mgpr"){
    stop("'object' must be either NULL or of type 'mgpr'")
  }
  
  train <- object
  
  hyper <- train$hyper
  X <- train$X
  Y <- train$Y
  N <- train$N
  idx <- train$idx
  Cov <- train$Cov
  mu <- train$mu
  meanModel <- train$meanModel
  meanLinearModel <- train$meanLinearModel

  
  if(!is.null(DataObs)){
    N <- length(DataObs$input)
    X <- as.matrix(unlist(DataObs$input))
    
    if(!is.matrix(DataObs$response[[1]])){
      DataObs$response <- lapply(DataObs$response, as.matrix)
    }
  }

  X.new <- as.matrix(unlist(DataNew$input))
  ns.new <- sapply(DataNew$input, length)
  idx.new <- c(unlist(sapply(1:N, function(i) rep(i, ns.new[i]))))

  ns <- sapply(DataObs$input, length)
  nsTest <- sapply(DataNew$input, length)
  idx <- c(unlist(sapply(1:N, function(i) rep(i, ns[i]))))
  Y <- Reduce('rbind', DataObs$response)
  
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
      newtrend <- data.frame(xxx=DataObs$input[[j]])
      meanList[[j]] <- predict(meanLinearModel[[j]], newdata=newtrend)
    }
    meanY <- do.call(cbind, replicate(nrep, unlist(meanList), simplify=FALSE))
  }
  Y <- Y - meanY

  hp <- TransfToNatScaleCGP(hyper, N)
  
  Q <- 1
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
  Knm <- KCGPnm(X=X, Xp = X.new, idx=idx, idx_new = idx.new, va0s=va0s, 
                va1s=va1s, A0s=A0s, A1s=A1s, sig=0)
  
  Kstar <- KCGP(X=X.new, idx=idx.new, va0s=va0s, va1s=va1s, A0s=A0s, A1s=A1s, 
                sig=sig)

  invPsi <- chol2inv(chol(Psi))
  QR <- invPsi%*%Y
  
  if(meanModel==0){
    meanList <- list()
    for(j in 1:N){
      meanList[[j]] <- rep(0, nsTest[j])
    }
    mu <- do.call(cbind, replicate(nrep, unlist(meanList), simplify=FALSE))
  }
  if(meanModel=='t'){
    meanList <- list()
    for(j in 1:N){
      newtrend <- data.frame(xxx=DataNew$input[[j]])
      meanList[[j]] <- predict(meanLinearModel[[j]], newdata=newtrend)
    }
    mu <- do.call(cbind, replicate(nrep, unlist(meanList), simplify=FALSE))
  }
  
  pred.mu. <- t(Knm)%*%QR + mu
  
  if(noiseFreePred){
    sigma2 <- diag(Kstar)-diag(t(Knm)%*%invPsi%*%Knm) - sig^2
  }else{
    sigma2 <- diag(Kstar)-diag(t(Knm)%*%invPsi%*%Knm)
  }
  pred.sd. <- sqrt(sigma2)

  pred.mean <- list()
  pred.sd <- list()
  for(j in 1:N){
    pred.mean[[j]] <- pred.mu.[idx.new==j,,drop=F]
    pred.sd[[j]] <- pred.sd.[idx.new==j]
  }

  result=c(list('pred.mean'=pred.mean,
                'pred.sd'=pred.sd,
                'noiseFreePred'=noiseFreePred))
  class(result)='mgpr'
  
  return(result)
  
}



#' Calculate a multivariate Gaussian processes covariance matrix given a
#' vector of hyperparameters
#'
#' @inheritParams mgpr
#' @param hp Vector of hyperparameters
#' @references Shi, J. Q., and Choi, T. (2011), ``Gaussian Process Regression
#'  Analysis for Functional Data'', CRC Press.
#'  
#' @return Covariance matrix
#' @export
#' @examples
#' ## See examples in vignette:
#' # vignette("mgpr", package = "GPFDA")
mgpCovMat <- function(Data, hp){

  N <- length(Data$input)
  X <- as.matrix(unlist(Data$input))
  ns <- sapply(Data$input, length)
  
  idx <- c(unlist(sapply(1:N, function(i) rep(i, ns[i]))))
  hp <- TransfToNatScaleCGP(hp, N)
  
  Q <- 1
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
  
  Q <- 1
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







#' Plot predictions of MGPR model
#' 
#' Plot predictons of each of the N elements of the multivariate Gaussian Process 
#' for a given an object of class 'mgpr'.
#'
#' @param x An object of class 'mgpr'.
#' @param DataObs List of observed data.
#' @param DataNew List of test data.
#' @param nCol Number of columns that the grid of plots should have.
#' @param xlab Title for the x-axis.
#' @param ylabs Title for the y-axis.
#' @param ylims Ranges of values for y-axis for each of the N response variables. It must be a matrix with N rows and 2 columns.
#' @param alpha Significance level used for 'fitted' or 'prediction'. Default is 0.05.
#' @param colourData Colour for training data.
#' @param colourPred Colour for predictive mean for the new curves.
#' @param lwd Graphical parameter.
#' @param cex.points Graphical parameter.
#' @param cex.lab Graphical parameter.
#' @param cex.axis Graphical parameter.
#' @param ... Not used here.
#' @import  ggplot2
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange
#' @return A plot showing predictions of each element of the multivariate
#'   process.
#' @export
#' @examples 
#' ## See examples in vignette:
#' # vignette("mgpr", package = "GPFDA")
plot.mgpr <- function(x, DataObs, DataNew, nCol=NULL,
                      xlab='t', ylabs=NULL, ylims=NULL,
                      alpha=0.05, 
                      colourData="black", colourPred="red", lwd=0.5, cex.points=2,
                      cex.lab=10, cex.axis=10, ...){
  
  if(class(x) != "mgpr"){
    stop("'object' must be of type mgpr")
  }

  z <- stats::qnorm(1-alpha/2)
  
  if(!is.matrix(DataObs$response[[1]])){
    DataObs$response <- lapply(DataObs$response, as.matrix)
  }
  
  predCGP <- predict.mgpr(object=x, 
                         DataObs=DataObs,
                         DataNew=DataNew)
  
  N <- length(predCGP$pred.mean)
  
  outList <- vector('list', N)
  
  for(variable in 1:N){
    
    predMean <- predCGP$pred.mean[[variable]]
    upper <- predMean+z*predCGP$pred.sd[[variable]]
    lower <- predMean-z*predCGP$pred.sd[[variable]]
    
    df <- data.frame(t=DataNew$input[[variable]], y=predMean)
    meltdf <- melt(df, id="t")
    
    out <- ggplot(data=meltdf, aes(x=t,y=value,colour=colourPred,group=variable),
                  lwd=lwd, linetype="dashed",
                  color=colourPred) +
      geom_line(lwd=lwd, colour=colourPred) + theme(legend.position="none")
    
    df2 <- data.frame(t=DataNew$input[[variable]], lower=lower, upper=upper)
    out <- out + geom_ribbon(data=df2, aes(x=t, ymin=lower, ymax=upper),
                             alpha=0.3, size=0, inherit.aes=FALSE)
    
    df3 <- data.frame(t=DataObs$input[[variable]], y=DataObs$response[[variable]])
    meltdf3 <- melt(df3, id="t")
    out <- out + geom_point(data=meltdf3, shape = 16, size = cex.points, 
                            color = colourData)
    
    ylab_i <- NULL
    if(is.null(ylabs)){
      ylab_i <- bquote(x[.(variable)])
    }else{
      if(length(ylabs)!=N)stop("length(ylabs) must be equal to the number of responses
                              to be plotted.")
      ylab_i <- ylabs[variable]
      
    }
    
    ylim_i <- NULL
    if(is.null(ylims)){
      ylim_i <- range(c(lower, upper))
    }else{
      stopMessage <- "'ylims' must  be a matrix with two columns; the number of rows must be equal to the number of responses to be plotted."
      if(!is.matrix(ylims)){
        stop(stopMessage)
      }else{
        if(nrow(ylims)!=N | ncol(ylims)!=2){
          stop(stopMessage)
        }
      }
      ylim_i <- ylims[[variable]]
    }
    
    if(!is.null(ylab_i)){
      out <- out + labs(y = ylab_i)
    }
    
    if(!is.null(ylim_i)){
      out <- out + coord_cartesian(ylim = ylim_i)
    }
    
    if(!is.null(cex.lab)){
      out <- out + theme(axis.title=element_text(size=cex.lab))
    }
    if(!is.null(cex.axis)){
      out <- out + theme(axis.text=element_text(size=cex.axis))
    }
    
    outList[[variable]] <- out
  }
  
  if(is.null(nCol)){
    nCol <- ifelse(N<4, N, ceiling(sqrt(N)))
  }
  
  outPlot <- do.call("grid.arrange", c(outList, ncol=nCol))

}






#' Plot auto- or cross-covariance function of a multivariate Gaussian process
#'
#' @param type Character. It can be either 'Cov' (for covariance function) or
#'   'Cor' (for corresponding correlation function).
#' @param output Integer identifying one element of the multivariate process.
#' @param outputp Integer identifying one element of the multivariate process.
#'   If 'output' and 'outputp' are the same, the auto-covariance function will
#'   be plotted. Otherwise, the cross-covariance function between 'output' and
#'   'outputp' will be plotted.
#' @param Data List containing data.
#' @param hp Vector of hyperparameters
#' @param idx Index vector identifying to which output the elements of
#'   concatenated vectors correspond to.
#' @param xlab Graphical parameter.
#' @param ylab Graphical parameter.
#' @param ylim Graphical parameter.
#' @param xlim Graphical parameter.
#' @param cex.lab Graphical parameter.
#' @param cex.axis Graphical parameter.
#' @import  ggplot2
#' @importFrom reshape2 melt
#' @return A plot.
#' @export
#' @examples
#' ## See examples in vignette:
#' # vignette("mgpr", package = "GPFDA")
plotmgpCovFun <- function(type="Cov", output, outputp, Data, hp, idx, 
                          xlab=NULL, ylab=NULL, 
                          ylim=NULL, xlim=NULL, cex.lab=10, cex.axis=10){
  
  Psi <- mgpCovMat(Data=Data, hp=hp)
  
  if(type=="Cor"){
    Psi <- cov2cor(Psi)
  }
  
  toPlot <- Psi[idx==output, idx==outputp][,1]
  tp0 <- Data$input[[outputp]][1]
  
  df <- data.frame(t=Data$input[[output]] - tp0, y=toPlot)
  meltdf <- melt(df, id="t")
  
  out <- ggplot(data=meltdf, aes(x=t,y=value,group=variable)) +
    geom_line() + theme(legend.position="none")
  
  
  if(is.null(ylab)){
    ylab <- bquote(.(type)*"["*X[.(output)]*"(t),"~X[.(outputp)]*"("*.(tp0)*")]")
  }
  if(is.null(xlab)){
    xlab <- bquote("t-("*.(tp0)*")")
  }
  
  
  if(!is.null(ylab)){
    out <- out + labs(y = ylab)
  }
  if(!is.null(xlab)){
    out <- out + labs(x = xlab)
  }
  if(!is.null(ylim)){
    out <- out + coord_cartesian(ylim = ylim)
  }
  if(!is.null(xlim)){
    out <- out + coord_cartesian(xlim = xlim)
  }
  if(!is.null(cex.lab)){
    out <- out + theme(axis.title=element_text(size=cex.lab))
  }
  if(!is.null(cex.axis)){
    out <- out + theme(axis.text=element_text(size=cex.axis))
  }
  
  suppressWarnings(print(out))
  
}




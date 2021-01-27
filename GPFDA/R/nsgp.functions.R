
#' Estimation of a nonseparable and/or nonstationary covariance structure (NSGPR
#' model)
#'
#' Estimate the covariance structure of a zero-mean Gaussian Process with
#' Q-dimensional input coordinates (covariates). \cr  \cr Multiple realisations
#' for the response variable can be used, provided they are observed on the same
#' grid of dimension n_1 x n_2 x ... x n_Q.\cr \cr Let n = n_1 x n_2 x ... x n_Q
#' and let nSamples be the number of realisations.

#' @param response Response variable. This should be a (n x nSamples) matrix
#'   where each column is a realisation
#' @param input List of Q input variables (see Details).
#' @param inputSubsetIdx A list identifying a subset of the input values to be
#'   used in the estimation (see Details).
#' @param corrModel Correlation function specification used for g(.). It can be
#'   either "pow.ex" or "matern".
#' @param gamma Power parameter used in powered exponential kernel function. It
#'   must be 0<gamma<=2.
#' @param nu Smoothness parameter of the Matern class. It must be a positive
#'   value.
#' @param whichTau Logical vector of dimension Q identifying which input
#'   coordinates the parameters are function of. For example, if Q=2 and
#'   parameters change only with respect to the first coordinate, then we set
#'   whichTau=c(T,F).
#' @param nBasis Number of B-spline basis functions in each coordinate direction
#'   along which parameters change.
#' @param cyclic Logical vector of dimension Q which defines which covariates
#'   are cyclic (periodic). For example, if basis functions should be cyclic
#'   only in the first coordinate direction, then cyclic=c(T,F). cyclic must
#'   have the same dimension of whichTau. If cyclic is TRUE for some coordinate
#'   direction, then cyclic B-spline functions will be used and the varying
#'   parameters (and their first two derivatives) will match at the boundaries
#'   of that coordinate direction.
#' @param unitSignalVariance Logical. TRUE if we assume realisations have
#'   variance 1. This is useful when we want to estimate an NSGP correlation
#'   function.
#' @param zeroNoiseVariance Logical. TRUE if we assume the realisations are
#'   noise-free.
#' @param sepCov Logical. TRUE only if we fix to zero all off-diagonal elements
#'   of the varying anisotropy matrix. Default to FALSE, allowing for a
#'   separable covariance function.
#' @param nInitCandidates number of initial hyperparameter vectors which are
#'   used to evaluate the log-likelihood function at a first step. After
#'   evaluating the log-likelihood using these 'nInitCandidates' vectors, the
#'   optimisation via nlminb() begins with the best of these vectors.
#' @param absBounds lower and upper boundaries for B-spline coefficients (if
#'   wanted).
#'
#' @importFrom mgcv cSplineDes
#' @importFrom splines bs
#' @details The input argument for Q=2 can be constructed as follows: \describe{
#'   \item{}{n1 <- 10} \item{}{n2 <- 1000} \item{}{input <- list()}
#'   \item{}{input[[1]] <- seq(0,1,length.out = n1)} \item{}{input[[2]] <-
#'   seq(0,1,length.out = n2)} } If we want to use every third lattice point in
#'   the second input variable (using Subset of Data), then we can set
#'   \describe{ \item{}{inputSubsetIdx <- list()} \item{}{inputSubsetIdx[[1]] <-
#'   1:n1} \item{}{inputSubsetIdx[[2]] <- seq(1,n2, by=3)} }
#'
#' @references Konzen, E., Shi, J. Q. and Wang, Z. (2020) "Modeling
#'   Function-Valued Processes with Nonseparable and/or Nonstationary Covariance
#'   Structure" <arXiv:1903.09981>.
#'
#' @return A list containing:  \describe{ \item{MLEsts}{Maximum likelihood
#'   estimates of B-spline coefficients and noise variance.}
#'   \item{response}{Matrix of response.} \item{inputMat}{Input coordinates in a
#'   matrix form} \item{corrModel}{Correlation function specification used for
#'   g(.)} }
#' @export
#' @examples
#' ## See examples in vignette:
#' # vignette("nsgpr", package = "GPFDA")
nsgpr <- function( response,
                   input,
                   corrModel="pow.ex",
                   gamma=2,
                   nu=1.5,
                   whichTau=NULL,
                   nBasis=5,
                   cyclic=NULL,
                   unitSignalVariance=F,
                   zeroNoiseVariance=F, 
                   sepCov=F,
                   nInitCandidates=300,
                   absBounds=6,
                   inputSubsetIdx=NULL){
  
  if(!is.list(input)){
    stop("The argument 'input' must be a list with Q elements")
  }
  n <- prod(sapply(input, length))
  
  response <- as.matrix(response)
  
  # number of realisations
  nSamples <- ncol(response)
  if(nrow(response)!=n){
    stop("The argument 'response' must be a vector of n elements or a matrix 
         with n rows, where n = prod(sapply(input, length))")
  }
  
  
  # number of input variables
  Q <- length(input)

  if(length(whichTau)!=Q){
    stop("whichTau must be a vector with length(input) elements")
  }
  if(length(cyclic)!=Q){
    stop("cyclic must be a vector with length(input) elements")
  }
  
  
  if(is.null(inputSubsetIdx)){
    inputSubset <- input
    inputMat <- as.matrix(expand.grid(input))
    
    inputIdx <- lapply(input, function(i) 1:length(i))
    inputSubsetIdx <- inputIdx
    inputIdxMat <- expand.grid(inputIdx)
  }else{
    
    inputSubset <- list()
    whichSubsetList <- list()
    for(q in 1:Q){
      inputSubset[[q]] <- input[[q]][ inputSubsetIdx[[q]] ]
      whichSubsetList[[q]] <- inputIdx[[q]]%in%inputSubsetIdx[[q]]
    }
    n <- prod(sapply(inputSubset, length))
    
    whichSubsetMat <- expand.grid(whichSubsetList)
    whichSubset <- apply(whichSubsetMat, 1, function(j){sum(j)==Q})
    
    inputMat <- as.matrix(expand.grid(inputSubset))
    
    inputIdx <- lapply(inputSubset, function(i) 1:length(i))
    inputIdxMat <- expand.grid(inputSubsetIdx)
    response <- response[whichSubset,]
    
  }

  if( is.null(whichTau) ){
  cat("\nPlease specify 'whichTau' 
    (i.e. the logical vector of length Q saying which are the 'tau' dimensions).
    \n")}
  
  cat(paste0("Parameters are varying in coordinate direction(s): "), 
      which(whichTau==T), " \n")
  
  if(nBasis<5)stop("nBasis must be >= 5")
  
  nTaus <- sum(whichTau)
  if(!(nTaus%in%c(1,2))){
    stop("The code works only if dimension of tau is 1 or 2")
  }
  
  if(is.null(cyclic)){
    cyclic <- rep(F, Q)
    cat(paste0("'cyclic' not used in any coordinate direction \n"))
  }else{
    if(length(cyclic)!=length(whichTau)){
      stop("'cyclic' and 'whichTau' must have same length")
    }
  }
  
  #===========================================================================
  # Specify lower, upper, and initial parameter values for nlminb()
  #===========================================================================
  
  num_betas <- nBasis^nTaus
  
  num_coeffs_omegas <- num_betas*Q*(Q+1)/2
  total_num_coeffs <- num_coeffs_omegas + num_betas 
  # num omegas + num betas for logsig2
  loc_coeffs <- matrix(1:total_num_coeffs, byrow=T, ncol=num_betas) 
  # last row is for logsig2
  
  coeffs_omegas_LB <- rep(-absBounds, num_coeffs_omegas)
  coeffs_omegas_UB <- rep(absBounds, num_coeffs_omegas)
  if(sepCov){
    whichZero <- (num_coeffs_omegas-num_betas+1):num_coeffs_omegas
    coeffs_omegas_LB[whichZero] <- 0
    coeffs_omegas_UB[whichZero] <- 0
  }
  
  coeffs_logsig2_LB <- rep(-absBounds, num_betas)
  coeffs_logsig2_UB <- rep(absBounds, num_betas)
  if(unitSignalVariance){
    coeffs_logsig2_LB <- rep(0, num_betas)
    coeffs_logsig2_UB <- rep(0, num_betas)
  }

  log_vareps_LB <- log(1e-05)
  log_vareps_UB <- log(10)
  if(zeroNoiseVariance){
    log_vareps_LB <- log(1e-10)
    log_vareps_UB <- log(1e-10)
  }
  lower <- c(coeffs_omegas_LB, coeffs_logsig2_LB, log_vareps_LB)
  upper <- c(coeffs_omegas_UB, coeffs_logsig2_UB, log_vareps_UB)
  
  
  if(nTaus==1){
    
   tauVec <- input[[which(whichTau==T)]]
   tauVecSub <- inputSubset[[which(whichTau==T)]]
   lengthTaus <- length(tauVecSub)
   if(cyclic[which(whichTau==T)]){
     numIntKnots <- nBasis-1
     quantiles_tau <- quantile(tauVec, 
                               probs=seq(0, 1, length.out=2+numIntKnots))
     bspl <- cSplineDes(x=tauVecSub, knots=quantiles_tau, ord=4, derivs=0)
   }else{
     numIntKnots <- nBasis-4
     quantiles_tau <- quantile(tauVec, 
                               probs=seq(0, 1, length.out=2+numIntKnots))
     quantiles_tau <- quantiles_tau[c(-1, -length(quantiles_tau))]
     bspl <- bs(tauVecSub, knots=quantiles_tau, intercept=T)
   }
   
  }else{
    
    bspl_j <- list()
    lengthTaus <- rep(NA, nTaus)
    for(j in 1:nTaus){
      tauVec <- input[which(whichTau==T)][[j]]
      tauVecSub <- inputSubset[which(whichTau==T)][[j]]
      lengthTaus[j] <- length(tauVecSub)
      if(cyclic[j]){
        numIntKnots <- nBasis-1
        quantiles_tau <- quantile(tauVec, 
                                  probs=seq(0, 1, length.out=2+numIntKnots))
        bspl_j[[j]] <- cSplineDes(x=tauVecSub, knots=quantiles_tau, ord=4, 
                                  derivs=0)
      }else{
        numIntKnots <- nBasis-4
        quantiles_tau <- quantile(tauVec, 
                                  probs=seq(0, 1, length.out=2+numIntKnots))
        quantiles_tau <- quantiles_tau[c(-1, -length(quantiles_tau))]
        bspl_j[[j]] <- bs(tauVecSub, knots=quantiles_tau, intercept=T)
      }
    }
    
    bspl <- kronecker(bspl_j[[2]], bspl_j[[1]])
  }
  
  tauVecIdx <- 1:prod(lengthTaus)
  
  n_hp <- length(lower)
  
  candidates <- matrix(0, nInitCandidates, n_hp)
  for(ipar in 1:n_hp){
    candidates[,ipar] <- runif(n=nInitCandidates, min=lower[ipar], 
                               max=upper[ipar] )
  }

  resCand <- apply(candidates, 1, function(x) 
    LogLikNSGP(hp=x, response=response, inputMat=inputMat, 
               inputIdxMat=inputIdxMat,
           inputSubsetIdx=inputSubsetIdx, bspl=bspl, 
           loc_coeffs=loc_coeffs, tauVecIdx=tauVecIdx, 
           corrModel=corrModel, gamma=gamma, nu=nu, whichTau=whichTau))

  init_opt <- candidates[which.min(resCand),]
  
  MLEs <- nlminb( start=init_opt,
                  objective=LogLikNSGP,
                  lower=lower,
                  upper=upper,
                  response=response, inputMat=inputMat, inputIdxMat=inputIdxMat, 
                  inputSubsetIdx=inputSubsetIdx, bspl=bspl, 
                  loc_coeffs=loc_coeffs, tauVecIdx=tauVecIdx, 
                  corrModel=corrModel, gamma=gamma, nu=nu, whichTau=whichTau)
  
  output <- list( MLEsts=MLEs$par,
                  response=response,
                  inputMat=inputMat,
                  corrModel=corrModel)

  return(output)
}





#' Calculate a NSGP covariance matrix given a vector of hyperparameters
#'
#'
#' @inheritParams nsgpr
#' @param hp Vector of hyperparameters estimated by function nsgpr.
#' @param calcCov Logical. Calculate covariance matrix or not. If FALSE, time or
#'   spatially-varying parameters are still provided.
#' @importFrom mgcv cSplineDes
#' @importFrom splines bs
#' @references Konzen, E., Shi, J. Q. and Wang, Z. (2020) "Modeling
#'   Function-Valued Processes with Nonseparable and/or Nonstationary Covariance
#'   Structure" <arXiv:1903.09981>.
#' @return A list containing  \describe{ 
#' \item{Cov}{Covariance matrix}
#' \item{vareps}{Noise variance}
#' \item{As_perTau}{List of varying anisotropy matrix over the input space}
#' \item{sig2_perTau}{Vector of signal variance over the input space}
#' }
#' @export
#' @examples
#' ## See examples in vignette:
#' # vignette("nsgpr", package = "GPFDA")
nsgpCovMat <- function(hp, input, inputSubsetIdx=NULL, nBasis=5, 
                       corrModel=corrModel, gamma=NULL, nu=NULL, cyclic=NULL, 
                       whichTau=NULL, calcCov=T){
  
  if(!is.list(input)){
    stop("The argument 'input' must be a list with Q elements")
  }
  n <- prod(sapply(input, length))
  
  Q <- length(input)
  
  if(is.null(inputSubsetIdx)){
    
    inputSubset <- input
    inputMat <- as.matrix(expand.grid(input))
    
    inputIdx <- lapply(input, function(i) 1:length(i))
    inputIdxMat <- expand.grid(inputIdx)
    
    inputSubsetIdx <- inputIdx
  }else{
    
    inputSubset <- list()
    whichSubsetList <- list()
    for(q in 1:Q){
      inputSubset[[q]] <- input[[q]][ inputSubsetIdx[[q]] ]
      whichSubsetList[[q]] <- inputIdx[[q]]%in%inputSubsetIdx[[q]]
    }
    n <- prod(sapply(inputSubset, length))
    
    whichSubsetMat <- expand.grid(whichSubsetList)
    whichSubset <- apply(whichSubsetMat, 1, function(j){sum(j)==Q})
    
    inputMat <- as.matrix(expand.grid(inputSubset))
    
    inputIdx <- lapply(inputSubset, function(i) 1:length(i))
    inputIdxMat <- expand.grid(inputSubsetIdx)
    
  }
  
  
  if( is.null(whichTau) ){
    cat("\nPlease specify 'whichTau' 
    (i.e. the logical vector of length Q saying which are the 'tau' dimensions).
        \n")}

  
  if(nBasis<5)stop("nBasis must be >= 5")
  
  
  nTaus <- sum(whichTau)
  if(!(nTaus%in%c(1,2))){
    stop("The code works only if dimension of tau is 1 or 2")
  }
  
  if(is.null(cyclic)){
    cyclic <- rep(F, Q)
    cat(paste0("'cyclic' not used in any coordinate direction \n"))
  }else{
    if(length(cyclic)!=length(whichTau)){
      stop("'cyclic' and 'whichTau' must have same length")
    }
  }
  
  
  if(nTaus==1){
    
    tauVec <- input[[which(whichTau==T)]]
    tauVecSub <- inputSubset[[which(whichTau==T)]]
    lengthTaus <- length(tauVecSub)
    if(cyclic[which(whichTau==T)]){
      numIntKnots <- nBasis-1
      quantiles_tau <- quantile(tauVec, 
                                probs=seq(0, 1, length.out=2+numIntKnots))
      bspl <- cSplineDes(x = tauVecSub, knots=quantiles_tau, ord=4, derivs=0)
    }else{
      numIntKnots <- nBasis-4
      quantiles_tau <- quantile(tauVec, 
                                probs=seq(0, 1, length.out=2+numIntKnots))
      quantiles_tau <- quantiles_tau[c(-1, -length(quantiles_tau))]
      bspl <- bs(tauVecSub, knots=quantiles_tau, intercept=T)
    }
    
  }else{
    
    bspl_j <- list()
    lengthTaus <- rep(NA, nTaus)
    for(j in 1:nTaus){
      tauVec <- input[which(whichTau==T)][[j]]
      tauVecSub <- inputSubset[which(whichTau==T)][[j]]
      lengthTaus[j] <- length(tauVecSub)
      if(cyclic[j]){
        numIntKnots <- nBasis-1
        quantiles_tau <- quantile(tauVec, 
                                  probs=seq(0, 1, length.out=2+numIntKnots))
        bspl_j[[j]] <- cSplineDes(x=tauVecSub, knots=quantiles_tau, 
                                  ord=4, derivs=0)
      }else{
        numIntKnots <- nBasis-4
        quantiles_tau <- quantile(tauVec, 
                                  probs=seq(0, 1, length.out=2+numIntKnots))
        quantiles_tau <- quantiles_tau[c(-1, -length(quantiles_tau))]
        bspl_j[[j]] <- bs(tauVecSub, knots=quantiles_tau, intercept=T)
      }
    }
    
    bspl <- kronecker(bspl_j[[2]], bspl_j[[1]])
  }
  
  tauVecIdx <- 1:prod(lengthTaus)
  
###################################
  vareps <- exp(hp[length(hp)])

  num_betas <- nBasis^nTaus
  
  num_coeffs_omegas <- num_betas*Q*(Q+1)/2
  total_num_coeffs <- num_coeffs_omegas + num_betas 
  # num omegas + num betas for logsig2
  loc_coeffs <- matrix(1:total_num_coeffs, byrow=T, ncol=num_betas) 
  # last row is for logsig2
  
  if(Q==1){
    omega1 <- bspl%*%hp[loc_coeffs[1,]]
    logsig2 <- bspl%*%hp[loc_coeffs[2,]]
    As_perTau <- list()
    for(jtau in seq_along(tauVecIdx)){
      As_perTau[[jtau]] <- as.matrix(exp(omega1[jtau]))
    }
  }
  
  if(Q==2){
    omega1 <- bspl%*%hp[loc_coeffs[1,]] 
    omega2 <- bspl%*%hp[loc_coeffs[2,]]
    omega3 <- bspl%*%hp[loc_coeffs[3,]]
    logsig2 <- bspl%*%hp[loc_coeffs[4,]]
    
    As_perTau <- list()
    for(jtau in seq_along(tauVecIdx)){
      As_perTau[[jtau]] <- CalcA_Q2(theta=c(omega1[jtau], omega2[jtau], 
                                            omega3[jtau]))  
    }
  }
  
  if(Q==3){
    omega1 <- bspl%*%hp[loc_coeffs[1,]]
    omega2 <- bspl%*%hp[loc_coeffs[2,]]
    omega3 <- bspl%*%hp[loc_coeffs[3,]]
    omega4 <- bspl%*%hp[loc_coeffs[4,]]
    omega5 <- bspl%*%hp[loc_coeffs[5,]]
    omega6 <- bspl%*%hp[loc_coeffs[6,]]
    logsig2 <- bspl%*%hp[loc_coeffs[7,]]
    
    As_perTau <- list()
    for(jtau in seq_along(tauVecIdx)){
      As_perTau[[jtau]] <- CalcA_Q3(theta=c(omega1[jtau], omega2[jtau], 
                                            omega3[jtau], omega4[jtau], 
                                            omega5[jtau], omega6[jtau]))  
    }
  }
  
  obs_variances_perTau <- exp(logsig2)
  
  if(calcCov==T){
    nTaus <- sum(whichTau)
    if(nTaus==1){
      A_List <- list()
      obs_variance <- rep(0, n)
      for(n_i in 1:n){
        
        whitau <- which(whichTau==T)
        k <- inputIdxMat[n_i,whitau]
        k <- which(k==inputSubsetIdx[[whitau]])
        # k <- which(k==unique(inputIdxMat[,which(whichTau==T)]))
        A_List[[n_i]] <- As_perTau[[k]]
        obs_variance[n_i] <- c(obs_variances_perTau[k])
      }
    }
    if(nTaus==2){
      if(Q!=2)stop("Q must equal 2 if nTaus=2")
      A_List <- As_perTau
      obs_variance <- c(obs_variances_perTau)
    }
    
    ScaleDistMats <- calcScaleDistMats(A_List = A_List, coords = inputMat)
    
    Scale.mat <- ScaleDistMats$Scale.mat
    Dist.mat <- ScaleDistMats$Dist.mat
    
    UnscalCorr <- unscaledCorr( Dist.mat=Dist.mat, corrModel=corrModel, 
                                gamma=gamma, nu=nu)

    NS.corr <- Scale.mat*UnscalCorr
    
    
    Cov <- diag( sqrt(obs_variance) ) %*% NS.corr %*% diag(sqrt(obs_variance))
    diag(Cov) <- diag(Cov) + vareps + 1e-8
  }else{
    Cov=NULL
  }
  
  return(list(Cov=Cov, vareps=vareps, As_perTau=As_perTau, 
              sig2_perTau=obs_variances_perTau))
}




LogLikNSGP <- function(hp, response, inputMat, inputIdxMat, inputSubsetIdx, 
                       bspl, loc_coeffs, tauVecIdx, corrModel, gamma, nu, 
                       whichTau){
  
  n <- nrow(inputMat)
  Q <- ncol(inputMat)
  nrep <- length(response)/n
  
  vareps <- exp(hp[length(hp)])

  if(Q==1){
    omega1 <- bspl%*%hp[loc_coeffs[1,]]
    logsig2 <- bspl%*%hp[loc_coeffs[2,]]
    As_perTau <- list()
    for(jtau in seq_along(tauVecIdx)){
      As_perTau[[jtau]] <- as.matrix(exp(omega1[jtau]))
    }
  }
  
  if(Q==2){
    omega1 <- bspl%*%hp[loc_coeffs[1,]] 
    omega2 <- bspl%*%hp[loc_coeffs[2,]]
    omega3 <- bspl%*%hp[loc_coeffs[3,]]
    logsig2 <- bspl%*%hp[loc_coeffs[4,]]
    
    As_perTau <- list()
    for(jtau in seq_along(tauVecIdx)){
      As_perTau[[jtau]] <- CalcA_Q2(theta=c(omega1[jtau], omega2[jtau], 
                                            omega3[jtau]))  
    }
  }
  
  if(Q==3){
    omega1 <- bspl%*%hp[loc_coeffs[1,]]
    omega2 <- bspl%*%hp[loc_coeffs[2,]]
    omega3 <- bspl%*%hp[loc_coeffs[3,]]
    omega4 <- bspl%*%hp[loc_coeffs[4,]]
    omega5 <- bspl%*%hp[loc_coeffs[5,]]
    omega6 <- bspl%*%hp[loc_coeffs[6,]]
    logsig2 <- bspl%*%hp[loc_coeffs[7,]]
    
    As_perTau <- list()
    for(jtau in seq_along(tauVecIdx)){
      As_perTau[[jtau]] <- CalcA_Q3(theta=c(omega1[jtau], omega2[jtau], 
                                            omega3[jtau], omega4[jtau], 
                                            omega5[jtau], omega6[jtau]))  
    }
  }

  obs_variances_perTau <- exp(logsig2)

  nTaus <- sum(whichTau)
  if(nTaus==1){
    A_List <- list()
    obs_variance <- rep(0, n)
    for(n_i in 1:n){
      
      whitau <- which(whichTau==T)
      k <- inputIdxMat[n_i,whitau]
      k <- which(k==inputSubsetIdx[[whitau]])
      
      A_List[[n_i]] <- As_perTau[[k]]
      obs_variance[n_i] <- c(obs_variances_perTau[k])
    }
  }
  if(nTaus==2){
    if(Q!=2)stop("Q must equal 2 if nTaus=2")
    A_List <- As_perTau
    obs_variance <- c(obs_variances_perTau)
  }
  
  ScaleDistMats <- calcScaleDistMats(A_List = A_List, coords = inputMat)
  
  Scale.mat <- ScaleDistMats$Scale.mat
  Dist.mat <- ScaleDistMats$Dist.mat
  
  UnscalCorr <- unscaledCorr( Dist.mat=Dist.mat, corrModel=corrModel, 
                              gamma=gamma, nu=nu)
  NS.corr <- Scale.mat*UnscalCorr
  
  
  Cov <- diag( sqrt(obs_variance) ) %*% NS.corr %*% diag(sqrt(obs_variance))
  diag(Cov) <- diag(Cov) + vareps + 1e-8
  
  cholA <- chol(Cov)
  yt.invK.y <- t(response)%*%chol2inv(cholA)%*%response
  logdetK <- 2*sum(log(diag(cholA)))
  
  if(nrep==1){
    fX <- 0.5*logdetK + 0.5*yt.invK.y + 0.5*n*log(2*pi)
  }else{
    fX <- nrep*0.5*logdetK + 0.5*sum(diag( yt.invK.y )) + nrep*0.5*n*log(2*pi)
  }
  fX <- as.numeric(fX)
  
  return(fX)
}




#' Prediction of NSGPR model
#'
#' @inheritParams nsgpr
#' @param hp Vector of hyperparameters estimated by function nsgpr.
#' @param inputNew List of Q test set input variables.
#' @param noiseFreePred Logical.  If TRUE, predictions will be noise-free.
#'
#' @references Konzen, E., Shi, J. Q. and Wang, Z. (2020) "Modeling
#'   Function-Valued Processes with Nonseparable and/or Nonstationary Covariance
#'   Structure" <arXiv:1903.09981>.
#' @return A list containing  \describe{ 
#' \item{pred.mean}{Mean of predictions for the test set.}
#' \item{pred.sd}{Standard deviation of predictions for the test set.}
#' \item{noiseFreePred}{Logical. If TRUE, predictions are noise-free.}
#' }
#' 
#' @export
#' @examples
#' ## See examples in vignette:
#' # vignette("nsgpr", package = "GPFDA")
nsgprPredict <- function(hp, response, input, inputNew, noiseFreePred=F, 
                           nBasis=nBasis, corrModel=corrModel, gamma=gamma, 
                           nu=nu, cyclic=cyclic, whichTau=whichTau){
  
  if(is.null(inputNew)){
    inputNew <- input
  }
  
  Kobs <- nsgpCovMat(hp=hp, input=input, inputSubsetIdx=NULL,
                        nBasis=nBasis, corrModel=corrModel, gamma=gamma, nu=nu,
                        cyclic=cyclic, whichTau=whichTau, calcCov=T)$Cov
  invQ <- chol2inv(chol(Kobs))
  
  Q1 <- nsgpCovMatAsym(hp=hp, input=input, inputNew=inputNew, nBasis=nBasis, 
                        corrModel=corrModel, gamma=gamma, nu=nu, cyclic=cyclic, 
                        whichTau=whichTau)
  # response is a (n x nSamples) matrix
  mu <- t(Q1)%*%invQ%*%response
  
  Qstar <- nsgpCovMat(hp=hp, input=inputNew, inputSubsetIdx=NULL,
                         nBasis=nBasis, corrModel=corrModel, gamma=gamma, nu=nu,
                         cyclic=cyclic, whichTau=whichTau, calcCov=T)$Cov
  
  if(noiseFreePred){
    sigma2 <- diag(Qstar) - diag(t(Q1)%*%invQ%*%Q1)
  }else{
    sigma2 <- diag(Qstar) - diag(t(Q1)%*%invQ%*%Q1) + exp(hp[length(hp)])
  }
  pred.sd <- sqrt(sigma2)
  
  result <- c(list('pred.mean'=mu,
                   'pred.sd'=pred.sd,
                   'noiseFreePred'=noiseFreePred))
  
  return(result)
}



#' Calculate an asymmetric NSGP covariance matrix
#' 
#' @inheritParams nsgpCovMat 
#' 
#' @param inputNew  List of Q test set input variables.
#'
#' @importFrom mgcv cSplineDes
#' @importFrom splines bs
#' 
#' @references Konzen, E., Shi, J. Q. and Wang, Z. (2020) "Modeling
#'   Function-Valued Processes with Nonseparable and/or Nonstationary Covariance
#'   Structure" <arXiv:1903.09981>.
#' @return An asymmetric covariance matrix
#' @export
nsgpCovMatAsym <- function(hp, input, inputNew, nBasis=5, corrModel=corrModel, 
                            gamma=NULL, nu=NULL, cyclic=NULL, whichTau=NULL){
  
  if(!is.list(input)){
    stop("The argument 'input' must be a list with Q elements")
  }
  if(!is.list(inputNew)){
    stop("The argument 'inputNew' must be a list with Q elements")
  }
  n <- prod(sapply(input, length))
  nStar <- prod(sapply(inputNew, length))
  
  # number of input variables
  Q <- length(input)
  
  inputSubset <- input
  inputSubsetStar <- inputNew
  inputMat <- as.matrix(expand.grid(input))
  inputMatStar <- as.matrix(expand.grid(inputNew))
  
  inputIdx <- lapply(input, function(i) 1:length(i))
  inputIdxMat <- expand.grid(inputIdx)
  inputIdxStar <- lapply(inputNew, function(i) 1:length(i))
  inputIdxMatStar <- expand.grid(inputIdxStar)
  
  
  if( is.null(whichTau) ){
    cat("\nPlease specify 'whichTau' 
    (i.e. the logical vector of length Q saying which are the 'tau' dimensions).
        \n")}
  
  if(nBasis<5)stop("nBasis must be >= 5")
  
  nTaus <- sum(whichTau)
  if(!(nTaus%in%c(1,2))){
    stop("The code works only if dimension of tau is 1 or 2")
  }
  
  if(is.null(cyclic)){
    cyclic <- rep(F, Q)
    cat(paste0("'cyclic' not used in any coordinate direction \n"))
  }else{
    if(length(cyclic)!=length(whichTau)){
      stop("'cyclic' and 'whichTau' must have same length")
    }
  }
  
  vareps <- exp(hp[length(hp)])
  
  #######################################################################
  ### Calculate obs_variance and A_List
  
  if(nTaus==1){
    
    tauVec <- input[[which(whichTau==T)]]
    tauVecSub <- inputSubset[[which(whichTau==T)]]
    lengthTaus <- length(tauVecSub)
    if(cyclic[which(whichTau==T)]){
      numIntKnots <- nBasis-1
      quantiles_tau <- quantile(tauVec, 
                                probs=seq(0, 1, length.out=2+numIntKnots))
      bspl <- cSplineDes(x = tauVecSub, knots=quantiles_tau, ord=4, derivs=0)
    }else{
      numIntKnots <- nBasis-4
      quantiles_tau <- quantile(tauVec, 
                                probs=seq(0, 1, length.out=2+numIntKnots))
      quantiles_tau <- quantiles_tau[c(-1, -length(quantiles_tau))]
      bspl <- bs(tauVecSub, knots=quantiles_tau, intercept=T)
    }
    
  }else{
    
    bspl_j <- list()
    lengthTaus <- rep(NA, nTaus)
    for(j in 1:nTaus){
      tauVec <- input[which(whichTau==T)][[j]]
      tauVecSub <- inputSubset[which(whichTau==T)][[j]]
      lengthTaus[j] <- length(tauVecSub)
      if(cyclic[j]){
        numIntKnots <- nBasis-1
        quantiles_tau <- quantile(tauVec, 
                                  probs=seq(0, 1, length.out=2+numIntKnots))
        bspl_j[[j]] <- cSplineDes(x = tauVecSub, knots=quantiles_tau, 
                                  ord=4, derivs=0)
      }else{
        numIntKnots <- nBasis-4
        quantiles_tau <- quantile(tauVec, 
                                  probs=seq(0, 1, length.out=2+numIntKnots))
        quantiles_tau <- quantiles_tau[c(-1, -length(quantiles_tau))]
        bspl_j[[j]] <- bs(tauVecSub, knots=quantiles_tau, intercept=T)
      }
    }
    
    bspl <- kronecker(bspl_j[[2]], bspl_j[[1]])
  }
  
  tauVecIdx <- 1:prod(lengthTaus)
  
  ###################################
  
  
  num_betas <- nBasis^nTaus
  
  num_coeffs_omegas <- num_betas*Q*(Q+1)/2
  total_num_coeffs <- num_coeffs_omegas + num_betas
  loc_coeffs <- matrix(1:total_num_coeffs, byrow=T, ncol=num_betas)
  
  if(Q==1){
    omega1 <- bspl%*%hp[loc_coeffs[1,]]
    logsig2 <- bspl%*%hp[loc_coeffs[2,]]
    As_perTau <- list()
    for(jtau in seq_along(tauVecIdx)){
      As_perTau[[jtau]] <- as.matrix(exp(omega1[jtau]))
    }
  }
  
  if(Q==2){
    omega1 <- bspl%*%hp[loc_coeffs[1,]]
    omega2 <- bspl%*%hp[loc_coeffs[2,]]
    omega3 <- bspl%*%hp[loc_coeffs[3,]]
    logsig2 <- bspl%*%hp[loc_coeffs[4,]]
    
    As_perTau <- list()
    for(jtau in seq_along(tauVecIdx)){
      As_perTau[[jtau]] <- CalcA_Q2(theta=c(omega1[jtau], omega2[jtau], 
                                            omega3[jtau]))  
    }
  }
  
  if(Q==3){
    omega1 <- bspl%*%hp[loc_coeffs[1,]]
    omega2 <- bspl%*%hp[loc_coeffs[2,]]
    omega3 <- bspl%*%hp[loc_coeffs[3,]]
    omega4 <- bspl%*%hp[loc_coeffs[4,]]
    omega5 <- bspl%*%hp[loc_coeffs[5,]]
    omega6 <- bspl%*%hp[loc_coeffs[6,]]
    logsig2 <- bspl%*%hp[loc_coeffs[7,]]
    
    As_perTau <- list()
    for(jtau in seq_along(tauVecIdx)){
      As_perTau[[jtau]] <- CalcA_Q3(theta=c(omega1[jtau], omega2[jtau], 
                                            omega3[jtau], omega4[jtau], 
                                            omega5[jtau], omega6[jtau]))  
    }
  }
  
  obs_variances_perTau <- exp(logsig2)

  
  nTaus <- sum(whichTau)
  if(nTaus==1){
    A_List <- list()
    obs_variance <- rep(0, n)
    for(n_i in 1:n){
      
      whitau <- which(whichTau==T)
      k <- inputIdxMat[n_i,whitau]
      k <- which(k==inputIdx[[whitau]])
      # k <- which(k==unique(inputIdxMat[,which(whichTau==T)]))
      A_List[[n_i]] <- As_perTau[[k]]
      obs_variance[n_i] <- c(obs_variances_perTau[k])
    }
  }
  if(nTaus==2){
    if(Q!=2)stop("Q must equal 2 if nTaus=2")
    A_List <- As_perTau
    obs_variance <- c(obs_variances_perTau)
  }
  
  
  ### finish calculation of obs_variance and A_List
  #######################################################################
  
  
  #######################################################################
  ### Calculate obs_varianceList and Astar_List
  
  
  if(nTaus==1){
    
    tauVec <- inputNew[[which(whichTau==T)]]
    tauVecSub <- inputSubsetStar[[which(whichTau==T)]]
    lengthTaus <- length(tauVecSub)
    if(cyclic[which(whichTau==T)]){
      numIntKnots <- nBasis-1
      quantiles_tau <- quantile(tauVec, 
                                probs=seq(0, 1, length.out=2+numIntKnots))
      bspl <- cSplineDes(x = tauVecSub, knots=quantiles_tau, ord=4, derivs=0)
    }else{
      numIntKnots <- nBasis-4
      quantiles_tau <- quantile(tauVec, 
                                probs=seq(0, 1, length.out=2+numIntKnots))
      quantiles_tau <- quantiles_tau[c(-1, -length(quantiles_tau))]
      bspl <- bs(tauVecSub, knots=quantiles_tau, intercept=T)
    }
    
  }else{
    
    bspl_j <- list()
    lengthTaus <- rep(NA, nTaus)
    for(j in 1:nTaus){
      tauVec <- inputNew[which(whichTau==T)][[j]]
      tauVecSub <- inputSubsetStar[which(whichTau==T)][[j]]
      lengthTaus[j] <- length(tauVecSub)
      if(cyclic[j]){
        numIntKnots <- nBasis-1
        quantiles_tau <- quantile(tauVec, 
                                  probs=seq(0, 1, length.out=2+numIntKnots))
        bspl_j[[j]] <- cSplineDes(x=tauVecSub, knots=quantiles_tau, ord=4, 
                                  derivs=0)
      }else{
        numIntKnots <- nBasis-4
        quantiles_tau <- quantile(tauVec, 
                                  probs=seq(0, 1, length.out=2+numIntKnots))
        quantiles_tau <- quantiles_tau[c(-1, -length(quantiles_tau))]
        bspl_j[[j]] <- bs(tauVecSub, knots=quantiles_tau, intercept=T)
      }
    }
    
    bspl <- kronecker(bspl_j[[2]], bspl_j[[1]])
  }
  
  tauVecIdx <- 1:prod(lengthTaus)
  
  ###################################
  
  num_betas <- nBasis^nTaus
  
  num_coeffs_omegas <- num_betas*Q*(Q+1)/2
  total_num_coeffs <- num_coeffs_omegas + num_betas
  loc_coeffs <- matrix(1:total_num_coeffs, byrow=T, ncol=num_betas)
  
  if(Q==1){
    omega1 <- bspl%*%hp[loc_coeffs[1,]]
    logsig2 <- bspl%*%hp[loc_coeffs[2,]]
    As_perTau <- list()
    for(jtau in seq_along(tauVecIdx)){
      As_perTau[[jtau]] <- as.matrix(exp(omega1[jtau]))
    }
  }
  
  if(Q==2){
    omega1 <- bspl%*%hp[loc_coeffs[1,]] 
    omega2 <- bspl%*%hp[loc_coeffs[2,]]
    omega3 <- bspl%*%hp[loc_coeffs[3,]]
    logsig2 <- bspl%*%hp[loc_coeffs[4,]]
    
    As_perTau <- list()
    for(jtau in seq_along(tauVecIdx)){
      As_perTau[[jtau]] <- CalcA_Q2(theta=c(omega1[jtau], omega2[jtau], 
                                            omega3[jtau]))  
    }
  }
  
  if(Q==3){
    omega1 <- bspl%*%hp[loc_coeffs[1,]]
    omega2 <- bspl%*%hp[loc_coeffs[2,]]
    omega3 <- bspl%*%hp[loc_coeffs[3,]]
    omega4 <- bspl%*%hp[loc_coeffs[4,]]
    omega5 <- bspl%*%hp[loc_coeffs[5,]]
    omega6 <- bspl%*%hp[loc_coeffs[6,]]
    logsig2 <- bspl%*%hp[loc_coeffs[7,]]
    
    As_perTau <- list()
    for(jtau in seq_along(tauVecIdx)){
      As_perTau[[jtau]] <- CalcA_Q3(theta=c(omega1[jtau], omega2[jtau], 
                                            omega3[jtau], omega4[jtau], 
                                            omega5[jtau], omega6[jtau]))  
    }
  }
  
  obs_variances_perTau <- exp(logsig2)
  
  
  nTaus <- sum(whichTau)
  if(nTaus==1){
    Astar_List <- list()
    obs_varianceStar <- rep(0, nStar)
    for(n_i in 1:nStar){
      
      whitau <- which(whichTau==T)
      k <- inputIdxMatStar[n_i,whitau]
      k <- which(k==inputIdxStar[[whitau]])
      # k <- which(k==unique(inputIdxMatStar[,which(whichTau==T)]))
      Astar_List[[n_i]] <- As_perTau[[k]]
      obs_varianceStar[n_i] <- c(obs_variances_perTau[k])
    }
  }
  if(nTaus==2){
    if(Q!=2)stop("Q must equal 2 if nTaus=2")
    Astar_List <- As_perTau
    obs_varianceStar <- c(obs_variances_perTau)
  }
  
  
  ### finish calculation of obs_varianceStar and Astar_List
  #######################################################################
  
  # ScaleDistMats <- calcScaleDistMats(A_List = A_List, coords = inputMat)
  
  ScaleDistMats <- calcScaleDistMatsAsym(A_List=A_List, Astar_List=Astar_List, 
                                         coords=inputMat, 
                                         coordsStar=inputMatStar)
  Scale.mat <- ScaleDistMats$Scale.mat
  Dist.mat <- ScaleDistMats$Dist.mat
  
  UnscalCorr <- unscaledCorr( Dist.mat=Dist.mat, corrModel=corrModel, 
                              gamma=gamma, nu=nu)
  NS.corr <- Scale.mat*UnscalCorr
  
  
  Cov <- diag( sqrt(obs_variance) ) %*% NS.corr %*% diag(sqrt(obs_varianceStar))
  diag(Cov) <- diag(Cov) + vareps + 1e-8
  
  return(Cov=Cov)
}


#' Calculate an unscaled NSGP correlation matrix
#'
#' @inheritParams nsgpr
#' @param Dist.mat Distance matrix
#'
#' @references Konzen, E., Shi, J. Q. and Wang, Z. (2020) "Modeling
#'   Function-Valued Processes with Nonseparable and/or Nonstationary Covariance
#'   Structure" <arXiv:1903.09981>.
#' @return A matrix
#' @export
#' @examples
#' ## See examples in vignette:
#' # vignette("nsgpr", package = "GPFDA")
unscaledCorr <- function(Dist.mat, corrModel, gamma=NULL, nu=NULL){
  
  if(!(corrModel%in%c("pow.ex", "matern"))){
    stop("corrModel must be either 'pow.ex' or 'matern'")
  }
  if(corrModel=="pow.ex"){
    if(!(gamma >0 & gamma<=2)){
      stop("Parameter 'gamma' must be in (0,2]")
    }
  }
  if(corrModel=="matern"){
    if(!(nu >0)){
      stop("Parameter 'nu' must be positive")
    }
  }
  
  if(corrModel=="pow.ex"){
    Cor <- exp(-(Dist.mat)^gamma)  
  }
  if(corrModel=="matern"){
    besselMod <- besselK(x=Dist.mat, nu=nu)
    Cor <- (Dist.mat^nu)*besselMod/(gamma(nu)*(2^(nu-1)))
    diag(Cor) <- 1
  }
  
  return(Cor)
}

CalcA_Q2 <- function(theta){
  
  l1 <- c(exp(theta[1]), 0)
  l2 <- c(exp(theta[2]), pi*exp(theta[3])/(1+exp(theta[3])))
  
  L1 <- c(l1[1]*cos(l1[2]), 0)
  L2 <- c(l2[1]*cos(l2[2]), l2[1]*sin(l2[2]))
  
  myL <- cbind(L1,L2, deparse.level = 0)
  A <- crossprod(myL)
  diag(A) <- diag(A) + 1e-4
  
  return(A)
}

CalcA_Q3 <- function(theta){
  
  l1 <- c(exp(theta[1]),0,0)
  l2 <- c(exp(theta[2]),pi*exp(theta[4])/(1+exp(theta[4])),0)
  l3 <- c(exp(theta[3]),pi*exp(theta[5])/(1+exp(theta[5])),
          pi*exp(theta[6])/(1+exp(theta[6])))
  
  L1 <- c(l1[1]*cos(l1[2]), 0, 0)
  L2 <- c(l2[1]*cos(l2[2]), l2[1]*sin(l2[2])*cos(l2[3]), 0)
  L3 <- c(l3[1]*cos(l3[2]), l3[1]*sin(l3[2])*cos(l3[3]), 
          l3[1]*sin(l3[2])*sin(l3[3]))
  
  myL <- cbind(L1,L2,L3, deparse.level = 0)
  A <- crossprod(myL)
  diag(A) <- diag(A) + 1e-4
  
  return(A)
} 

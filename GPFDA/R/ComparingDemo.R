

    # rm(list = ls())
    
    library("GPFDA")
    require("MASS")
    
    library("useful")
    library("rbenchmark")
    library("Rcpp")
    library("inline")
    
    
    n <- 1000 # training set size + test set size
    nSamples <- 1
    # hp <- list('linear.a'=log(10), 'pow.ex.v'=log(5),'pow.ex.w'=log(10),'vv'=log(1))
    hp <- list('pow.ex.v'=log(2),'pow.ex.w'=log(100),'vv'=log(0.3))
    h <- seq(0,1,len=n)
    idx <- sort(sample(1:n,floor(n/2)))
    X <- as.matrix(h[idx])
    # Y <- (mvrnorm(n=n,mu=h-h,Sigma=(cov.linear(hp,h)+cov.pow.ex(hp,h)))[,1]) # *0.1+sin(h*6)
    Sigma <- cov.pow.ex(hp,h)
    diag(Sigma) <- diag(Sigma) + exp(hp$vv)
    
    if(nSamples==1){
      Y <- matrix(mvrnorm(n=nSamples,mu=h-h,Sigma=Sigma), ncol=1)
    }else{
      Y <- t(mvrnorm(n=nSamples,mu=h-h,Sigma=Sigma))
    }
    str(Y)
    str(X)
    
    
    # plot(h, Y)
    
    Y <- as.matrix(Y[idx,])
    
    str(Y)
    str(X)
    
    hyper <- list('pow.ex.v'=log(8),'pow.ex.w'=log(8),'vv'=log(0.02))
    
    x <- as.matrix(seq(0,1,by=0.01))
    t0 <- proc.time()
    a <- gpr(X,Y,c('pow.ex'),hyper=hyper,trace=2)
    GPFDA_time <- proc.time() - t0
    b <- gppredict(a,x)
    
    
    
    true_pars <- exp(unlist(hp))
    GPFDA_estimates <- round(exp(unlist(a$hyper)), 4)
    
    
    
    source("xixj_staNEW_sq.R")
    source("xixj_staNEW.R")
    source("xixjNEW_sq.R")
    source("xixjNEW.R")
    
    source("gp.functions5NEW.R")
    # source("gp.functions5NEW_TryCleanLogLik.R")
    
    
    # Data=X; response=Y; Cov=c('pow.ex'); # Cov=c('linear','pow.ex');
    # hyper=hyper;trace=2
    # NewHyper=NULL; mean=0; gamma=1;itermax=100;reltol=8e-10
    
    t0 <- proc.time()
    aNew <- gprNEW(Data=X, response=Y, Cov=c('pow.ex'),hyper=hyper,trace=2)
    newGPFDA_time <- proc.time() - t0
    
    exp(unlist(hp)) # true
    newGPFDA_estimates <- round(exp(unlist(aNew$hyper)), 4)
    
    bNew <- gppredictNEW(aNew,x)
    
    GPFDA:::plot.gpr(b)
    plot(bNew)
    text(0.8, 3, "new results", cex=2)
    
    max(abs(b$pred.mean - bNew$pred.mean))
    max(abs(b$pred.sd - bNew$pred.sd))
    
    
    
    m <- 100
    t0 <- proc.time()
    aNewSubset <- gprNEW(Data=X, response=Y, Cov=c('pow.ex'), m = m,
                   hyper=hyper, trace=2)
    Subset_newGPFDA_time <- proc.time() - t0
    
    exp(unlist(hp)) # true
    Subset_newGPFDA_estimates <- round(exp(unlist(aNewSubset$hyper)), 4)
    
    # train <- aNewSubset
    # Data.new <- x
    bNewSubset <- gppredictNEW(aNewSubset,x, noiseFreePred = F)
    plot(bNewSubset)
    
    bNewSubsetNoiseFree <- gppredictNEW(aNewSubset,x, noiseFreePred = T)
    plot(bNewSubsetNoiseFree)
    
    
    
    GPFDA:::plot.gpr(b)
    plot(bNew); mtext("new code", side = 3, at=0.85, line = 1, cex=2, col=2)
    plot(bNewSubset); mtext("new code: subset", side = 3, at=0.85, line = 1, cex=1.5, col=2)
    plot(bNewSubsetNoiseFree); mtext("new code: subset", side = 3, at=0.85, line = 1, cex=1.5, col=2)
    
    GPFDA_estimates
    newGPFDA_estimates
    Subset_newGPFDA_estimates
    
    GPFDA_time
    newGPFDA_time
    Subset_newGPFDA_time



# fixInNamespace("DCov.vv", pos="package:GPFDA")
# DCov.vv=function(hyper,Alpha,invQ){
#   # Dfvv=-sum(diag(invQ))*exp(hyper$vv) + t(Alpha)%*%Alpha*exp(hyper$vv)
#   Dfvv=0.5*sum(diag(Alpha%*%t(Alpha) - invQ))*exp(hyper$vv)
#   return(Dfvv)
# }
# Dfvv=-sum(diag(invQ))*exp(hyper$vv) + t(Alpha)%*%Alpha*exp(hyper$vv) ## OLD

# plot(a)
# plot(b)

# upper=b$pred.mean+1.96*b$pred.sd
# lower=b$pred.mean-1.96*b$pred.sd
# plot(-100,-100,col=0,xlim=range(x[,1]),ylim=c(min(upper,lower,Y)-0.1*abs(min(upper,lower,Y)),max(upper,lower,Y)+0.1*abs(max(upper,lower,Y))),main="Prediction", xlab="input ( x )",ylab="response")
# polygon(c(x[,1], rev(x[,1])), c(upper, rev(lower)),col = "grey60", border = NA)
# points(X[,1],Y,pch=4,col=2,cex=.8)
# lines(x[,1],b$pred.mean,col=4,lwd=1.5)

# SigmaT <- cov.pow.ex(hp,h)
# plot(h, SigmaT[1,], type="l", ylab="Cov[Y(0), Y(h)]", main="True covariance function of f")
# SigmaTN <- SigmaT
# diag(SigmaTN) <- diag(SigmaTN) + exp(hp$vv)
# plot(h, SigmaTN[1,], type="l", ylab="Cov[Y(0), Y(h)]", main="True covariance function of X")
# SigmaEst <- cov.pow.ex(a$hyper,h)
# plot(h, SigmaEst[1,], type="l", ylab="Cov[Y(0), Y(h)]", main="Estimated covariance function of f")
# SigmaEstN <- SigmaEst
# diag(SigmaEstN) <- diag(SigmaEstN) + exp(a$hyper$vv)
# plot(h, SigmaEstN[1,], type="l", ylab="Cov[Y(0), Y(h)]", main="Estimated covariance function of X")












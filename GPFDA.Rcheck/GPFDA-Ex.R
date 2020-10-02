pkgname <- "GPFDA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "GPFDA-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('GPFDA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("GPFDA-package")
### * GPFDA-package

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: GPFDA-package
### Title: Gaussian Process Regression for Functional Data Analysis
### Aliases: GPFDA-package GPFDA
### Keywords: package

### ** Examples

  ## Optional simple examples of the most important functions
  ## Use \dontrun{} around code to be shown but not executed



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("GPFDA-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("betaPar")
### * betaPar

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: betaPar
### Title: Create an fdPar object
### Aliases: betaPar

### ** Examples

library(GPFDA)
beta1=betaPar()
beta2=betaPar(list(nbasis=7,lambda=0.01))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("betaPar", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gpfr")
### * gpfr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gpfr
### Title: Gaussian Process for functional data.
### Aliases: gpfr

### ** Examples

library(GPFDA)

traindata=vector('list',20)
for(i in 1:20) traindata[[i]]=i
n=50
traindata=lapply(traindata,function(i) {
  x=seq(-3,3,len=n)
  y=sin(x^2)-x+0.2*rnorm(n,runif(1,-3,3),runif(1,0.5,3))
  x1=0.5*x^3+exp(x)+rnorm(n,runif(1,-3,3),runif(1,0.5,5))
  x2=cos(x^3)+0.2*rnorm(n,runif(1,-3,3),runif(1,0.5,5))
  mat=cbind(x,x1,x2,y)
  colnames(mat)=c('time','x1','x2','y')
  scale=t(c(2*(mean(y)>0.25)-1,(var(y)>3.6)*2-1,(sd(y)-sd(x)>1.4)*2-1))
  i=list(mat,scale)
})

lx=do.call('rbind',lapply(traindata,function(i)i[[2]]))
fx1=do.call('rbind',lapply(traindata,function(i)i[[1]][,2]))
fx2=do.call('rbind',lapply(traindata,function(i)i[[1]][,3]))
fy1=do.call('rbind',lapply(traindata,function(i)i[[1]][,4]))
time_old=traindata[[1]][[1]][,1]

## NOT RUN
# system.time(a1<-gpfr(response=(fy1),lReg=lx,fReg=NULL,gpReg=list(fx1,fx2)
#                ,fyList=list(nbasis=23,lambda=0.1),fbetaList_l=
#                list(list(lambda=.01,nbasi=17)),hyper=NULL,
#                Cov=c('pow.ex','linear'),fitting=TRUE,
#                time=seq(-3,3,len=50),rPreIdx=TRUE,concurrent=TRUE))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gpfr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gpfrpred")
### * gpfrpred

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gpfrpred
### Title: Prediction of the Gaussian Process using functional regression
### Aliases: gpfrpred

### ** Examples

library(GPFDA)
traindata <- vector('list',20)
for(i in 1:20) traindata[[i]]=i
n <- 50
traindata <- lapply(traindata,function(i) {
  x <- seq(-3,3,len=n)
  y <- sin(x^2)-x+0.2*rnorm(n,0,3)
  x1 <- 0.5*x^3+exp(x)+rnorm(n,0,3)
  x2 <- cos(x^3)+0.2*rnorm(n,0,3)
  mat <- cbind(x,x1,x2,y)
  colnames(mat) <- c('time','x1','x2','y')
  scale <- t(c(2*(mean(y)>0.25)-1,(var(y)>3.6)*2-1,(sd(y)-sd(x)>1.4)*2-1))
  i <- list(mat,scale)
})
n <- 800 #test input
x <- seq(-3,3,len=n)
y <- sin(x^2)-x+0.2*rnorm(n,0,3)
x1 <- 0.5*x^3+exp(x)+rnorm(n,0,3)
x2 <- cos(x^3)+0.2*rnorm(n,0,3)
mat <- cbind(x,x1,x2,y)
colnames(mat) <- c('time','x1','x2','y')
scale <- t(c(2*(mean(y)>0.25)-1,(var(y)>3.6)*2-1,(sd(y)-sd(x)>1.4)*2-1))
n <- 100 # test new points
xt <- seq(1,3,len=n)
yt <- sin(xt^2)-xt+0.2*rnorm(n,0,3)
xt1 <- 0.5*xt^3+exp(xt)+rnorm(n,0,3)
xt2 <- cos(xt^3)+0.2*rnorm(n,0,3)
mat_t <- cbind(xt,xt1,xt2)
colnames(mat_t) <- c('time','xt1','xt2')
td <- list(mat,scale,mat_t)
lx=do.call('rbind',lapply(traindata,function(i)i[[2]]))
fx1=do.call('rbind',lapply(traindata,function(i)i[[1]][,2]))
fx2=do.call('rbind',lapply(traindata,function(i)i[[1]][,3]))
fy1=do.call('rbind',lapply(traindata,function(i)i[[1]][,4]))
time_old=traindata[[1]][[1]][,1]
pfx=td[[1]][,c(2,3)]
pfy=td[[1]][,4]
ptime=td[[1]][,1]
time_new=td[[3]][,1]
tfx=td[[3]][,c(2,3)]
tx=td[[2]]
## comment out because running time is a bit long
# system.time(a1<-gpfr(response=(fy1),lReg=lx,fReg=NULL,gpReg=list(fx1,fx2),
# fyList=list(nbasis=23,lambda=0.1),fbetaList_l=list(list(lambda=100,
# nbasi=17)),hyper=NULL,Cov=c('pow.ex','linear'),fitting=TRUE,
# time=seq(-3,3,len=50),rPreIdx=TRUE,concurrent=TRUE))

# type I prediction
# system.time(b1<-gpfrpred(a1,TestData=(tfx),NewTime=time_new,lReg=tx,
# fReg=NULL,gpReg=list('response'=(pfy),'input'=(pfx),'time'=ptime)))

# type II prediction
# system.time(b2<-gpfrpred(a1,TestData=(tfx),NewTime=time_new,lReg=tx,
# fReg=NULL,gpReg=NULL))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gpfrpred", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat2fd")
### * mat2fd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat2fd
### Title: Create an fd object from a matrix
### Aliases: mat2fd

### ** Examples

ry=rnorm(20,sd=10)
y1=matrix(0,ncol=100,nrow=20)
for(i in 1:20)  y1[i,]=sin(seq(-1,pi,len=100))*ry[i]

y1fd=mat2fd(y1)
y1fd=mat2fd(y1,list(lambda=1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat2fd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

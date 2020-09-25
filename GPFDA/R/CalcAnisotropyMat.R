 

# Calculating anisotropy matrix given values for the unconstrained parameters

# Q=2 setting
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

# Q=3 setting
CalcA_Q3 <- function(theta){
  
  l1 <- c(exp(theta[1]),0,0)
  l2 <- c(exp(theta[2]),pi*exp(theta[4])/(1+exp(theta[4])),0)
  l3 <- c(exp(theta[3]),pi*exp(theta[5])/(1+exp(theta[5])),pi*exp(theta[6])/(1+exp(theta[6])))
  
  L1 <- c(l1[1]*cos(l1[2]), 0, 0)
  L2 <- c(l2[1]*cos(l2[2]), l2[1]*sin(l2[2])*cos(l2[3]), 0)
  L3 <- c(l3[1]*cos(l3[2]), l3[1]*sin(l3[2])*cos(l3[3]), l3[1]*sin(l3[2])*sin(l3[3]))
  
  myL <- cbind(L1,L2,L3, deparse.level = 0)
  A <- crossprod(myL)
  diag(A) <- diag(A) + 1e-4
  
  return(A)
} 

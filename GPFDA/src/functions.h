#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// List of square matrices for NSGP covariance function
List calcScaleDistMats(List A_List, NumericMatrix coords);

// List of non-square matrices for NSGP covariance function
List calcScaleDistMatsAsym(List A_List, List Astar_List,
                             NumericMatrix coords, NumericMatrix coordsStar);

// Square MGP covariance matrix
arma::mat KCGP(NumericMatrix X, NumericVector idx, NumericVector va0s, 
                   NumericVector va1s, List A0s, List A1s, NumericVector sig);

// Non-square MGP covariance matrix
arma::mat KCGPnm(NumericMatrix X, NumericMatrix Xp, NumericVector idx, 
     	              	NumericVector idx_new, NumericVector va0s, NumericVector va1s, 
           	      	List A0s, List A1s, NumericVector sig);

// Square covariance matrix using Matern class
arma::mat CovMaternCpp(NumericMatrix input, NumericMatrix inputNew, NumericVector cc, 
                       NumericMatrix A, NumericVector nu);

// Non-square covariance matrix using Matern class
arma::mat CovMaternCppSq(NumericMatrix input, NumericVector cc, NumericMatrix A, NumericVector nu);


// Generalised distances (square and non-square matrices)
arma::mat distMat(NumericMatrix input, NumericMatrix inputNew, NumericMatrix A, NumericVector power);
arma::mat distMatSq(NumericMatrix input, NumericMatrix A, NumericVector power);
arma::mat distMatLinear(NumericMatrix input, NumericMatrix inputNew, NumericMatrix A);
arma::mat distMatLinearSq(NumericMatrix input, NumericMatrix A);


#endif

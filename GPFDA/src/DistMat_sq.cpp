#include <math.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat DistMat_sq(NumericMatrix X, NumericMatrix A, NumericVector power) {
  
  arma::mat x = Rcpp::as<arma::mat>(X);
  int nn = x.n_rows;
  
  double powerr = Rcpp::as<double>(power);
  arma::mat AA = Rcpp::as<arma::mat>(A);
  
  arma::mat M = arma::zeros(nn,nn);
  
  for(int i=0;i<nn;i++){
    for(int j=0;j<(i+1);j++){
      arma::mat xi = x.row(i);
      arma::mat xj = x.row(j);
      
      arma::mat xsep = abs(xi - xj);
      xsep = pow(xsep, powerr/2);
      double AbsDist = as_scalar(xsep*AA*xsep.t());
      M(i,j) = AbsDist;
      M(j,i) = M(i,j);
    }
  }
  
  return(M);
}

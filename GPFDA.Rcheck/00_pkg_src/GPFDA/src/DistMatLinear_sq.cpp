#include <math.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat DistMatLinear_sq(NumericMatrix X, NumericMatrix A) {
  
  arma::mat x = Rcpp::as<arma::mat>(X);
  int nn = x.n_rows;
  arma::mat AA = Rcpp::as<arma::mat>(A);
  
  arma::mat M = arma::zeros(nn,nn);
  
  for(int i=0;i<nn;i++){
    for(int j=0;j<(i+1);j++){
      arma::mat xi = x.row(i);
      arma::mat xj = x.row(j);
      M(i,j) = as_scalar(xi*AA*xj.t());
      M(j,i) = M(i,j);
    }
  }
  
  return(M);
}

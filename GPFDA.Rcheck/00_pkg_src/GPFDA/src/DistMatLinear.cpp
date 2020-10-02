#include <math.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat DistMatLinear(NumericMatrix X, NumericMatrix Xnew, NumericMatrix A) {
  
  arma::mat x = Rcpp::as<arma::mat>(X);
  arma::mat xNew = Rcpp::as<arma::mat>(Xnew);
  int nn = x.n_rows;
  int nNew = xNew.n_rows;
  arma::mat AA = Rcpp::as<arma::mat>(A);
  
  arma::mat M = arma::zeros(nn,nNew);
  
  for(int i=0;i<nn;i++){
    for(int j=0;j<nNew;j++){
      arma::mat xi = x.row(i);
      arma::mat xj = xNew.row(j);
      M(i,j) = as_scalar(xi*AA*xj.t());
    }
  }
  
  return(M);
}

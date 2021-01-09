#include <math.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' @title Calculate generalised distances
//' @description Calculate the generalised distance between vectors t and t'
//' using an anisotropy matrix A.
//' 
//' \itemize{
//' \item \code{distMat} and \code{distMatSq} calculate:
//' \deqn{ [(t - t')^{p/2}]^T A  (t - t')^{p/2}  }
//' \item \code{distMatLinear} and \code{distMatLinearSq} calculate:
//' \deqn{ t^T A t' }
//' }
//' 
//' 
//' @name distanceMatrix
//' @param input Vector of the input coordinate t
//' @param inputNew Vector of the input coordinate t'
//' @param A Anisotropy matrix A
//' @param power Power value p
//' @details The \code{distMatSq} and \code{distMatLinearSq} functions are 
//' used when input vectors t and t' are identical, returning a symmetric matrix. \cr \cr
//' When \code{distMat} and \code{distMatSq} functions are used in 
//' powered exponential kernels, power=1 gives the exponential kernel and 
//' power=2 gives the squared exponential one. \cr \cr
//' \code{distMatLinear} and \code{distMatLinearSq} functions are used in the 
//' linear covariance kernel.
//' @return A matrix
// NULL


//' @rdname distanceMatrix
//' @export
// [[Rcpp::export]]
arma::mat distMat(NumericMatrix input, NumericMatrix inputNew, NumericMatrix A, NumericVector power) {
  
  arma::mat x = Rcpp::as<arma::mat>(input);
  arma::mat xNew = Rcpp::as<arma::mat>(inputNew);
  int nn = x.n_rows;
  int nNew = xNew.n_rows;
  
  double powerr = Rcpp::as<double>(power);
  arma::mat AA = Rcpp::as<arma::mat>(A);
  
  arma::mat M = arma::zeros(nn,nNew);
  
  for(int i=0;i<nn;i++){
    for(int j=0;j<nNew;j++){
      arma::mat xi = x.row(i);
      arma::mat xj = xNew.row(j);
      
      arma::mat xsep = abs(xi - xj);
      xsep = pow(xsep, powerr/2);
      double AbsDist = as_scalar(xsep*AA*xsep.t());
      M(i,j) = AbsDist;
    }
  }
  
  return(M);
}


//' @rdname distanceMatrix
//' @export
// [[Rcpp::export]]
arma::mat distMatSq(NumericMatrix input, NumericMatrix A, NumericVector power) {
  
  arma::mat x = Rcpp::as<arma::mat>(input);
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


//' @rdname distanceMatrix
//' @export
// [[Rcpp::export]]
arma::mat distMatLinear(NumericMatrix input, NumericMatrix inputNew, NumericMatrix A) {
  
  arma::mat x = Rcpp::as<arma::mat>(input);
  arma::mat xNew = Rcpp::as<arma::mat>(inputNew);
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

//' @rdname distanceMatrix
//' @export
// [[Rcpp::export]]
arma::mat distMatLinearSq(NumericMatrix input, NumericMatrix A) {
  
  arma::mat x = Rcpp::as<arma::mat>(input);
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


#include <math.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat CovMaternCppSq(NumericMatrix input, NumericVector cc, NumericMatrix A, NumericVector nu) {
  
  arma::mat x = Rcpp::as<arma::mat>(input);
  int nn = x.n_rows;
  
  double c = Rcpp::as<double>(cc);
  arma::mat AA = Rcpp::as<arma::mat>(A);
  double nuu = Rcpp::as<double>(nu);
  
  arma::mat CovMat = arma::zeros(nn,nn);
  
  for(int i=0;i<nn;i++){
    for(int j=0;j<(i+1);j++){
      
      arma::mat xi = x.row(i);
      arma::mat xj = x.row(j);
      arma::mat xsep = abs(xi - xj);
      
      double eucDist = sqrt(as_scalar(xsep*AA*xsep.t()));
      eucDist = sqrt(2*nuu)*eucDist;
      if(eucDist < 1e-10){
        CovMat(i,j) = c;
      }else{
        double besselMod = Rf_bessel_k(eucDist, nuu, 1);
        CovMat(i,j) = c*pow(eucDist, nuu)*besselMod/(tgamma(nuu)*pow(2, (nuu-1)));
      }
      
      CovMat(j,i) = CovMat(i,j);
    }
  }
  
  return(CovMat);
}

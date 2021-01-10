#include <math.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat KCGP(NumericMatrix X, NumericVector idx, NumericVector va0s, 
                   NumericVector va1s, List A0s, List A1s, NumericVector sig) {
  arma::mat x = Rcpp::as<arma::mat>(X);
  int nn = x.n_rows;
  int Q = x.n_cols;
  
  arma::vec idxvec = as<arma::vec>(idx);
  double sigma = as<double>(sig);
  arma::vec va0vec = as<arma::vec>(va0s);
  arma::vec va1vec = as<arma::vec>(va1s);
  
  Rcpp::List A0list(A0s);
  Rcpp::List A1list(A1s);
  
  arma::mat CovMat = arma::zeros(nn,nn);
  int mrow = CovMat.n_rows;
  
  double myJit = 1e-8;
  
  for(int i=0;i<mrow;i++){
    for(int j=0;j<(i+1);j++){ 
      
      arma::mat xi = x.row(i);
      arma::mat xj = x.row(j);
      arma::mat xsep = xi - xj;
      
      int idx_i = idxvec[i];
      int idx_j = idxvec[j];
      double vXidi = va0vec[idx_i-1];
      double vXidj = va0vec[idx_j-1];
      double vEta = va1vec[idx_i-1];
      
      arma::mat AXidi = as<arma::mat>(A0list[idx_i-1]);
      arma::mat AXidj = as<arma::mat>(A0list[idx_j-1]);
      arma::mat AEta = as<arma::mat>(A1list[idx_i-1]);
      arma::mat sumAXis = AXidi + AXidj;
      arma::mat Sig = AXidi*inv(sumAXis)*AXidj;
      
      if(idx_i == idx_j){
        CovMat(i,j) = (pow(M_PI,Q/2)*vXidi*vXidi/sqrt(det(AXidi)))*exp(-0.25*as_scalar(xsep*AXidi*xsep.t())) + 
          (pow(M_PI,Q/2)*vEta*vEta/sqrt(det(AEta)) )*exp(-0.25*as_scalar(xsep*AEta*xsep.t()));
        if(i==j){
          CovMat(i,j) = CovMat(i,j) + pow(sigma,2) + myJit;
        }
      }else{
        CovMat(i,j) = (pow(2*M_PI,Q/2)*vXidi*vXidj/sqrt(det(sumAXis)))*exp(-0.5*as_scalar(xsep*Sig*xsep.t()));
      }
      
      CovMat(j,i) = CovMat(i,j);
    }
  }

  return(CovMat);
}

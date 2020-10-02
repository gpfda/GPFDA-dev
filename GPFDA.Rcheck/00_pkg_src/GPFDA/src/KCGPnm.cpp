#include <math.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat KCGPnm(NumericMatrix X, NumericMatrix Xp, NumericVector idx, 
               NumericVector idx_new, NumericVector va0s, NumericVector va1s, 
               List A0s, List A1s, NumericVector sig) {

  arma::mat x = Rcpp::as<arma::mat>(X);
  arma::mat xp = Rcpp::as<arma::mat>(Xp);
  
  int Q = x.n_cols;
  int nrow = x.n_rows;
  int ncol = xp.n_rows;
  
  arma::vec idxvec = as<arma::vec>(idx);
  arma::vec idxvecstar = as<arma::vec>(idx_new);
  double sigma = as<double>(sig);
  arma::vec va0vec = as<arma::vec>(va0s);
  arma::vec va1vec = as<arma::vec>(va1s);
  
  Rcpp::List A0list(A0s);
  Rcpp::List A1list(A1s);
  
  arma::mat CovMat = arma::zeros(nrow,ncol);
  
  double myJit = 1e-8;
  
  for(int i=0;i<nrow;i++){
    for(int j=0;j<ncol;j++){
      
      arma::mat xi = x.row(i);
      arma::mat xj = xp.row(j);
      arma::mat xsep = xi - xj;
      
      int idx_i = idxvec[i];
      int idx_j = idxvecstar[j];
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
      
    }
  }

  return(CovMat);
}

#include <math.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' @title Calculate matrices for NSGP covariance function
//' @description Calculates matrices 'ScaleMat' and 'DistMat', which are used to 
//' obtain NSGP covariance matrices
//' @param A_List List of anisotropy matrices
//' @param coords Matrix of input coordinates (covariates)
//' 
//' @return A list of ScaleMat and DistMat matrices
//' 
//' @export
//' @examples
//' ## See examples in vignette:
//' # vignette("nsgpr", package = "GPFDA")
// [[Rcpp::export]]
List calcScaleDistMats(List A_List, NumericMatrix coords) {
  
  arma::mat coordins = Rcpp::as<arma::mat>(coords);
  Rcpp::List AList(A_List);
  int N = AList.size();
  
  arma::mat ScaleMat = arma::zeros(N,N);
  arma::mat DistMat = arma::zeros(N,N);
  
  for(int i=0;i<N;i++){
    
    arma::mat Kerneli = inv(as<arma::mat>(AList[i]));
    double det_i = det(Kerneli);
    
    ScaleMat(i,i) = 1;
    DistMat(i,i) = 0;
    
    
    if(i<N){
      for(int j=0; j<(i+1); j++){
          
        arma::mat Kernelj = inv(as<arma::mat>(AList[j]));
        double det_j = det(Kernelj);
        
        arma::mat avg_ij = 0.5 * (Kerneli + Kernelj);
        double det_ij = det(avg_ij);
        
        ScaleMat(i,j) = sqrt( sqrt(det_i*det_j) / det_ij );
        DistMat(i,j) = sqrt( as_scalar((coordins.row(i)-coordins.row(j)) * inv(avg_ij) * (coordins.row(i)-coordins.row(j) ).t()));
        
        ScaleMat(j,i) = ScaleMat(i,j);
        DistMat(j,i) = DistMat(i,j);
        
      }
    }
  }
  
  
  List ret;
  ret["Scale.mat"] = ScaleMat;
  ret["Dist.mat"] = DistMat;
  
  return(ret);
}

#include <math.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List calcScaleDistMatsAsym(List A_List, List Astar_List,
                                NumericMatrix coords, NumericMatrix coordsStar) {

  arma::mat coordins = Rcpp::as<arma::mat>(coords);
  arma::mat coordinsStar = Rcpp::as<arma::mat>(coordsStar);
  Rcpp::List AList(A_List);
  Rcpp::List AstarList(Astar_List);
  
  int N = AList.size();
  int Nstar = AstarList.size();
  
  arma::mat ScaleMat = arma::zeros(N,Nstar);
  arma::mat DistMat = arma::zeros(N,Nstar);
  
  for(int i=0;i<N;i++){
    
    arma::mat Kerneli = inv(as<arma::mat>(AList[i]));
    double det_i = det(Kerneli);
    
    for(int j=0; j<Nstar; j++){
      
      arma::mat Kernelj = inv(as<arma::mat>(AstarList[j]));
      double det_j = det(Kernelj);
      
      arma::mat avg_ij = 0.5 * (Kerneli + Kernelj);
      double det_ij = det(avg_ij);
      
      ScaleMat(i,j) = sqrt( sqrt(det_i*det_j) / det_ij );
      DistMat(i,j) = sqrt( as_scalar((coordins.row(i)-coordinsStar.row(j)) * inv(avg_ij) * (coordins.row(i)-coordinsStar.row(j) ).t()));
      
    }
    
  }
  
  
  List ret;
  ret["Scale.mat"] = ScaleMat;
  ret["Dist.mat"] = DistMat;
  
  return(ret);
}

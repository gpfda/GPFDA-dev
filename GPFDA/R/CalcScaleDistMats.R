
# Fast calculation of ScaleMat and DistMat for NSGP covariance functions

include='
#include <math.h>
#include <vector>
using namespace Rcpp;
using namespace arma;
'

body='
mat coordins = Rcpp::as<mat>(coords);
Rcpp::List AList(A_List);
int N = AList.size();

mat ScaleMat = zeros(N,N);
mat DistMat = zeros(N,N);

for(int i=0;i<N;i++){

  mat Kerneli = inv(as<mat>(AList[i]));
  double det_i = det(Kerneli);
  
  ScaleMat(i,i) = 1;
  DistMat(i,i) = 0;
  
  
    if(i<N){
      for(int j=0; j<(i+1); j++){
      
        mat Kernelj = inv(as<mat>(AList[j]));
        double det_j = det(Kernelj);
        
        mat avg_ij = 0.5 * (Kerneli + Kernelj);
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

return(wrap(ret));
'

CalcScaleDistMats <- cxxfunction(signature(A_List="List",
                               coords="matrix"),
                    body,plugin='RcppArmadillo',includes=include)




body='
mat coordins = Rcpp::as<mat>(coords);
mat coordinsStar = Rcpp::as<mat>(coordsStar);
Rcpp::List AList(A_List);
Rcpp::List AstarList(Astar_List);

int N = AList.size();
int Nstar = AstarList.size();

mat ScaleMat = zeros(N,Nstar);
mat DistMat = zeros(N,Nstar);

for(int i=0;i<N;i++){

  mat Kerneli = inv(as<mat>(AList[i]));
  double det_i = det(Kerneli);
  
  ScaleMat(i,i) = 1;
  DistMat(i,i) = 0;
  
    for(int j=0; j<Nstar; j++){
    
      mat Kernelj = inv(as<mat>(AstarList[j]));
      double det_j = det(Kernelj);
      
      mat avg_ij = 0.5 * (Kerneli + Kernelj);
      double det_ij = det(avg_ij);
      
      ScaleMat(i,j) = sqrt( sqrt(det_i*det_j) / det_ij );
      DistMat(i,j) = sqrt( as_scalar((coordins.row(i)-coordinsStar.row(j)) * inv(avg_ij) * (coordins.row(i)-coordinsStar.row(j) ).t()));
    
    }
    
}


List ret;
ret["Scale.mat"] = ScaleMat;
ret["Dist.mat"] = DistMat;

return(wrap(ret));
'

CalcScaleDistMatsAsym <- cxxfunction(signature(A_List="List",
                                               Astar_List="List",
                                           coords="matrix",
                                           coordsStar="matrix"),
                                 body,plugin='RcppArmadillo',includes=include)




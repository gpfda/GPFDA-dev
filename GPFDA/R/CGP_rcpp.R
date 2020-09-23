
include='
#include <math.h>
#include <vector>
using namespace Rcpp;
using namespace arma;
'

body='
mat x = Rcpp::as<mat>(X);
int nn = x.n_rows;
int Q = x.n_cols;

vec idxvec = as<vec>(idx);
double sigma = as<double>(sig);
vec va0vec = as<vec>(va0s);
vec va1vec = as<vec>(va1s);
int N = va0vec.size();

Rcpp::List A0list(A0s);
Rcpp::List A1list(A1s);

mat CovMat = zeros(nn,nn);
int mrow = CovMat.n_rows;
int mcol = CovMat.n_cols;

double myJit = 1e-8;

for(int i=0;i<mrow;i++){
  for(int j=0;j<(i+1);j++){ 

    mat xi = x.row(i);
    mat xj = x.row(j);
    mat xsep = xi - xj;
    
    int idx_i = idxvec[i];
    int idx_j = idxvec[j];
    double vXidi = va0vec[idx_i-1];
    double vXidj = va0vec[idx_j-1];
    double vEta = va1vec[idx_i-1];
    
    mat AXidi = as<mat>(A0list[idx_i-1]);
    mat AXidj = as<mat>(A0list[idx_j-1]);
    mat AEta = as<mat>(A1list[idx_i-1]);
    mat sumAXis = AXidi + AXidj;
    mat Sig = AXidi*inv(sumAXis)*AXidj;
    
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


return(wrap(CovMat));
'


KCGP <- cxxfunction(signature(X="matrix", idx="integer", 
                              va0s="numeric",va1s="numeric",
                              A0s="List", A1s="List", sig="numeric"),
                    body,plugin='RcppArmadillo',includes=include)



include='
#include <math.h>
#include <vector>
using namespace Rcpp;
using namespace arma;
'

body='
mat x = Rcpp::as<mat>(X);
mat xp = Rcpp::as<mat>(Xp);

int Q = x.n_cols;
int nrow = x.n_rows;
int ncol = xp.n_rows;

vec idxvec = as<vec>(idx);
vec idxvecstar = as<vec>(idx_new);
double sigma = as<double>(sig);
vec va0vec = as<vec>(va0s);
vec va1vec = as<vec>(va1s);
int N = va0vec.size();

Rcpp::List A0list(A0s);
Rcpp::List A1list(A1s);

mat CovMat = zeros(nrow,ncol);

double dmax = 100;
double myJit = 1e-8;

for(int i=0;i<nrow;i++){
  for(int j=0;j<ncol;j++){
  
    mat xi = x.row(i);
    mat xj = xp.row(j);
    mat xsep = xi - xj;
    
    int idx_i = idxvec[i];
    int idx_j = idxvecstar[j];
    double vXidi = va0vec[idx_i-1];
    double vXidj = va0vec[idx_j-1];
    double vEta = va1vec[idx_i-1];
    
    mat AXidi = as<mat>(A0list[idx_i-1]);
    mat AXidj = as<mat>(A0list[idx_j-1]);
    mat AEta = as<mat>(A1list[idx_i-1]);
    mat sumAXis = AXidi + AXidj;
    mat Sig = AXidi*inv(sumAXis)*AXidj;
    
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


return(wrap(CovMat));
'

KCGPnm <- cxxfunction(signature(X="matrix", Xp="matrix", idx="integer", idx_new="integer",
                              va0s="numeric",va1s="numeric",
                              A0s="List", A1s="List", sig="numeric"),
                    body,plugin='RcppArmadillo',includes=include)


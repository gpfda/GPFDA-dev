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

double c = Rcpp::as<double>(cc);
mat AA = Rcpp::as<mat>(A);
double nuu = Rcpp::as<double>(nu);

mat CovMat = zeros(nn,nn);

for(int i=0;i<nn;i++){
  for(int j=0;j<(i+1);j++){
  
    mat xi = x.row(i);
    mat xj = x.row(j);
    mat xsep = abs(xi - xj);

    double eucDist = sqrt(as_scalar(xsep*AA*xsep.t()));
    eucDist = sqrt(2*nuu)*eucDist;
    if(eucDist < 1e-10){
      eucDist = 1e-10;
    }

    double besselMod = Rf_bessel_k(eucDist, nuu, 1);
    CovMat(i,j) = c*pow(eucDist, nuu)*besselMod/(tgamma(nuu)*pow(2, (nuu-1)));
    CovMat(j,i) = CovMat(i,j);
  }
}

return(wrap(CovMat));
'
CovMaternCpp_sq <- cxxfunction(signature(X="matrix",cc="double", A="matrix", nu="double"),
                             body,plugin='RcppArmadillo',includes=include)


# mat xsep = xi - xj; 
# 
# double eucDist = sqrt(as_scalar(xsep*AA*xsep.t()));
# eucDist = sqrt(2*nuu)*eucDist;
# if(eucDist < 1e-10){
#   eucDist = 1e-10;
# }
# double besselMod = Rf_bessel_k(eucDist, nuu, 1);
# CovMat(i,j) = c*pow(eucDist, nuu)*besselMod/(tgamma(nuu)*pow(2, (nuu-1)));


# xsep = pow(xsep, 0.5);
# double AbsDist = as_scalar(xsep*AA*xsep.t());
# AbsDist = sqrt(2*nuu)*AbsDist;
# if(AbsDist < 1e-8){
#   AbsDist = 1e-8;
# }


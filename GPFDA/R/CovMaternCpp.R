include='
#include <math.h>
#include <vector>
using namespace Rcpp;
using namespace arma;
'
body='

mat x = Rcpp::as<mat>(X);
mat xNew = Rcpp::as<mat>(Xnew);
int nn = x.n_rows;
int nNew = xNew.n_rows;
int Q = x.n_cols;

double c = Rcpp::as<double>(cc);
mat AA = Rcpp::as<mat>(A);
double nuu = Rcpp::as<double>(nu);

mat CovMat = zeros(nn,nNew);

for(int i=0;i<nn;i++){
  for(int j=0;j<nNew;j++){
  
    mat xi = x.row(i);
    mat xj = xNew.row(j);
    mat xsep = abs(xi - xj);
    double eucDist = sqrt(as_scalar(xsep*AA*xsep.t()));
    eucDist = sqrt(2*nuu)*eucDist;
    if(eucDist < 1e-10){
      eucDist = 1e-10;
    }
    double besselMod = Rf_bessel_k(eucDist, nuu, 1);
    CovMat(i,j) = c*pow(eucDist, nuu)*besselMod/(tgamma(nuu)*pow(2, (nuu-1)));
  }
}

return(wrap(CovMat));
'
CovMaternCpp <- cxxfunction(signature(X="matrix", Xnew="matrix", cc="double", A="matrix", nu="double"),
                             body,plugin='RcppArmadillo',includes=include)

# mat xsep = xi - xj; 
# double eucDist = sqrt(as_scalar(xsep*AA*xsep.t()));
# eucDist = sqrt(2*nuu)*eucDist;
# if(eucDist < 1e-10){
#   eucDist = 1e-10;
# }
# double besselMod = Rf_bessel_k(eucDist, nuu, 1);
# CovMat(i,j) = c*pow(eucDist, nuu)*besselMod/(tgamma(nuu)*pow(2, (nuu-1)));

# bessel_j (index, expr)         Bessel function, 1st kind
# bessel_y (index, expr)         Bessel function, 2nd kind
# bessel_i (index, expr)         Modified Bessel function, 1st kind
# bessel_k (index, expr)         Modified Bessel function, 2nd kind


# xsep = pow(xsep, 0.5);
# double AbsDist = as_scalar(xsep*AA*xsep.t());
# AbsDist = sqrt(2*nuu)*AbsDist;
# if(AbsDist < 1e-8){
#   AbsDist = 1e-8;
# }

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

double powerr = Rcpp::as<double>(power);
mat AA = Rcpp::as<mat>(A);

mat M = zeros(nn,nNew);

for(int i=0;i<nn;i++){
  for(int j=0;j<nNew;j++){
    mat xi = x.row(i);
    mat xj = xNew.row(j);
    mat xsep = xi - xj;
    double eucDist = as_scalar(xsep*AA*xsep.t());
    M(i,j) = pow(eucDist, powerr);
  }
}

return(wrap(M));
'
xixj_staNEW <- cxxfunction(signature(X="matrix", Xnew="matrix", A="matrix", power="double"),
                             body,plugin='RcppArmadillo',includes=include)


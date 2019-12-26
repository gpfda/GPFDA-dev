
include='
#include <math.h>
#include <vector>
using namespace Rcpp;
using namespace arma;
'
body='

mat x = Rcpp::as<mat>(X);
int nn = x.n_rows;

double powerr = Rcpp::as<double>(power);
mat AA = Rcpp::as<mat>(A);

mat M = zeros(nn,nn);

for(int i=0;i<nn;i++){
  for(int j=0;j<(i+1);j++){
    mat xi = x.row(i);
    mat xj = x.row(j);
    mat xsep = xi - xj;
    double eucDist = as_scalar(xsep*AA*xsep.t());
    M(i,j) = pow(eucDist, powerr);
    M(j,i) = M(i,j);
  }
}

return(wrap(M));
'
xixj_staNEW_sq <- cxxfunction(signature(X="matrix", A="matrix", power="double"),
                             body,plugin='RcppArmadillo',includes=include)


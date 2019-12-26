
include='
#include <math.h>
#include <vector>
using namespace Rcpp;
using namespace arma;
'
body='
mat x = Rcpp::as<mat>(X);
int nn = x.n_rows;
mat AA = Rcpp::as<mat>(A);

mat M = zeros(nn,nn);

for(int i=0;i<nn;i++){
  for(int j=0;j<(i+1);j++){
    mat xi = x.row(i);
    mat xj = x.row(j);
    M(i,j) = as_scalar(xi*AA*xj.t());
    M(j,i) = M(i,j);
  }
}

return(wrap(M));
'
xixjNEW_sq <- cxxfunction(signature(X="matrix", A="matrix"),
                             body,plugin='RcppArmadillo',includes=include)


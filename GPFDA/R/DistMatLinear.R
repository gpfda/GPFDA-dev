
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
mat AA = Rcpp::as<mat>(A);

mat M = zeros(nn,nNew);

for(int i=0;i<nn;i++){
  for(int j=0;j<nNew;j++){
    mat xi = x.row(i);
    mat xj = xNew.row(j);
    M(i,j) = as_scalar(xi*AA*xj.t());
  }
}

return(wrap(M));
'
DistMatLinear <- cxxfunction(signature(X="matrix", Xnew="matrix", A="matrix"),
                             body,plugin='RcppArmadillo',includes=include)


# previously called xixjNEW
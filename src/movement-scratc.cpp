#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
NumericVector move(NumericVector x, Eigen::MatrixXd m) {
 
 VectorXd tmp(as<VectorXd>(x));
  
  
  tmp = m * tmp;
  
  SEXP out = Rcpp::wrap(tmp); // convert from eigen to Rcpp
  
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

x = c(1,10,3)

m = matrix(rep(c(1,0,0,0,.5,.5, 0,0,1)), nrow = 3, ncol = 3)

z = m %*% (x)

a = move(x, m)
*/

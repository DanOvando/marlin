#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::NumericMatrix sample_problem() {

  Eigen::MatrixXd x(2, 2), y(2, 2);
  x << 1,1,2,2;
  y << 3,3,4,4;

  // do some eigen matrix multiplication
  Eigen::MatrixXd z =  x * y;

  // what I'd like to be able to do somehow:
  // store the results of the eigen object z in
  // a NumericMatrix w
  // w = z;
  SEXP s = Rcpp::wrap(z);
  Rcpp::NumericMatrix w(s);

  return w;
}

/*** R
sample_problem()
*/

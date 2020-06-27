#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Eigen::MatrixXd move_matrix_eigen(Eigen::MatrixXd X, Eigen::MatrixXd y, int its) {

  // arma::rowvec z(y.n_cols);

  Eigen::MatrixXd z(y.rows(), y.cols());

  Eigen::MatrixXd mmat(X.rows(), X.cols());

  mmat = X.transpose();

  z = y;

  for (int i = 0; i < its; i++){


    z =  mmat * z;


  }


  return z;
} // close movement model



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
/*** R
m <- matrix(1:6, 3,3)

m <- m / rowSums(m)

v <- matrix(1:12, nrow = 3, ncol = 4)

d <- crossprod(m,v)

m %*% v

d <- crossprod(m,d)

d
a = move_matrix_eigen(m, d, 10)


a
*/

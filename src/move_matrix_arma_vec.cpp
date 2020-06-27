// #include <RcppArmadillo.h>
// // [[Rcpp::depends(RcppArmadillo)]]
// using namespace Rcpp;
// using namespace arma;
//
// // [[Rcpp::export]]
// arma::rowvec move_matrix_arma(const arma::mat& X, const arma::rowvec& y, int its) {
//
//   arma::rowvec z(y.n_cols);
//
//   z = y;
//
//   for (int i = 0; i < its; i++){
//
//
//     z = z * X;
//
//
//   }
//
//
//   return z;
// }
//
//
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically
// // run after the compilation.
// //
// /*** R
//
// m <- matrix(1:6, 3,3)
//
// m <- m / rowSums(m)
//
// v <- 1:3
//
// crossprod(v,m)
//
// a = move_matrix_arma(m, v, 10)
//
// sum(a)
// */

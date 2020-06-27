//
// #include <RcppArmadillo.h>
// // [[Rcpp::depends(RcppArmadillo)]]
// using namespace Rcpp;
// using namespace arma;
//
// // [[Rcpp::export]]
// arma::mat move_matrix_arma(const arma::mat x, arma::mat y, int its) {
//
//   arma::mat z(y.n_rows, y.n_cols);
//
//   z = y;
//
//   arma::mat mmat(x.n_rows, x.n_cols);
//
//   mmat = x.t();
//
//   for (int i = 0; i < its; i++){
//
//     z = mmat * z;
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
//   /*** R
//
// m <- matrix(1:6, 3,3)
//
// m <- m / rowSums(m)
//
// v <- matrix(1:12, nrow = 3, ncol = 4)
//
// d <- crossprod(m,v)
//
// m %*% v
//
// d <- crossprod(m,d)
//
// d
// a = move_matrix_arma(m, d, 10)
//
// sum(a)
// */

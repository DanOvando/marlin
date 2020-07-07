#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd sample_problem(Eigen::MatrixXd x, Eigen::MatrixXd y) {


  // do some eigen matrix multiplication
  Eigen::MatrixXd z =  x * y;

  // what I'd like to be able to do somehow:
  // store the results of the eigen object z in
  // a NumericMatrix w
  // w = z;
  // SEXP s = Rcpp::wrap(z);
  // Rcpp::NumericMatrix w(s);

  return z;
}

/*** R

x = matrix(c(1,1,2,2), nrow = 2, ncol = 2, byrow = TRUE)

y = matrix(c(3,3,4,4), nrow = 2, ncol = 2, byrow = TRUE)

sample_problem(x,y)



x <- fauna$bigeye$move_mat

# x <- x / rowSums(x)

y <- fauna$bigeye$unfished$n_p_a

b = x %*% y

d = t(y) %*% (x / rowSums(x))


y[2,] = 0

sum(y)
a = sample_problem(x, y)
sum(a)
sum(x[1,] * y[,1])

*/

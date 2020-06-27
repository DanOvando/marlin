// #include <RcppArmadillo.h>
// //[[Rcpp::depends(RcppArmadillo)]]
// using namespace Rcpp;
// using namespace arma;
// // [[Rcpp::export]]
//
// List fish_model(const arma::rowvec length_at_age,
//                 const arma::vec weight_at_age,
//                 const arma::vec maturity_at_age,
//                 arma::mat movement,
//                 arma::mat initpop,
//                 const int patches,
//                 const int sim_years,
//                 const int burn_years,
//                 const double steepness,
//                 const double r0,
//                 const double m
// ) {
//
// // this only needs to run for one year! you'll hopefully pass a matrix out of here, and put it in an array...
// // but we'll see
//   int years = sim_years + burn_years; // total number of years to run
//
//   int n_ages = length_at_age.n_cols;
//
//   arma::mat mmat(movement.n_rows, movement.n_cols);
//
//   mmat = movement.t();
//
//   arma::mat n_p_a(patches, n_ages); // numbers at age over time
//
//   // n_p_a = ones(patches, n_ages);
//
//   n_p_a = initpop;
//
//   // arma::cube n_p_a_y(patches, n_ages, years);
//
//   // test the idea here
//
//   // n_p_a_y.slice(0) = n_p_a;
//
//
//   for (int i = 0; i < years; i++){
//
//
//     n_p_a =  mmat * n_p_a;
//
//
//   }
//
//   return Rcpp::List::create(
//     Rcpp::Named("n_ages") = n_ages,
//     Rcpp::Named("n_p_a") = n_p_a
//   );
//
//
// } // close popmodel

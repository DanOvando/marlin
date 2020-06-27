#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::export]]
Rcpp::NumericMatrix fish_model(
    const Eigen::VectorXd length_at_age,
    const Eigen::VectorXd weight_at_age,
    const Eigen::VectorXd maturity_at_age,
    const Eigen::MatrixXd movement,
    Eigen::MatrixXd n_p_a,
    const int patches,
    const int sim_steps,
    const int burn_steps,
    const double steepness,
    const double r0,
    const double m){

  int steps = sim_steps + burn_steps; // total number of years to run

  // int n_ages = length_at_age.rows();

  // Eigen::MatrixXd mmat(movement.rows(), movement.cols());

  // Eigen::MatrixXd n_p_a(patches, n_ages);

  // n_p_a = initpop;

  // mmat = movement.transpose();

  for (int s = 0; s < steps; s++){

    // Rcpp::Rcout << s << std::endl;

    n_p_a =  movement * n_p_a;


  }


  SEXP tmp = Rcpp::wrap(n_p_a);

  Rcpp::NumericMatrix w(tmp);


  return(w);


  } // close fish model

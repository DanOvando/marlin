#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
List sim_fish_pop(
    const NumericVector length_at_age,
    const NumericVector weight_at_age,
    const NumericVector maturity_at_age,
    const Eigen::MatrixXd movement,
    Rcpp::NumericMatrix n_p_a,
    const int patches,
    const int sim_steps,
    const  int burn_steps,
    const  double steepness,
    const  double r0,
    const  double ssb0,
    const  double m){

  int steps = sim_steps; // total number of years to run

  int ages = length_at_age.length();

  MatrixXd  tmpmat(as<MatrixXd>(n_p_a));

  // NumericMatrix n_p_a(patches, ages); // biomass at age over time

  NumericMatrix b_p_a(patches, ages); // spawning stock biomass at age over time

  NumericMatrix ssb_p_a(patches, ages); // spawning stock biomass at age over time

  NumericMatrix c_a(patches, ages); // catch at age over time

  NumericVector rec_s(steps); // catch at age over time

  //////////////////// move ////////////////////////

  // for (int s = 0; s < steps; s++){

    // Rcpp::Rcout << s << std::endl;

    tmpmat =  movement * tmpmat; // matrix multiplication of numbers at age by movement matrix


  // }

  SEXP tmp = Rcpp::wrap(tmpmat); // convert from eigen to Rcpp

  n_p_a = tmp;

  //////////////////// grow ////////////////////////

  NumericVector plus_group = n_p_a(_,ages - 1) * exp(-m);

  for (int a = 1; a < ages; a++){

    n_p_a(_,a) =  n_p_a(_,a - 1) * exp(-m);

  }

  n_p_a(_, ages - 1) = n_p_a(_, ages - 1) + plus_group;

  for (int p = 0;p < patches; p++){

    b_p_a(p,_) =  n_p_a(p,_) * weight_at_age;

    ssb_p_a(p,_) =  b_p_a(p,_) * maturity_at_age;

  }

  //////////////////// spawn / recruit ////////////////////////

    if (ssb0 == -999){

      n_p_a(_,0) = rep(r0 / patches, patches);

    } else {

      // global density dependence, distribute recruits evenly
      double ssb = sum(ssb_p_a);

      n_p_a(_,0) = rep(((0.8 * r0 * steepness * ssb) / (0.2 * ssb0 * (1 - steepness) + (steepness - 0.2) * ssb)) / patches, patches);

    }

  //////////////////// process results ////////////////////////

  return Rcpp::List::create(
    Rcpp::Named("n_p_a") = n_p_a,
    Rcpp::Named("b_p_a") = b_p_a,
    Rcpp::Named("ssb_p_a") = ssb_p_a);
} // close fish model

// [[Rcpp::export]]
List est_ssb0( NumericVector length_at_age,
                 NumericVector weight_at_age,
                 NumericVector maturity_at_age,
                Eigen::MatrixXd movement,
                 Rcpp::NumericMatrix n_p_a,
                  int patches,
                  int sim_steps,
                  int burn_steps,
                  double steepness,
                 double r0,
                  double ssb0,
                  double m) {

// you can definitely build this into the main function
  Rcpp::List tmp(1);

NumericMatrix tmp_n_p_a = n_p_a;

for (int b = 0; b < burn_steps; b++){


  tmp = sim_fish_pop(
    length_at_age,
    weight_at_age,
    maturity_at_age,
    movement,
    tmp_n_p_a,
    patches,
    sim_steps,
    burn_steps,
    steepness,
    r0,
    ssb0 = -999,
    m);


  NumericMatrix tmp_n_p_a = tmp["n_p_a"];

}

return(tmp);

} // close est_ssb0





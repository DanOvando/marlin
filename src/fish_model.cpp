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
    const Rcpp::NumericMatrix last_n_p_a,
    const int patches,
    const int burn_steps,
    const double steepness,
    const double r0,
    double ssb0,
    const double m,
    bool tune_unfished){


  int ages = length_at_age.length();

  MatrixXd  tmpmat(as<MatrixXd>(last_n_p_a));

  NumericMatrix n_p_a(patches, ages); // biomass at age over time

  NumericMatrix b_p_a(patches, ages); // spawning stock biomass at age over time

  NumericMatrix ssb_p_a(patches, ages); // spawning stock biomass at age over time

  NumericMatrix c_a(patches, ages); // catch at age over time

  //////////////////// tune things ////////////////////////
  NumericMatrix tmp_n_p_a = clone(last_n_p_a);

  Rcpp::List tmppop;

  // Rcpp::Rcout<< ssb0 << std::endl;

  // Rcpp::Rcout<< Rcpp::traits::is_na<REALSXP>(ssb0) << std::endl;

  if (Rcpp::traits::is_na<REALSXP>(ssb0) == 1){


    for (int b = 0; b < burn_steps; b++){

      // these HAVE to be in the same order
      // as function call: get lots of warnings
      // if you try and name them
      // Rcpp::List tmp;

      // tmppop = b;

      // Rcpp::Rcout << tmp_n_p_a(1,10) << std::endl;

      tmppop = sim_fish_pop(
        length_at_age,
        weight_at_age,
        maturity_at_age,
        movement,
        tmp_n_p_a,
        patches,
        0,
        steepness,
        r0,
        -999,
        m,
        1);

      tmp_n_p_a = wrap(tmppop["n_p_a"]);

     // NumericMatrix tmp_n_p_a = (tmppop["n_p_a"]);

      // Rcpp::List tmppop;

    }

    NumericMatrix tmp_ssb_p_a = tmppop["ssb_p_a"];

    ssb0 = sum(tmp_ssb_p_a);

    // return(tmp);


  }



  //////////////////// move ////////////////////////

  // for (int s = 0; s < steps; s++){

    // Rcpp::Rcout << s << std::endl;

    tmpmat =  movement * tmpmat; // matrix multiplication of numbers at age by movement matrix


  // }

  SEXP tmp = Rcpp::wrap(tmpmat); // convert from eigen to Rcpp

  n_p_a = tmp;

  //////////////////// grow ////////////////////////

  // NumericMatrix last_n_p_a = n_p_a;

  NumericVector plus_group = n_p_a(_,ages - 1) * exp(-m);

  // Rcpp::Rcout << "last is" << last_n_p_a(0,0) * exp(-m) << std::endl;

  // Rcpp::Rcout << "before is" << n_p_a(0,1)<< std::endl;

  // Rcpp::Rcout << "exp is" <<exp(-m)<< std::endl;

  // std::cout << "exp is" << last_n_p_a(_,11)* exp(-m)<< std::endl;

  for (int a = 1; a < ages; a++){

    // NumericMatrix::Column wtf = last_n_p_a(_,a - 1);

    // Rcpp::Rcout << "a is" << a - 1<< std::endl;

    n_p_a(_,a) =  last_n_p_a(_,a - 1) * exp(-m);


  }

  // Rcpp::Rcout << "new is" << n_p_a(0,1)<< std::endl;

  n_p_a(_, ages - 1) = n_p_a(_, ages - 1) + plus_group;

  // n_p_a(1, _) = rep(999, ages);

  // weight_at_age = weight_at_age + 999;


  for (int p = 0;p < patches; p++){

    b_p_a(p,_) =  n_p_a(p,_) * weight_at_age;

    ssb_p_a(p,_) =  b_p_a(p,_) * maturity_at_age;

  }

  //////////////////// spawn / recruit ////////////////////////

    if (tune_unfished == 1){

      n_p_a(_,0) = rep(r0 / patches, patches);

    } else {

      // global density dependence, distribute recruits evenly
      double ssb = sum(ssb_p_a);

      // Rcpp::Rcout << ssb << std::endl;

      n_p_a(_,0) = rep(((0.8 * r0 * steepness * ssb) / (0.2 * ssb0 * (1 - steepness) + (steepness - 0.2) * ssb)) / patches, patches);

    }

  //////////////////// process results ////////////////////////

  return Rcpp::List::create(
    Rcpp::Named("n_p_a") = n_p_a,
    Rcpp::Named("b_p_a") = b_p_a,
    Rcpp::Named("ssb_p_a") = ssb_p_a,
    Rcpp::Named("ssb0") = ssb0,
    Rcpp::Named("tmppop") = tmppop);
} // close fish model






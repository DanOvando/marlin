#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List popmodel(NumericVector length_at_age,
              NumericVector weight_at_age,
              NumericVector maturity_at_age,
              NumericVector selectivity_at_age,
              NumericVector rec_devs,
              NumericVector age_vector,
              int sim_years,
              int n_ages,
              int burn_years,
              int rec_form,
              double steepness,
              double r0,
              double m,
              double f
) {


  int years = sim_years + burn_years; // total number of years to run


  NumericMatrix n_t(years, n_ages); // numbers at age over time

  NumericMatrix b_t(years, n_ages); // biomass at age over time

  NumericMatrix ssb_t(years, n_ages); // spawning stock biomass at age over time

  NumericMatrix c_t(years, n_ages); // catch at age over time

  NumericVector rec_t(years); // catch at age over time

  // initilize population

  n_t(0,_) = r0 * exp(-m * age_vector);

  n_t(0,n_ages - 1) = n_t(0,n_ages - 2)* exp(-m) / (1 - exp(-m));

  b_t(0,_) =  n_t(0,_) * weight_at_age;

  ssb_t(0,_) =  b_t(0,_) * maturity_at_age;

  rec_t(0) = r0 * rec_devs(0);

  // burn in population to get average ssb0 in case you allow for recruitment deviates
  for (int t = 1; t < burn_years; t++){

    n_t(t,0) = rec_t(t - 1);

    for (int a = 1; a < n_ages; a++){

      n_t(t,a) = n_t(t - 1, a - 1) * exp(-m);

    }

    n_t(t, n_ages -1) = n_t(t, n_ages -1) + n_t(t - 1, n_ages -1) * exp(-m);

    b_t(t,_) =  n_t(t,_) * weight_at_age;

    ssb_t(t,_) =  b_t(t,_) * maturity_at_age;

    rec_t(t) = r0 * rec_devs(t);

  }

  // calcualte mean ssb0

  NumericMatrix burn_ssb = ssb_t(Range(0,fmax(0,burn_years - 1)),_);

  double ssb0 = mean(rowSums(burn_ssb)); // ssb0, allowing for potential recruitment variation around r0

  // project forward now with catches

  int startyear = fmax(burn_years,1);

  for (int t = (startyear); t < years; t++){

    n_t(t,0) = rec_t(t - 1); // recruit

    for (int a = 1; a < n_ages; a++){ //grow and die

      n_t(t,a) = n_t(t - 1, a - 1) * exp(-(m + f * selectivity_at_age(a - 1)));

      // c_t(t - 1, a - 1) =  (f * selectivity_at_age(a - 1)) / (f * selectivity_at_age(a - 1) + m) * n_t(t - 1, a - 1) * (1 - exp(-(m + f * selectivity_at_age(a - 1)))) * weight_at_age(a - 1);

    } // close grow and die

    c_t(t - 1,_) = (f * selectivity_at_age) / (f * selectivity_at_age + m) * n_t(t-1, _) * (1 - exp(-(m + f * selectivity_at_age))) * weight_at_age;

    n_t(t, n_ages -1) = n_t(t, n_ages -1) + n_t(t - 1, n_ages -1) * exp(-(m + f * selectivity_at_age(n_ages -1)));

    // add in plus group

    b_t(t,_) =  n_t(t,_) * weight_at_age;

    ssb_t(t,_) =  b_t(t,_) * maturity_at_age;

    double tempssb = sum(ssb_t(t,_));

    if (rec_form == 1){ // calculate recruitment

      rec_t(t) = ((0.8 * r0 * steepness * tempssb) / (0.2 * ssb0 * (1 - steepness) + (steepness - 0.2) * tempssb)) * rec_devs(t);

    } else {

      // double beta = log(5 * steepness) / (0.8 * ssb_max);
      //
      // alpha = log(r_max / ssb_max) + 0.8 * beta * ssb_max;
      //
      // rec_t(t) = tempssb .* exp(alpha - beta * ssb); //.* exp(alpha - beta .* ssb);
      //
    }
  }

  // fill in final catches
  c_t(years - 1,_) = (f * selectivity_at_age) / (f * selectivity_at_age + m) * n_t(years -1, _) * (1 - exp(-(m + f * selectivity_at_age))) * weight_at_age;


  return Rcpp::List::create(
    Rcpp::Named("n_t") = n_t,
    Rcpp::Named("b_t") = b_t,
    Rcpp::Named("ssb_t") = ssb_t,
    Rcpp::Named("c_t") = c_t,
    Rcpp::Named("burn_ssb") = burn_ssb,
    Rcpp::Named("ssb0") = ssb0,
    Rcpp::Named("rec_t") = rec_t);
} // close popmodel

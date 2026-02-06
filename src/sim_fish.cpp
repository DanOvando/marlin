#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
NumericVector cpp_seq(int steps, double step_size) {
  NumericVector y(steps);
  y(0) = 1;

  for (int i = 1; i < y.length(); i++) {
    y(i) = y(i - 1) + step_size;
  }

  return y;
}


// [[Rcpp::export]]
List sim_fish(
    const NumericVector length_at_age,
    const NumericVector weight_at_age,
    const NumericVector fec_at_age,
    const NumericVector maturity_at_age,
    bool semelparous, // introduce semelparous behavior
    const NumericMatrix f_p_a, // fishing mortality by patch and age
    const List movement_matrix,
    const List movement_seasons,
    Eigen::SparseMatrix<double> recruit_movement_matrix,
    Rcpp::NumericMatrix last_n_p_a, // last numbers by patch and age
    const int patches,
    const int burn_steps, // number of burn steps if burn is in effect
    const double time_step, // the time step in portions of a year
    int season, // what season you're currently in
    const double steepness,
    const NumericVector r0s,
    double ssb0, // unfished spawning stock biomass across entire system
    NumericVector ssb0_p, // unfished spawning stock biomass in each patch
    const NumericVector m_at_age,
    bool tune_unfished, //0 = use a spawner recruit relationship, 1 = don't
    const String rec_form,
    const NumericVector spawning_seasons,
    const NumericVector rec_devs) {

  const int ages = length_at_age.length();

  // --- outputs (these remain R matrices for easy downstream compatibility) ---
  NumericMatrix b_p_a(patches, ages);   // biomass in patch p age a
  NumericMatrix ssb_p_a(patches, ages); // spawning stock biomass in patch p age a

  // catch at age: compute in Eigen first, wrap once
  Eigen::MatrixXd c_p_a_eig = Eigen::MatrixXd::Zero(patches, ages);

  // numbers at age: compute in Eigen (fast), wrap once
  Eigen::MatrixXd n_p_a_eig = Eigen::MatrixXd::Zero(patches, ages);

  // seasonal movement matrix (dense, for step 1)
  // Eigen::MatrixXd movement = Eigen::MatrixXd::Zero(patches, patches);

  Eigen::SparseMatrix<double> movement(patches, patches);

  NumericVector recruits(patches);      // recruits by patch
  NumericVector zero_recruits(patches); // helper

  VectorXd tmp_rec(as<VectorXd>(recruits));

  //////////////////// tune things ////////////////////////
  NumericMatrix tmp_n_p_a = clone(last_n_p_a);
  Rcpp::List tmppop; // annoying step related to memory pointers

  // find the correct movement matrix for the season that you're in
  for (int s = 0; s < movement_seasons.length(); s++) {
    NumericVector tmp_seas = movement_seasons[s];
    bool b = Rcpp::any(tmp_seas == season).is_true();

    if (b) {
      // explicit conversion keeps intent clear
      movement = Rcpp::as<Eigen::SparseMatrix<double>>(movement_matrix[s]);
      // movement = Rcpp::as<Eigen::MatrixXd>(movement_matrix[s]);
      s = movement_seasons.length() + 1; // stop loop once found
    }
  }

  if (Rcpp::traits::is_na<REALSXP>(ssb0) == 1) { // if is.na(ssb0), tune unfished conditions recursively

    int stopper = 0;
    int b = 0;
    double season = 1;
    int seasons = round(1 / time_step);

    // keep burn loop going until it ends in a spawning season
    while (stopper == 0) {

      for (int season = 1; season <= seasons; season++) {

        tmppop = sim_fish(
          length_at_age,
          weight_at_age,
          fec_at_age,
          maturity_at_age,
          semelparous,
          f_p_a,
          movement_matrix,
          movement_seasons,
          recruit_movement_matrix,
          tmp_n_p_a,
          patches,
          0,
          time_step,
          season,
          steepness,
          r0s,
          -999,
          ssb0_p,
          m_at_age,
          1, // this HAS to be 1 inside burn loop
          rec_form,
          spawning_seasons,
          rec_devs);

        b++; // increment b

        tmp_n_p_a = Rcpp::wrap(tmppop["n_p_a"]);

        if (Rcpp::any(season == spawning_seasons).is_true()) {

          NumericMatrix tmp_ssb_p_a = tmppop["ssb_p_a"];

          // collect unfished spawning stock biomass
          ssb0 = sum(tmp_ssb_p_a);
          ssb0_p = rowSums(tmp_ssb_p_a);

        } // keep track of SSB0 in the last spawning season

      } // close season loop

      if (b >= burn_steps) {
        stopper = 1;
      }

    } // close while loop
  } // close SSB0 tuner


  //////////////////// age and die ////////////////////////
  // Explicit loops avoid Rcpp sugar temporaries; still very readable.
  for (int p = 0; p < patches; ++p) {

    // plus group survival and catch
    const int a_last = ages - 1;
    {
      const double Zp = m_at_age[a_last] + f_p_a(p, a_last);
      const double surv_p = std::exp(-time_step * Zp);

      const double N_last = last_n_p_a(p, a_last);
      const double plus_survivors = N_last * surv_p;

      // Baranov catch in biomass units (matches original expression)
      c_p_a_eig(p, a_last) =
        (f_p_a(p, a_last) / Zp) * N_last * (1.0 - surv_p) * weight_at_age[a_last];

      // start plus group in n_p_a_eig last age
      n_p_a_eig(p, a_last) += plus_survivors;
    }

    // ages 1..(A-1) get survivors from age-1
    for (int a = 1; a < ages; ++a) {
      const double Z = m_at_age[a - 1] + f_p_a(p, a - 1);
      const double surv = std::exp(-time_step * Z);
      const double N0 = last_n_p_a(p, a - 1);

      n_p_a_eig(p, a) += N0 * surv;

      // Baranov catch in biomass units (matches original expression)
      c_p_a_eig(p, a - 1) =
        weight_at_age[a] * (f_p_a(p, a - 1) / Z) * N0 * (1.0 - surv);
    }
  }


  //////////////////// move ////////////////////////
  // Big win: keep in Eigen, do one multiply, no wrap/clone ping-pong.
  n_p_a_eig = movement * n_p_a_eig;

  // Wrap once for downstream Rcpp code (recruitment / semelparous blocks)
  NumericMatrix n_p_a = Rcpp::wrap(n_p_a_eig);
  NumericMatrix c_p_a = Rcpp::wrap(c_p_a_eig);


  //////////////////// biomass and ssb ////////////////////////
  for (int p = 0; p < patches; p++) {
    for (int a = 0; a < ages; a++) {
      b_p_a(p, a) = n_p_a(p, a) * weight_at_age[a];
      ssb_p_a(p, a) = n_p_a(p, a) * fec_at_age[a] * maturity_at_age[a];
    }
  }


  //////////////////// spawn / recruit ////////////////////////
  double ssb = sum(ssb_p_a);               // total spawning stock biomass system wide
  NumericVector ssb_p = rowSums(ssb_p_a);  // spawning stock biomass by patch

  if (tune_unfished == 1) { // turn off stock recruitment relationship

    if (Rcpp::any(season == spawning_seasons).is_true()) {
      recruits = r0s;
    } else {
      recruits = zero_recruits;
    }

  } else { // if stock recruitment relationship is in effect

    if (Rcpp::any(season == spawning_seasons).is_true()) {

      if (rec_form == "global_habitat") { // global Beverton-Holt density dependence; distribute recruits by habitat

        recruits = (((0.8 * sum(r0s) * steepness * ssb) /
          (0.2 * ssb0 * (1 - steepness) + (steepness - 0.2) * ssb))) * (r0s / sum(r0s));

      } else if (rec_form == "local_habitat") { // local Beverton-Holt; r0 set by local habitat

        recruits = ((0.8 * r0s * steepness * ssb_p) /
          (0.2 * (ssb0_p + 1e-6) * (1 - steepness) + (steepness - 0.2) * (ssb_p + 1e-6)));

      } else if (rec_form == "pre_dispersal") { // local Beverton-Holt then disperse recruits

        recruits = ((0.8 * r0s * steepness * ssb_p) /
          (0.2 * (ssb0_p + 1e-6) * (1 - steepness) + (steepness - 0.2) * (ssb_p + 1e-6)));

        tmp_rec = as<VectorXd>(clone(recruits));
        tmp_rec = recruit_movement_matrix * tmp_rec;
        recruits = Rcpp::wrap(tmp_rec);

      } else if (rec_form == "post_dispersal") { // disperse larvae then recruit locally per Beverton-Holt

        tmp_rec = as<VectorXd>(clone(ssb_p));
        tmp_rec = recruit_movement_matrix * tmp_rec;

        NumericVector huh(patches);
        huh = Rcpp::wrap(tmp_rec);

        recruits = ((0.8 * r0s * steepness * huh) /
          (0.2 * (ssb0_p + 1e-6) * (1 - steepness) + (steepness - 0.2) * (huh + 1e-6)));

      } else if (rec_form == "global_ssb") {

        // global density dependence, redistribute recruits by spawning biomass
        const double denom =
          0.2 * ssb0 * (1.0 - steepness) +
          (steepness - 0.2) * ssb;

        recruits =
          ((0.8 * sum(r0s) * steepness * ssb) / denom) *
          (ssb_p / sum(ssb_p));

      }


    } else {
      recruits = zero_recruits;
    }
  }

  // assign recruits
  n_p_a(_, 0) = recruits * rec_devs;

  // update biomass/ssb for age-0 (matches original)
  b_p_a(_, 0) = n_p_a(_, 0) * weight_at_age(0);
  ssb_p_a(_, 0) = b_p_a(_, 0) * fec_at_age(0);


  // add in senescence for semelparous species
  if (semelparous == 1) {

    for (int p = 0; p < patches; p++) {

      n_p_a(p, _) = n_p_a(p, _) * (1 - maturity_at_age);

      b_p_a(p, _) = n_p_a(p, _) * weight_at_age;

      ssb_p_a(p, _) = n_p_a(p, _) * fec_at_age * maturity_at_age;

    }

  }


  //////////////////// process results ////////////////////////
  return Rcpp::List::create(
    Rcpp::Named("season") = season,
    Rcpp::Named("n_p_a") = n_p_a,
    Rcpp::Named("b_p_a") = b_p_a,
    Rcpp::Named("ssb_p_a") = ssb_p_a,
    Rcpp::Named("ssb0") = ssb0,
    Rcpp::Named("ssb0_p") = ssb0_p,
    Rcpp::Named("tmppop") = tmppop,
    Rcpp::Named("c_p_a") = c_p_a);
}

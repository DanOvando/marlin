#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
List sim_fish(
    const NumericVector length_at_age, 
    const NumericVector weight_at_age,
    const NumericVector maturity_at_age,
    const NumericMatrix f_p_a, // fishing mortality by patch and age
    const Eigen::MatrixXd movement, // movement matrix
    Rcpp::NumericMatrix last_n_p_a, // last numbers by patch and age
    const int patches,
    const int burn_steps, // number of burn steps if burn is in effect
    const double time_step, // the time step in portions of a year
    const double steepness,
    const double r0,
    double ssb0, // unfished spawning stock biomass across entire system
    NumericVector ssb0_p, // unfished spawning stock biomass in each patch
    const NumericVector m_at_age,
    bool tune_unfished, //0 = use a spawner recruit relationship, 1 = don't
    const int rec_form) // recruitment form, one of ....
  { 

  int ages = length_at_age.length();

  MatrixXd  tmpmat(as<MatrixXd>(last_n_p_a)); // create temporary matrix for movement

  NumericMatrix n_p_a(patches, ages); // numbers in patch p age a

  NumericMatrix b_p_a(patches, ages); // biomass in patch p age a

  NumericMatrix ssb_p_a(patches, ages); // spawning stock biomass in patch p age a

  NumericMatrix c_p_a(patches, ages); // catch at age 

  //////////////////// tune things ////////////////////////
  NumericMatrix tmp_n_p_a = clone(last_n_p_a);

  Rcpp::List tmppop; // annoying step related to memory pointeds

  if (Rcpp::traits::is_na<REALSXP>(ssb0) == 1){ // basically, if is.na(ssb0), tune unfished conditions


    for (int b = 0; b < burn_steps; b++){

      // these HAVE to be in the same order
      // as function call: get lots of warnings
      // if you try and name them

      tmppop = sim_fish(
        length_at_age,
        weight_at_age,
        maturity_at_age,
        f_p_a,
        movement,
        tmp_n_p_a,
        patches,
        0,
        time_step,
        steepness,
        r0,
        -999,
        ssb0_p,
        m_at_age,
        1, // this HAS to be 1
        rec_form);

      tmp_n_p_a = wrap(tmppop["n_p_a"]);

    }

    NumericMatrix tmp_ssb_p_a = tmppop["ssb_p_a"];

    // collect unfished spawning stock biomass 
    ssb0 = sum(tmp_ssb_p_a); 
    
    ssb0_p = rowSums(tmp_ssb_p_a);

  }

  //////////////////// move ////////////////////////


  tmpmat =  tmpmat.transpose() * movement; // matrix multiplication of numbers at age by movement matrix

  SEXP tmp = Rcpp::wrap(tmpmat.transpose()); // convert from eigen to Rcpp

  last_n_p_a = tmp; // set last population to post-movement last population

  //////////////////// grow ////////////////////////


  NumericVector plus_group = last_n_p_a(_,ages - 1) * exp(-(m_at_age(ages - 1) +f_p_a(_,ages - 1))); // calculate numbers in the oldest group that survive
  
  c_p_a(_,ages - 1) =   (f_p_a(_,ages - 1) / (m_at_age(ages - 1) * time_step + f_p_a(_,ages - 1))) * last_n_p_a(_,ages - 1) * (1 - exp(-(m_at_age(ages - 1) * time_step + f_p_a(_,ages - 1))));
  
  // age and die
  for (int a = 1; a < ages; a++){

    n_p_a(_,a) =  last_n_p_a(_,a - 1) * exp(-(m_at_age(a - 1) + f_p_a(_,a - 1)));
    
    c_p_a(_,a - 1) =   (f_p_a(_,a - 1) / (m_at_age(a - 1) * time_step + f_p_a(_,a - 1))) * last_n_p_a(_,a - 1) * (1 - exp(-(m_at_age(a - 1) * time_step + f_p_a(_,a - 1))));

  }

  // add numbers that survived in the oldest group to numbers that grew into the oldest group
  n_p_a(_, ages - 1) = n_p_a(_, ages - 1) + plus_group;

  // calculate biomass and spawning stock biomass at age
  for (int p = 0;p < patches; p++){
    
    c_p_a(p,_) = c_p_a(p,_) * weight_at_age;

    b_p_a(p,_) =  n_p_a(p,_) * weight_at_age;

    ssb_p_a(p,_) =  b_p_a(p,_) * maturity_at_age;

  }

  //////////////////// spawn / recruit ////////////////////////

  double ssb = sum(ssb_p_a); // total spawning stock biomass system wide
  
  NumericVector ssb_p = rowSums(ssb_p_a); // spawning stock biomass by patch
  
  NumericVector recruits(patches); // blank vector for recruits by patch
  
    if (tune_unfished == 1){ // turn off stock recruitment relationship

      if (rec_form == 0){
      
      recruits = rep(r0 / patches, patches);
      
      } else if (rec_form == 1){
        
        // recruits = r0 * ssb_p / sum(ssb_p);
        
        recruits = rep(r0 / patches, patches);
        
      }
      
    } else { // if stock recruitment relationship is in effect


      if (rec_form == 0){  // global beverton-holt density dependence, distribute recruits evenly

      recruits = rep(((0.8 * r0 * steepness * ssb) / (0.2 * ssb0 * (1 - steepness) + (steepness - 0.2) * ssb)) / patches, patches);
        
      } else if (rec_form == 1){ // local beverton-holt density dependence, same r0 everywhere
        
        recruits = ((0.8 * (r0 / patches) * steepness * ssb_p) / (0.2 * ssb0_p * (1 - steepness) + (steepness - 0.2) * ssb_p));
        
      } // close recruitment ifs
    }
    
    // assign recruits
    n_p_a(_,0) = recruits;
    
    b_p_a(_,0) =   n_p_a(_,0) * weight_at_age(0);
    
    ssb_p_a(_,0) =   b_p_a(_,0) * maturity_at_age(0);
    

  //////////////////// process results ////////////////////////

  return Rcpp::List::create(
    Rcpp::Named("n_p_a") = n_p_a,
    Rcpp::Named("b_p_a") = b_p_a,
    Rcpp::Named("ssb_p_a") = ssb_p_a,
    Rcpp::Named("ssb0") = ssb0,
    Rcpp::Named("ssb0_p") = ssb0_p,
    Rcpp::Named("tmppop") = tmppop,
    Rcpp::Named("c_p_a") = c_p_a);
} // close fish model


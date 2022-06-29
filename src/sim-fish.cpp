#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
NumericVector cpp_seq(int steps,double step_size) {
  
  NumericVector y(steps);
  
  y(0) = 1;
  
  for (int i = 1; i < y.length(); i++){
    
    y(i) = y(i - 1) + step_size;
    
  }
  
  return y ;
}


// [[Rcpp::export]]
List sim_fish(
    const NumericVector length_at_age, 
    const NumericVector weight_at_age,
    const NumericVector fec_at_age,
    const NumericVector maturity_at_age,
    const NumericMatrix f_p_a, // fishing mortality by patch and age
    const List seasonal_movement,
    const List movement_seasons,
    Eigen::MatrixXd recruit_movement,
    Rcpp::NumericMatrix last_n_p_a, // last numbers by patch and age
    const int patches,
    const int burn_steps, // number of burn steps if burn is in effect
    const double time_step, // the time step in portions of a year
    double season, // what season you're currently in
    const double steepness,
    const NumericVector r0s,
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

  Eigen::MatrixXd movement = Eigen::MatrixXd::Zero(patches, patches);
  
  NumericVector recruits(patches); // blank vector for recruits by patch
  
  VectorXd tmp_rec(as<VectorXd>(recruits)); 
  
  
  //////////////////// tune things ////////////////////////
  NumericMatrix tmp_n_p_a = clone(last_n_p_a);

  Rcpp::List tmppop; // annoying step related to memory pointers

  
  // find the correct movement matrix for the season that you're in
  
  for (int s  = 0; s < movement_seasons.length(); s++){ // find which season block you're in
    
    NumericVector tmp_seas = movement_seasons[s];
    
    bool b = Rcpp::any( tmp_seas == season).is_true();
    
    if (b == 1){
      
      movement = seasonal_movement[s];
      
      s = movement_seasons.length() + 1; // stop loop once you've found what season you're in
    } 
    
  } // close seasonal movement finder loop
  
  if (Rcpp::traits::is_na<REALSXP>(ssb0) == 1){ // basically, if is.na(ssb0), tune unfished conditions

    NumericVector burn_seq = cpp_seq(burn_steps,time_step);
    
    for (int b = 0; b < burn_steps; b++){

      season = burn_seq(b) - floor(burn_seq(b)); // determine what season you're in
      
       
      // these HAVE to be in the same order
      // as function call: get lots of warnings
      // if you try and name them
      
      tmppop = sim_fish(
        length_at_age,
        weight_at_age,
        fec_at_age,
        maturity_at_age,
        f_p_a,
        seasonal_movement,
        movement_seasons,
        recruit_movement,
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
        rec_form);

      tmp_n_p_a = Rcpp::wrap(tmppop["n_p_a"]);

    }

    NumericMatrix tmp_ssb_p_a = tmppop["ssb_p_a"];

    // collect unfished spawning stock biomass 
    ssb0 = sum(tmp_ssb_p_a); 
    
    ssb0_p = rowSums(tmp_ssb_p_a);

  }

  //////////////////// move ////////////////////////

  tmpmat =  movement * tmpmat; // matrix multiplication of numbers at age by movement matrix
  
  SEXP tmp = Rcpp::wrap(tmpmat); // convert from eigen to Rcpp

  Rcpp::NumericMatrix tmp2(tmp); // attempt to resolve weird issue with random erros based on this https://stackoverflow.com/questions/62586950/how-do-you-convert-object-of-class-eigenmatrixxd-to-class-rcppnumericmatrix
  
  last_n_p_a = clone(tmp2); // set last population to post-movement last population

  //////////////////// grow ////////////////////////


  NumericVector plus_group = last_n_p_a(_,ages - 1) * exp(-time_step * (m_at_age(ages - 1) +f_p_a(_,ages - 1))); // calculate numbers in the oldest group that survive
  
  c_p_a(_,ages - 1) =   (f_p_a(_,ages - 1) / (m_at_age(ages - 1) + f_p_a(_,ages - 1))) * last_n_p_a(_,ages - 1) * (1 - exp(-time_step * (m_at_age(ages - 1) + f_p_a(_,ages - 1))));
  
  // age and die
  for (int a = 1; a < ages; a++){

    n_p_a(_,a) =  last_n_p_a(_,a - 1) * exp(-time_step * (m_at_age(a - 1) + f_p_a(_,a - 1)));
    
    c_p_a(_,a - 1) =   (f_p_a(_,a - 1) / (m_at_age(a - 1) + f_p_a(_,a - 1))) * last_n_p_a(_,a - 1) * (1 - exp(-time_step * (m_at_age(a - 1) + f_p_a(_,a - 1))));

  }

  // add numbers that survived in the oldest group to numbers that grew into the oldest group
  n_p_a(_, ages - 1) = n_p_a(_, ages - 1) + plus_group;

  // calculate biomass and spawning stock biomass at age
  for (int p = 0;p < patches; p++){
    
    c_p_a(p,_) = c_p_a(p,_) * weight_at_age;

    b_p_a(p,_) =  n_p_a(p,_) * weight_at_age;

    // ssb_p_a(p,_) =  b_p_a(p,_) * maturity_at_age;

    ssb_p_a(p,_) = n_p_a(p,_) * fec_at_age * maturity_at_age;
    
  }

  //////////////////// spawn / recruit ////////////////////////

  double ssb = sum(ssb_p_a); // total spawning stock biomass system wide
  
  NumericVector ssb_p = rowSums(ssb_p_a); // spawning stock biomass by patch
  
    if (tune_unfished == 1){ // turn off stock recruitment relationship

      recruits = r0s;
      
    } else { // if stock recruitment relationship is in effect


      if (rec_form == 0){  // global beverton-holt density dependence, distribute recruits according to recruitment habitat

      recruits = (((0.8 * sum(r0s) * steepness * ssb) / (0.2 * ssb0 * (1 - steepness) + (steepness - 0.2) * ssb))) * (r0s / sum(r0s));
        
      } else if (rec_form == 1){ // local beverton-holt density dependence, r0 set by local habitat
        
        recruits = ((0.8 * r0s * steepness * ssb_p) / (0.2 * (ssb0_p + 1e-6) * (1 - steepness) + (steepness - 0.2) * (ssb_p + 1e-6)));
        
      } else if (rec_form == 2) { // local beverton-holt then disperse recruits
        
        recruits = ((0.8 * r0s * steepness * ssb_p) / (0.2 * (ssb0_p + 1e-6) * (1 - steepness) + (steepness - 0.2) * (ssb_p + 1e-6)));
        
        tmp_rec = as<VectorXd>(clone(recruits));
        
        tmp_rec = recruit_movement * tmp_rec;

        recruits = Rcpp::wrap(tmp_rec); // convert from eigen to Rcpp

        // Rcpp::Rcout << "hello" << std::endl;

      } else if (rec_form == 3){ // disperse larvae then recruit locally per beverton-holt
        
        tmp_rec = as<VectorXd>(clone(ssb_p));
        
        tmp_rec = recruit_movement * tmp_rec;
        
        NumericVector huh(patches); // no idea why I have to do this
        
        huh = Rcpp::wrap(tmp_rec); // convert from eigen to Rcpp
        
        recruits = ((0.8 * r0s * steepness * huh) / (0.2 * (ssb0_p + 1e-6) * (1 - steepness) + (steepness - 0.2) * (huh + 1e-6)));
        
      }
        
    }   // close recruitment ifs
    
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


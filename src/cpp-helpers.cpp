#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Eigen::MatrixXd move_test(List seasonal_movement,List movement_seasons, double season) {
  
  Eigen::MatrixXd movement = Eigen::MatrixXd::Zero(3, 3);
  

  for (int s  = 0; s < movement_seasons.length(); s++){ // find which season block you're in
    
    NumericVector tmp_seas = movement_seasons[s];
  
    bool b = Rcpp::any( tmp_seas == season).is_true();
    
    Rcpp::Rcout << tmp_seas << std::endl;
    
    Rcpp::Rcout << season << std::endl;
    
    Rcpp::Rcout << "b is" << b << std::endl;
    
    if (b == 1){
      
      Rcpp::Rcout << "b is" << b << std::endl;
      
      movement = seasonal_movement[s];
      
      s = movement_seasons.length() + 1; // stop loop once you've found what season you're in
    } 
    
  }
  // return Rcpp::List::create(
  //   Rcpp::Named("movement") = movement,
  //   Rcpp::Named("b") = b);
  //   
  return(movement);

}

/*** R

seasonal_movement = list(matrix(1, nrow = 4, ncol = 2), matrix(2, nrow = 4, ncol = 2))

movement_seasons = list(c(0.,.25), c(0.5, .75))


# seasonal_movement = list(matrix(1, nrow = 4, ncol = 2))
# 
# movement_seasons = list(c(0.,.25))

move_test(seasonal_movement,movement_seasons, .75)



*/
// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// sample_problem
Eigen::MatrixXd sample_problem(Eigen::MatrixXd x, Eigen::MatrixXd y);
RcppExport SEXP _marlin_sample_problem(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(sample_problem(x, y));
    return rcpp_result_gen;
END_RCPP
}
// move_test
Eigen::MatrixXd move_test(List seasonal_movement, List movement_seasons, double season);
RcppExport SEXP _marlin_move_test(SEXP seasonal_movementSEXP, SEXP movement_seasonsSEXP, SEXP seasonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type seasonal_movement(seasonal_movementSEXP);
    Rcpp::traits::input_parameter< List >::type movement_seasons(movement_seasonsSEXP);
    Rcpp::traits::input_parameter< double >::type season(seasonSEXP);
    rcpp_result_gen = Rcpp::wrap(move_test(seasonal_movement, movement_seasons, season));
    return rcpp_result_gen;
END_RCPP
}
// cpp_seq
NumericVector cpp_seq(int steps, double step_size);
RcppExport SEXP _marlin_cpp_seq(SEXP stepsSEXP, SEXP step_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< double >::type step_size(step_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_seq(steps, step_size));
    return rcpp_result_gen;
END_RCPP
}
// sim_fish
List sim_fish(const NumericVector length_at_age, const NumericVector weight_at_age, const NumericVector maturity_at_age, const NumericMatrix f_p_a, const List seasonal_movement, const List movement_seasons, Rcpp::NumericMatrix last_n_p_a, const int patches, const int burn_steps, const double time_step, double season, const double steepness, const NumericVector r0s, double ssb0, NumericVector ssb0_p, const NumericVector m_at_age, bool tune_unfished, const int rec_form);
RcppExport SEXP _marlin_sim_fish(SEXP length_at_ageSEXP, SEXP weight_at_ageSEXP, SEXP maturity_at_ageSEXP, SEXP f_p_aSEXP, SEXP seasonal_movementSEXP, SEXP movement_seasonsSEXP, SEXP last_n_p_aSEXP, SEXP patchesSEXP, SEXP burn_stepsSEXP, SEXP time_stepSEXP, SEXP seasonSEXP, SEXP steepnessSEXP, SEXP r0sSEXP, SEXP ssb0SEXP, SEXP ssb0_pSEXP, SEXP m_at_ageSEXP, SEXP tune_unfishedSEXP, SEXP rec_formSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type length_at_age(length_at_ageSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type weight_at_age(weight_at_ageSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type maturity_at_age(maturity_at_ageSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type f_p_a(f_p_aSEXP);
    Rcpp::traits::input_parameter< const List >::type seasonal_movement(seasonal_movementSEXP);
    Rcpp::traits::input_parameter< const List >::type movement_seasons(movement_seasonsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type last_n_p_a(last_n_p_aSEXP);
    Rcpp::traits::input_parameter< const int >::type patches(patchesSEXP);
    Rcpp::traits::input_parameter< const int >::type burn_steps(burn_stepsSEXP);
    Rcpp::traits::input_parameter< const double >::type time_step(time_stepSEXP);
    Rcpp::traits::input_parameter< double >::type season(seasonSEXP);
    Rcpp::traits::input_parameter< const double >::type steepness(steepnessSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type r0s(r0sSEXP);
    Rcpp::traits::input_parameter< double >::type ssb0(ssb0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ssb0_p(ssb0_pSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type m_at_age(m_at_ageSEXP);
    Rcpp::traits::input_parameter< bool >::type tune_unfished(tune_unfishedSEXP);
    Rcpp::traits::input_parameter< const int >::type rec_form(rec_formSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_fish(length_at_age, weight_at_age, maturity_at_age, f_p_a, seasonal_movement, movement_seasons, last_n_p_a, patches, burn_steps, time_step, season, steepness, r0s, ssb0, ssb0_p, m_at_age, tune_unfished, rec_form));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_marlin_sample_problem", (DL_FUNC) &_marlin_sample_problem, 2},
    {"_marlin_move_test", (DL_FUNC) &_marlin_move_test, 3},
    {"_marlin_cpp_seq", (DL_FUNC) &_marlin_cpp_seq, 2},
    {"_marlin_sim_fish", (DL_FUNC) &_marlin_sim_fish, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_marlin(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

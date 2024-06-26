# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

sample_problem <- function(x, y) {
    .Call(`_marlin_sample_problem`, x, y)
}

move_test <- function(seasonal_movement, movement_seasons, season) {
    .Call(`_marlin_move_test`, seasonal_movement, movement_seasons, season)
}

move <- function(x, m) {
    .Call(`_marlin_move`, x, m)
}

cpp_seq <- function(steps, step_size) {
    .Call(`_marlin_cpp_seq`, steps, step_size)
}

sim_fish <- function(length_at_age, weight_at_age, fec_at_age, maturity_at_age, f_p_a, movement_matrix, movement_seasons, recruit_movement_matrix, last_n_p_a, patches, burn_steps, time_step, season, steepness, r0s, ssb0, ssb0_p, m_at_age, tune_unfished, rec_form, spawning_seasons, rec_devs) {
    .Call(`_marlin_sim_fish`, length_at_age, weight_at_age, fec_at_age, maturity_at_age, f_p_a, movement_matrix, movement_seasons, recruit_movement_matrix, last_n_p_a, patches, burn_steps, time_step, season, steepness, r0s, ssb0, ssb0_p, m_at_age, tune_unfished, rec_form, spawning_seasons, rec_devs)
}


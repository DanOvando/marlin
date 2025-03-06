#' Title
#'
#' @param fmult multiplier on fishing mortality
#' @param quota target quota
#' @param fauna fauna object
#' @param current_season current season
#' @param movement movement matrix
#' @param f_p_a fishing mortality rate at patch and age
#' @param last_n_p_a last numbers at patch and age
#' @param f_p_a_fl fishing mortality rate by patch age and fleet
#'
#' @return sum of squares of catch relative to quota

quota_finder <- function(fmult, quota, fauna, current_season, movement, f_p_a, last_n_p_a, f_p_a_fl, f, rec_devs, patches, ages, fleets) {
  tmp_pop <-
    fauna[[f]]$swim(
      season = current_season,
      adult_movement = movement,
      f_p_a = fmult * f_p_a,
      last_n_p_a = last_n_p_a,
      rec_devs = rec_devs
    )

  tmp_catch <-
    f_p_a_fl * array(
      tmp_pop$c_p_a,
      dim = c(patches, ages, length(fleets)),
      dimnames = list(1:patches, fauna[[f]]$ages, names(fleets))
    )

  log_ss <- (sum(tmp_catch, na.rm = TRUE) - quota)^2
}

#' Objective Function for Quota-Constrained Fishing Mortality
#'
#' @description
#' Computes the squared difference between realised total catch and a target
#' quota under a uniformly scaled fishing mortality. Used as the objective
#' for \code{optim} inside \code{\link{simmar}} whenever a species-level
#' catch quota is active and would be exceeded.
#'
#' @details
#' When a quota is binding, \code{\link{simmar}} calls \code{optim} with this
#' function to find the scalar \code{fmult} in [0, 1] that reduces the
#' fishing mortality matrix \code{f_p_a} (and hence catch) to exactly the
#' quota level. All fleets experience the same proportional reduction.
#'
#' @param fmult Numeric in [0, 1]. Effort multiplier applied uniformly to
#'   \code{f_p_a} before simulating the population.
#' @param quota Numeric. Target total catch (in numbers) that must not be
#'   exceeded.
#' @param fauna Named list of fauna objects.
#' @param current_season Integer. Current season index.
#' @param movement List of movement matrices (one per season block).
#' @param f_p_a Numeric matrix \code{[patches, ages]} of baseline fishing
#'   mortality before quota reduction.
#' @param last_n_p_a Numeric matrix \code{[patches, ages]} of numbers at the
#'   start of the current step.
#' @param f_p_a_fl 3-D array \code{[patches, ages, fleets]} of proportional
#'   fishing mortality share by fleet.
#' @param critter Character. Name of the target species.
#' @param patches Integer. Number of patches.
#' @param ages Integer. Number of age classes.
#' @param fleets Named list of fleet objects.
#' @param rec_devs Numeric vector of recruitment deviates (length = patches).
#'
#' @return Numeric scalar: squared deviation of total catch from \code{quota}.
#'   Minimised by \code{optim} to find the binding effort multiplier.
#'
#' @keywords internal
#' @seealso \code{\link{simmar}}
quota_finder <- function(fmult, quota, fauna, current_season, movement, f_p_a, last_n_p_a, f_p_a_fl, critter, rec_devs, patches, ages, fleets) {
  tmp_pop <-
    fauna[[critter]]$swim(
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
      dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
    )

  log_ss <- (sum(tmp_catch, na.rm = TRUE) - quota)^2
}

#' Objective function for depletion-based fleet tuning
#'
#' Given a vector of log fishing mortalities (one per species), sets
#' fleet catchabilities accordingly (via the logistic link function to
#' keep q in (0, 1)), runs a forward simulation, and returns log-ratio
#' residuals between achieved and target depletion.
#'
#' Used as the objective function for \code{\link[nleqslv]{nleqslv}}
#' inside \code{\link{tune_fleets}}.
#'
#' @param log_fs numeric vector of log instantaneous fishing mortality,
#'   one element per species in \code{fauna}
#' @param fauna a fauna object
#' @param fleets a fleet object (already set to constant effort)
#' @param e_fl numeric vector of effective effort per fleet
#' @param years number of years to simulate
#'
#' @return numeric vector of residuals (log achieved depletion - log target
#'   depletion), one per species
#' @keywords internal
#' @export
fleet_tuner <- function(log_fs, fauna, fleets, e_fl, years = 50) {

  fs <- exp(log_fs)

  tfleets <- clone_fleet_metiers(fleets)

  # Set catchabilities for every fleet x species combination
  # using the link function to keep q in (0, 1)
  for (f in seq_along(tfleets)) {
    for (ff in seq_along(fauna)) {

      f_critter <- fs[ff]
      f_metier <- tfleets[[f]]$metiers[[ff]]$p_explt * f_critter
      raw_q <- f_metier / e_fl[f]

      if (!is.finite(raw_q)) raw_q <- 0

      set_metier_catchability(tfleets[[f]]$metiers[[ff]], raw_q, use_link = TRUE)
    }
  }

  storage <- simmar(
    fauna = fauna,
    fleets = tfleets,
    years = years
  )

  tmp <- purrr::map(
    storage[[length(storage)]],
    ~ as.data.frame(.x$ssb_p_a)
  ) |>
    purrr::list_rbind(names_to = "fauna")

  b_p <- rowSums(tmp[, 2:ncol(tmp)], na.rm = TRUE)

  tmp <- data.frame(fauna = tmp$fauna, ssb = b_p) |>
    dplyr::group_by(fauna) |>
    dplyr::summarise(ssb = sum(ssb)) |>
    dplyr::ungroup() |>
    dplyr::arrange(fauna)

  ssb0s <- purrr::map_dbl(fauna, "ssb0")
  ssb0s <- ssb0s[sort(names(ssb0s))]

  target_depletion <- purrr::map_dbl(fauna, "depletion")
  target_depletion <- target_depletion[sort(names(target_depletion))]

  tmp$depletion <- tmp$ssb / ssb0s

  r <- log(pmax(tmp$depletion, 1e-12)) - log(target_depletion)

  # If something pathological happened, return large residuals instead of stopping
  if (any(!is.finite(r))) {
    r <- rep(1e6, length(r))
  }

  r
}

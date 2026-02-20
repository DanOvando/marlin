#' Assign MSY-Based Reference Points to Each Species
#'
#' @description
#' Finds the effort level that maximises yield for each species in \code{fauna}
#' by calling \code{optim} on \code{\link{find_msy}}, then stores the
#' resulting MSY, SSBmsy, Bmsy, and Umsy in each critter's \code{ref_points}
#' field.
#'
#' @details
#' The optimisation searches over a scalar effort multiplier applied to the
#' first fleet's \code{base_effort} (range 0.001--10). MSY is estimated for
#' each species independently (i.e. holding all else equal; no joint
#' multi-species MSY is computed). Results are stored in each critter's
#' \code{ref_points} data frame with columns \code{ssb_msy}, \code{b_msy},
#' \code{n_msy}, \code{msy}, \code{u_msy}, \code{base_e_msy_mult}, and
#' \code{base_e_msy}.
#'
#' @param fauna Named list of fauna objects from \code{\link{create_critter}}.
#' @param fleets Named list of fleet objects from \code{\link{create_fleet}},
#'   already tuned with \code{\link{tune_fleets}}.
#'
#' @return A copy of \code{fauna} with the \code{ref_points} field populated
#'   for each species.
#'
#' @seealso \code{\link{find_msy}}, \code{\link{create_critter}},
#'   \code{\link{tune_fleets}}
#'
#' @export
#' Assign References Points
#'
#' @param fauna a list of critters
#' @param fleets a list of cleets
#'
#' @return a fauna object but with MSY based reference points included for each critter
#' @export
#'
assign_ref_points <- function(fauna, fleets) {
  for (f in seq_along(fauna)) {
    msy_mult <- optim(
      1e-3,
      find_msy,
      lower = 0,
      upper = 10,
      fauna = fauna,
      fleets = fleets,
      opt = TRUE,
      target_critter = names(fauna)[f],
      method = "L-BFGS-B"
    )

    msy_state <- find_msy(msy_mult$par, fauna = fauna, fleets = fleets, opt = FALSE, target_critter = names(fauna)[f])

    ref_points <- msy_state$fauna %>%
      dplyr::filter(critter == names(fauna)[f]) %>%
      dplyr::summarise(
        ssb_msy = sum(ssb),
        b_msy = sum(b),
        n_msy = sum(n),
        msy = sum(c),
        u_msy = sum(c) / sum(b)
      )

    base_e_msy <- mean(purrr::map_dbl(fleets, "base_effort")) * msy_mult$par

    ref_points$base_e_msy_mult <- msy_mult$par

    ref_points$base_e_msy <- base_e_msy

    fauna[[f]]$ref_points <- ref_points
  }

  return(fauna)
}

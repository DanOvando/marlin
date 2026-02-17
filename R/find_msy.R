#' Find MSY Conditions for a Single Species
#'
#' @description
#' Estimates maximum sustainable yield (MSY) conditions for one species by
#' running a simulation under a scaled effort level and returning either the
#' negative total yield (for minimisation by \code{optim}) or the full
#' \code{\link{process_marlin}} output at that effort level.
#'
#' Mostly an internal helper for \code{\link{assign_ref_points}}, but
#' exported for use in custom MSY calculations.
#'
#' @param effort_mult Numeric. Scalar multiplier applied to the first fleet's
#'   \code{base_effort} in \code{fleets}. The optimiser searches over this
#'   value to find the effort that maximises yield.
#' @param fauna Named list of fauna objects from \code{\link{create_critter}}.
#' @param fleets Named list of fleet objects from \code{\link{create_fleet}},
#'   already tuned with \code{\link{tune_fleets}}.
#' @param opt Logical. If \code{TRUE} (default), returns negative total yield
#'   (scalar) for use as an \code{optim} objective. If \code{FALSE}, returns
#'   the full \code{\link{process_marlin}} output at the given effort level.
#' @param target_critter Character. Name of the species whose yield is
#'   maximised. Must match a name in \code{fauna}.
#'
#' @return When \code{opt = TRUE}: a negative numeric scalar (negative total
#'   yield for the target critter, for minimisation). When \code{opt = FALSE}:
#'   a \code{\link{process_marlin}} output list at the given effort level.
#'
#' @seealso \code{\link{assign_ref_points}}, \code{\link{simmar}},
#'   \code{\link{process_marlin}}
#'
#' @export
#' Find MSY
#' estimates MSY for an individual critter
#'
#' Note that this is mostly an internal function to \code{assign_ref_points}
#'
#' @param effort_mult the multiplier on the base effort in fleet (assumes fleet model is already tuned)
#' @param fauna a list of critters
#' @param fleets a list of fleets
#' @param opt TRUE = optimize, FALSE = return MSY conditions
#' @param target_critter the name of the critter to find MSY for
#'
#' @return MSY conditions
#' @export
#'
find_msy <- function(effort_mult, fauna, fleets, opt = TRUE, target_critter) {
  tmp_fleets <- fleets

  # for (f in seq_along(fleets)){
  #
  #   tmp_fleets[[f]]$base_effort <- tmp_fleets[[f]]$base_effort * effort_mult
  #
  # }

  tmp_fleets <- purrr::modify_in(tmp_fleets, list(1, "base_effort"), ~ .x * effort_mult)

  sim <- simmar(
    fauna = fauna,
    fleets = tmp_fleets,
    years = years
  )


  processed_marlin <-
    process_marlin(
      sim = sim,
      time_step = time_step,
      steps_to_keep = last(names(sim)),
      keep_age = FALSE
    )

  fauna_yield <- processed_marlin$fleets %>%
    dplyr::filter(critter == target_critter) %>%
    dplyr::summarise(yield = sum(catch))

  out <- -sum(fauna_yield$yield)

  if (opt == FALSE) {
    out <- processed_marlin
  }
  # yield
  #
  return(out)
}

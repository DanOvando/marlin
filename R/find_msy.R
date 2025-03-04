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

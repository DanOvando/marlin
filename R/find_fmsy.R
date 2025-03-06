#' Find MSY Referemce Points
#'
#' Deprecated and likely not to come back, only really possible with one fleet, one critter, and no space, so in other words not what `marlin` is intended for
#'
#' @param fauna a fauna object
#' @param fleets a fleet object
#'
#' @return nothing at the moment
#' @export
#'
find_msy <- function(log_emult, fauna, fleets, years) {
  if (length(fleets > 1) | length(fauna) > 1) {
    stop("find_msy only works with one fleet and one critter at the moment")
  }

  emult <- exp(log_emult)

  sim <- simmar(
    fauna = hyper_critter,
    fleets = hyper_fleet,
    years = years
  )
}

#' Calculate Fishable Revenue by Patch and Fleet
#'
#' @description
#' Computes the revenue standing available to each fleet in each patch: the
#' intersection of biomass at age, the fleet's vulnerability at age, and price.
#' Used by \code{\link{allocate_effort}} when a fleet's
#' \code{spatial_allocation} is \code{"fishable_revenue"}.
#'
#' @details
#' For fleet \eqn{f} and patch \eqn{p}, summing over species \eqn{s} and ages
#' \eqn{a}:
#' \deqn{FR_{p,f} = \sum_s price_{s,f} \sum_a B_{p,a,s} \, v_{p,a,s,f}}
#' where \eqn{B} is biomass at age and \eqn{v} is the metier's \code{vul_p_a}
#' (spatial catchability times selectivity at age).
#'
#' ## Why this exists
#' Every other spatial-allocation objective is a function of realised fishing:
#' \code{\link{go_fish}} computes what a fleet \emph{got} given the effort it
#' actually deployed, so the signal reaching \code{\link{allocate_effort}} at
#' step \eqn{t} necessarily describes step \eqn{t-1}. When habitat changes on a
#' timescale close to the time step -- seasonal migration, say -- the fleet is
#' then always responding to where fish were rather than where they are, and no
#' value of \code{responsiveness} fixes it, because the lag is structural rather
#' than a matter of step size.
#'
#' Fishable revenue takes no effort argument. That is the entire point: because
#' it does not depend on how effort is distributed, it can be evaluated against
#' the biomass a step is \emph{about} to fish rather than the one it just
#' fished, with no circularity and no fixed point to solve. It is also cheaper
#' than the realised objectives, since it needs no Baranov solve.
#'
#' It is additionally defined where effort is zero, unlike the per-unit-effort
#' objectives (\code{rpue}, \code{cpue}, \code{ppue}), which are \code{NA}
#' there. A fleet allocating on fishable revenue therefore sees unfished patches
#' honestly, rather than having to fish somewhere to learn anything about it.
#'
#' ## What it gives up
#' \itemize{
#'   \item \strong{No crowding.} The objective is standing stock, not return per
#'     unit effort, so within a step there is nothing to stop effort piling into
#'     the best patch. Depletion still pushes back, but only from the next step
#'     onward -- delayed negative feedback, which can oscillate. Consider
#'     \code{memory_halflife > 0} to damp it.
#'   \item \strong{No costs.} Cost depends on effort, so profit-side variants
#'     (\code{ppue}, \code{profit}) have no fishable-revenue analogue.
#'   \item \strong{No fleet interaction within a step.} Two fleets on the same
#'     stock do not see each other's effort; they feel each other only through
#'     next step's depleted biomass.
#'   \item \strong{Not an ideal free distribution.} \code{rpue} and the marginal
#'     rules relax toward equalised returns. This one relaxes toward effort
#'     tracking standing stock, which is a different equilibrium -- "fish where
#'     the fish are" rather than "fish where the next unit of effort pays the
#'     same everywhere". Neither is wrong; they are different fleets.
#' }
#'
#' @param b_p_a Named list of biomass-at-age matrices (patches x ages), one per
#'   species, named to match \code{fauna}. In \code{\link{simmar}} this is the
#'   \code{b_p_a} stored for the previous step, which is the biomass the current
#'   step is about to fish.
#' @param fauna Named list of fauna objects (e.g. from
#'   \code{\link{create_critter}}).
#' @param fleets Named list of fleet objects (e.g. from
#'   \code{\link{create_fleet}}), each with a \code{metiers} list.
#'
#' @return Numeric matrix, patches x fleets, with fleet names as column names.
#'
#' @seealso \code{\link{allocate_effort}}, \code{\link{go_fish}},
#'   \code{\link{calc_marginal_value}}, \code{\link{simmar}}
#'
#' @export
calc_fishable_revenue <- function(b_p_a, fauna, fleets) {
  fleet_names <- names(fleets)
  fauni <- names(fauna)

  patches <- nrow(b_p_a[[1]])

  out <- matrix(
    0,
    nrow = patches,
    ncol = length(fleets),
    dimnames = list(NULL, fleet_names)
  )

  for (l in seq_along(fleets)) {
    for (critter in fauni) {
      metier <- fleets[[l]]$metiers[[critter]]

      # a fleet need not have a metier for every critter
      if (is.null(metier)) {
        next
      }

      out[, l] <- out[, l] +
        rowSums(b_p_a[[critter]] * metier$vul_p_a) * metier$price
    }
  }

  out
}

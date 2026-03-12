#' Precompute Baranov Catch Equation Inputs
#'
#' Assembles the per-species selectivity, mortality, biomass, and price
#' matrices needed by the IFD effort-allocation solver.
#'
#' @param storage List of per-species population state (one element per
#'   species), each containing a \code{b_p_a} matrix (patches x ages).
#' @param fauna List of fauna objects (e.g. from \code{\link{create_critter}}).
#' @param fleets List of fleet objects (e.g. from \code{\link{create_fleet}}).
#' @param target_fleet Integer index of the fleet being optimised.
#' @param E_exo List of effort vectors for all other fleets (one per fleet),
#'   or \code{NULL} if no exogenous effort.
#' @param P Integer number of spatial patches.
#'
#' @return A list with elements \code{alpha_mats}, \code{other_mort_mats},
#'   \code{biomass_mats} (each a list of matrices, one per species) and
#'   \code{price_s} (numeric vector of prices).
#'
#' @keywords internal
precompute_baranov_inputs <- function(storage, fauna, fleets, target_fleet, E_exo, P) {
  n_species <- length(storage)
  n_fleet <- length(fleets)

  alpha_mats <- vector("list", n_species)
  other_mort_mats <- vector("list", n_species)
  biomass_mats <- vector("list", n_species)
  price_s <- numeric(n_species)

  has_exo <- !is.null(E_exo)
  if (has_exo) {
    stopifnot(is.list(E_exo), length(E_exo) == n_fleet)
  }

  for (s in seq_len(n_species)) {
    biom <- storage[[s]]$b_p_a
    n_age <- ncol(biom)

    met_t <- fleets[[target_fleet]]$metiers[[s]]
    price_s[s] <- met_t$price

    alpha_mat <- met_t$spatial_catchability %o% met_t$sel_at_age

    other_fish <- matrix(0, nrow = P, ncol = n_age)

    if (n_fleet > 1 && has_exo) {
      for (fleet_idx in seq_len(n_fleet)) {
        if (fleet_idx == target_fleet) next
        met_f <- fleets[[fleet_idx]]$metiers[[s]]
        other_fish <- other_fish +
          (met_f$spatial_catchability %o% met_f$sel_at_age) * E_exo[[fleet_idx]]
      }
    }

    nat_mort <- fauna[[s]]$m_at_age
    other_mort <- matrix(nat_mort, nrow = P, ncol = n_age, byrow = TRUE) + other_fish

    alpha_mats[[s]] <- alpha_mat
    other_mort_mats[[s]] <- other_mort
    biomass_mats[[s]] <- biom
  }

  list(
    alpha_mats = alpha_mats,
    other_mort_mats = other_mort_mats,
    biomass_mats = biomass_mats,
    price_s = price_s
  )
}

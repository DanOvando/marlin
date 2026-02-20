# ============================================================
# marlin_softmax_patchwise_ragged_FAST.R
#
# End-to-end "fast swaps" applied:
#   1) Analytic marginal (no finite differences):
#        marginal_unconstrained_ragged_analytic_fast()
#   2) Vectorized swap-marginal in O(P):
#        swap_marginal_from_m_fast()
#
# Changes in this version:
#   - precompute_baranov_inputs_softmax() stores outputs by species name (critter)
#   - fixed typo: v_mats[[critter]] (was [[critter]])
#   - marginal_unconstrained_ragged_analytic_fast() iterates by species names
#     and pulls prices by name (price_s[[critter]]) for extra robustness.
#
# Notes:
# - Ragged ages by species => v_mats/other_mats/b_mats are lists of [P x A_s]
# - metier$vul_p_a is assumed precomputed and consistent with storage/fauna
#
# ============================================================

# ----------------------------
# 1) Precompute (name-safe)
# ----------------------------

compute_baranov_inputs <- function(storage,
                                              fauna,
                                              fleets,
                                              target_fleet,
                                              E_exo = NULL,
                                              P) {
  n_species <- length(storage)

  critters <- names(fauna)

  if (any(names(fauna) != names(storage))){
    stop("fauna and storage must have the same names in the same order")
  }

  stopifnot(!is.null(critters), length(critters) == n_species)

  n_fleet <- length(fleets)

  v_mats     <- setNames(vector("list", n_species), critters)
  other_mats <- setNames(vector("list", n_species), critters)
  b_mats     <- setNames(vector("list", n_species), critters)

  price_s <- setNames(numeric(n_species), critters)

  for (critter in critters) {
    # biomass (P x A_s)
    b_pa <- storage[[critter]]$b_p_a
    A_s <- ncol(b_pa)
    stopifnot(nrow(b_pa) == P)

    # target fleet metier for this species
    met_t <- fleets[[target_fleet]]$metiers[[critter]]
    price_s[[critter]] <- met_t$price

    v_pa <- met_t$vul_p_a
    stopifnot(!is.null(v_pa), nrow(v_pa) == P, ncol(v_pa) == A_s)

    # exogenous fishing mortality from other fleets (as "other_fish_pa")
    other_fish_pa <- matrix(0, nrow = P, ncol = A_s)

    if (!is.null(E_exo) && n_fleet > 1) {
      for (fleet_idx in seq_len(n_fleet)) {
        if (fleet_idx == target_fleet) next

        met_f <- fleets[[fleet_idx]]$metiers[[critter]]
        v_pa_f <- met_f$vul_p_a
        stopifnot(!is.null(v_pa_f), nrow(v_pa_f) == P, ncol(v_pa_f) == A_s)

        eff_exo_p <- as.numeric(E_exo[[fleet_idx]])
        stopifnot(length(eff_exo_p) == P)

        # other_fish_pa += (v_pa_f) * eff_exo_p  (row-wise scaling)
        other_fish_pa <- other_fish_pa + v_pa_f * eff_exo_p
      }
    }

    # natural mortality for this species (length A_s)
    m_a <- fauna[[critter]]$m_at_age
    stopifnot(length(m_a) == A_s)
    m_pa <- matrix(m_a, nrow = P, ncol = A_s, byrow = TRUE)

    # store by species name
    v_mats[[critter]]     <- v_pa
    other_mats[[critter]] <- m_pa + other_fish_pa
    b_mats[[critter]]     <- b_pa
  }

  list(v_mats = v_mats, other_mats = other_mats, b_mats = b_mats, price_s = price_s)
}

# ----------------------------
# 2) FAST analytic marginal m[p] = d obj_p / d eff_p
# ----------------------------

calculate_returns <- function(
    eff_p,
    v_mats, other_mats, b_mats, price_s, dt,
    c0, gamma, travel_p,
    open_p,
    eps = 1e-12
) {
  # Returns per patch:
  # - revenue and profit (levels)
  # - marginal_revenue and marginal_profit (derivatives wrt eff_p)
  #
  # Baranov catch fraction:
  #   F = v * eff
  #   O = other
  #   Z = F + O
  #   g(F) = (F/Z) * (1 - exp(-dt*Z))
  #
  # d/dF g(F) = (O/Z^2)*(1-exp(-dt*Z)) + (F/Z)*dt*exp(-dt*Z)
  # d/deff g(F) = v * d/dF g(F)

  P <- length(eff_p)

  eff_use <- eff_p
  eff_use[!open_p] <- 0

  revenue <- rep(0, P)
  m_rev   <- rep(0, P)

  # iterate by species names to keep everything name-aligned
  for (critter in names(v_mats)) {
    v_pa     <- v_mats[[critter]]      # [P x A_s]
    other_pa <- other_mats[[critter]]  # [P x A_s]
    b_pa     <- b_mats[[critter]]      # [P x A_s]
    price    <- price_s[[critter]]     # scalar

    # F and Z
    F_pa <- v_pa * eff_use
    Z_pa <- pmax(F_pa + other_pa, eps)

    exp_term  <- exp(-dt * Z_pa)
    one_minus <- 1 - exp_term

    # --- levels: g(F) and revenue ---
    g_pa <- (F_pa / Z_pa) * one_minus
    revenue <- revenue + price * rowSums(b_pa * g_pa)

    # --- marginals: d/deff g(F) and marginal revenue ---
    term1 <- (other_pa / (Z_pa^2)) * one_minus
    term2 <- (F_pa / Z_pa) * (dt * exp_term)
    dgdF  <- term1 + term2

    dgdEff <- v_pa * dgdF
    m_rev  <- m_rev + price * rowSums(b_pa * dgdEff)
  }

  # cost levels:
  # cost_p = c0*(eff^gamma + travel*eff)
  eff_pos <- pmax(eff_use, 0)
  cost <- c0 * (eff_pos^gamma + travel_p * eff_pos)

  # marginal cost:
  # d/deff cost_p = c0*(gamma*eff^(gamma-1) + travel)
  m_cost <- c0 * (gamma * (eff_pos^(gamma - 1)) + travel_p)

  profit <- revenue - cost
  m_profit <- m_rev - m_cost

  # match your prior convention: NA on closed patches
  revenue[!open_p]   <- NA_real_
  profit[!open_p]    <- NA_real_
  m_rev[!open_p]     <- NA_real_
  m_cost[!open_p]    <- NA_real_
  m_profit[!open_p]  <- NA_real_
  cost[!open_p]      <- NA_real_

  list(
    marginal_revenue = m_rev,
    marginal_profit  = m_profit,
    revenue          = revenue,
    profit           = profit,
    marginal_cost    = m_cost,
    cost             = cost
  )
}


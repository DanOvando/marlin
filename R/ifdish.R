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

precompute_baranov_inputs_softmax <- function(storage,
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

marginal_unconstrained_ragged_analytic_fast <- function(
    eff_p,
    v_mats, other_mats, b_mats, price_s, dt,
    c0, gamma, travel_p,
    open_p,
    eps = 1e-12
) {
  # Analytic derivative of the Baranov catch fraction with respect to effort.
  #
  # For each species and age at patch p:
  #   F = v * eff
  #   O = other
  #   Z = F + O
  #   g(F) = (F/Z) * (1 - exp(-dt*Z))
  #
  # d/dF g(F) = (O/Z^2)*(1-exp(-dt*Z)) + (F/Z)*dt*exp(-dt*Z)
  # d/deff g(F) = v * d/dF g(F)
  #
  # Marginal revenue at patch p:
  #   sum_species price * rowSums( b * d/deff g(F) )
  #
  P <- length(eff_p)
  eff_use <- eff_p
  eff_use[!open_p] <- 0

  m_rev <- rep(0, P)

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

    # d/dF g(F)
    term1 <- (other_pa / (Z_pa^2)) * one_minus
    term2 <- (F_pa / Z_pa) * (dt * exp_term)
    dgdF  <- term1 + term2

    # d/deff g(F) = v * dgdF
    dgdEff <- v_pa * dgdF

    # price-weighted biomass contribution
    m_rev <- m_rev + price * rowSums(b_pa * dgdEff)
  }

  # marginal cost:
  # cost_p = c0*(eff^gamma + travel*eff) => d/deff = c0*(gamma*eff^(gamma-1) + travel)
  eff_pos <- pmax(eff_use, 0)
  m_cost  <- c0 * (gamma * (eff_pos^(gamma - 1)) + travel_p)

  m <- m_rev - m_cost
  m[!open_p] <- NA_real_
  m
}

# ----------------------------
# 3) FAST swap marginal u_swap[p] in O(P)
# ----------------------------

swap_marginal_from_m_fast <- function(eff_p, m_p, open_p, eps = 1e-12) {
  # Budget-neutral "proportional donor" swap marginal:
  # donors are all other open patches; donor funding shares proportional to donor effort.
  #
  # For open patch p:
  #   donor_mean_m(p) = (sum_j eff_j m_j - eff_p m_p) / (sum_j eff_j - eff_p)
  #   u_swap[p] = m_p - donor_mean_m(p)
  #
  u_swap <- rep(NA_real_, length(eff_p))

  idx <- which(open_p & is.finite(m_p))
  if (length(idx) <= 1) return(u_swap)

  eff_i <- eff_p[idx]
  m_i   <- m_p[idx]

  Eff_tot_open <- sum(eff_i)
  S <- sum(eff_i * m_i)

  denom <- pmax(Eff_tot_open - eff_i, eps)
  donor_mean <- (S - eff_i * m_i) / denom

  u_swap[idx] <- m_i - donor_mean
  u_swap
}

# ----------------------------
# 4) Softmax update (unchanged)
# ----------------------------

softmax_stable <- function(x) {
  z <- x - max(x, na.rm = TRUE)
  ez <- exp(z)
  ez / sum(ez, na.rm = TRUE)
}

normalize_util <- function(u, open_p, method = c("iqr", "sd"), eps = 1e-12) {
  method <- match.arg(method)
  idx <- which(open_p & is.finite(u))

  u_c <- u - stats::median(u[idx])
  scale <- if (method == "iqr") stats::IQR(u_c[idx]) else stats::sd(u_c[idx])
  if (!is.finite(scale) || scale < eps) scale <- 1

  list(u_scaled = u_c / scale, scale = scale)
}

update_effort_softmax <- function(
    eff_p, E_tot, u_swap, open_p,
    beta = 2.5, rho = 0.1, norm = "iqr",
    cap_frac = 0
) {
  eff_p[!open_p] <- 0

  nu <- normalize_util(u_swap, open_p, method = norm)
  u_star <- nu$u_scaled

  idx <- which(open_p & is.finite(u_swap))
  shares <- rep(0, length(eff_p))
  shares[idx] <- softmax_stable(beta * u_star[idx])

  eff_target <- E_tot * shares
  eff_new <- (1 - rho) * eff_p + rho * eff_target

  if (cap_frac > 0) {
    cap <- cap_frac * E_tot
    eff_new <- pmax(eff_p - cap, pmin(eff_p + cap, eff_new))
  }

  # project back to constraints
  eff_new[!open_p] <- 0
  s <- sum(eff_new)

  if (s > 0) {
    eff_new <- eff_new * (E_tot / s)
  } else {
    eff_new <- rep(0, length(eff_p))
    eff_new[open_p] <- E_tot / sum(open_p)
  }

  list(eff_p = eff_new, shares = shares, util_scale = nu$scale, u_scaled = u_star)
}

# ----------------------------
# 5) Inner-loop allocator (fast marginal + fast swap)
# ----------------------------

allocate_until_stable_patchwise <- function(
    eff_init_p, E_tot,
    v_mats, other_mats, b_mats, price_s, dt,
    c0, gamma, travel_p, open_p,
    beta = 1, rho = 0.1,
    max_iter = 40, tol = 1e-5,
    norm = "iqr",
    cap_frac = 0,
    record = TRUE
) {
  P <- length(eff_init_p)
  open_idx <- which(open_p)

  # init on simplex
  eff_p <- pmax(eff_init_p, 0)
  eff_p[!open_p] <- 0
  if (sum(eff_p) > 0) {
    eff_p <- eff_p * (E_tot / sum(eff_p))
  } else {
    eff_p[open_idx] <- E_tot / length(open_idx)
  }

  history <- if (record) vector("list", max_iter) else NULL

  for (k in seq_len(max_iter)) {
    eff_old <- eff_p

    # 1) fast analytic marginal
    m_p <- marginal_unconstrained_ragged_analytic_fast(
      eff_p = eff_p,
      v_mats = v_mats, other_mats = other_mats, b_mats = b_mats,
      price_s = price_s, dt = dt,
      c0 = c0, gamma = gamma, travel_p = travel_p,
      open_p = open_p
    )

    # 2) fast swap marginal
    u_swap <- swap_marginal_from_m_fast(eff_p, m_p, open_p)

    # 3) softmax update + inertia + projection
    upd <- update_effort_softmax(
      eff_p = eff_p, E_tot = E_tot,
      u_swap = u_swap, open_p = open_p,
      beta = beta, rho = rho,
      norm = norm, cap_frac = cap_frac
    )
    eff_p <- upd$eff_p

    rel_change <- sum(abs(eff_p - eff_old)) / E_tot

    if (record) {
      history[[k]] <- tibble::tibble(
        iter = k,
        patch = seq_len(P),
        effort = eff_p,
        m = m_p,
        u_swap = u_swap,
        open = open_p
      )
    }

    if (rel_change < tol) {
      if (record) history <- history[seq_len(k)]
      break
    }
  }

  list(eff_p = eff_p, history = history)
}

# ============================================================
# End of core suite
# ============================================================

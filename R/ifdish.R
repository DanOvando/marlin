# ============================================================
# marlin_softmax_patchwise_ragged_FAST.R
#
# End-to-end "fast swaps" applied:
#   1) Replace finite-difference marginal with analytic marginal:
#        marginal_unconstrained_ragged_analytic_fast()
#   2) Replace O(P^2) swap-marginal with vectorized O(P):
#        swap_marginal_from_m_fast()
#
# Everything else left as-is, except allocate_until_stable_patchwise()
# now calls the new fast marginal + fast swap.
#
# Notes:
# - Ragged ages by species => v_mats/other_mats/b_mats are lists of [P x A_s]
# - metier$vul_p_a is assumed precomputed and consistent with storage/fauna
#
# ============================================================

# ----------------------------
# 1) Precompute (unchanged)
# ----------------------------

precompute_baranov_inputs_softmax <- function(storage,
                                              fauna,
                                              fleets,
                                              target_fleet,
                                              E_exo = NULL,
                                              P) {
  n_species <- length(storage)
  n_fleet <- length(fleets)

  v_mats <- vector("list", n_species)
  other_mats <- vector("list", n_species)
  b_mats <- vector("list", n_species)
  price_s <- numeric(n_species)

  for (s in seq_len(n_species)) {
    b_pa <- storage[[s]]$b_p_a
    A_s <- ncol(b_pa)
    stopifnot(nrow(b_pa) == P)

    met_t <- fleets[[target_fleet]]$metiers[[s]]
    price_s[s] <- met_t$price

    v_pa <- met_t$vul_p_a
    stopifnot(!is.null(v_pa), nrow(v_pa) == P, ncol(v_pa) == A_s)

    other_fish_pa <- matrix(0, nrow = P, ncol = A_s)

    if (!is.null(E_exo) && n_fleet > 1) {
      for (fleet_idx in seq_len(n_fleet)) {
        if (fleet_idx == target_fleet) next

        met_f <- fleets[[fleet_idx]]$metiers[[s]]
        alpha_pa_f <- met_f$vul_p_a
        stopifnot(!is.null(alpha_pa_f), nrow(alpha_pa_f) == P, ncol(alpha_pa_f) == A_s)

        E_p_f <- as.numeric(E_exo[[fleet_idx]])
        stopifnot(length(E_p_f) == P)

        other_fish_pa <- other_fish_pa + alpha_pa_f * E_p_f
      }
    }

    m_a <- fauna[[s]]$m_at_age
    stopifnot(length(m_a) == A_s)
    m_pa <- matrix(m_a, nrow = P, ncol = A_s, byrow = TRUE)

    v_mats[[s]] <- v_pa
    other_mats[[s]] <- m_pa + other_fish_pa
    b_mats[[s]] <- b_pa
  }

  list(v_mats = v_mats, other_mats = other_mats, b_mats = b_mats, price_s = price_s)
}

# ----------------------------
# 2) FAST analytic marginal m[p] = d obj_p / d e_p
# ----------------------------

marginal_unconstrained_ragged_analytic_fast <- function(
    e,
    v_mats, other_mats, b_mats, price_s, dt,
    c0, gamma, travel_p,
    open_p,
    eps = 1e-12
) {
  # Analytic derivative of Baranov term with respect to effort.
  #
  # For each species s and age a at patch p:
  #   F = v * e
  #   O = other
  #   Z = F + O
  #   g(F) = (F/Z) * (1 - exp(-dt*Z))
  #
  # d/dF g(F) = (O/Z^2)*(1-exp(-dt*Z)) + (F/Z)*dt*exp(-dt*Z)
  # d/de g(F) = v * d/dF g(F)
  #
  # Marginal revenue at patch p:
  #   sum_s price_s[s] * sum_a b * d/de g(F)
  #
  # Marginal cost at patch p:
  #   c0 * (gamma*e^(gamma-1) + travel_p)
  #
  P <- length(e)
  e_use <- e
  e_use[!open_p] <- 0

  m_rev <- rep(0, P)

  for (s in seq_along(v_mats)) {
    v_pa     <- v_mats[[s]]      # [P x A_s]
    other_pa <- other_mats[[s]]  # [P x A_s]  (this is O)
    b_pa     <- b_mats[[s]]      # [P x A_s]

    # F and Z
    F_pa <- v_pa * e_use
    Z_pa <- F_pa + other_pa
    Z_pa <- pmax(Z_pa, eps)

    exp_term <- exp(-dt * Z_pa)
    one_minus <- 1 - exp_term

    # d/dF g(F)
    term1 <- (other_pa / (Z_pa^2)) * one_minus
    term2 <- (F_pa / Z_pa) * (dt * exp_term)
    dgdF <- term1 + term2

    # d/de g(F) = v * dgdF
    dgdE <- v_pa * dgdF

    # add price-weighted biomass contribution (rowSums is fast)
    m_rev <- m_rev + price_s[s] * rowSums(b_pa * dgdE)
  }

  # marginal cost
  # gamma*e^(gamma-1) is well-defined at e=0 for gamma>=1 (0^(>0)=0); use pmax just in case
  e_pos <- pmax(e_use, 0)
  m_cost <- c0 * (gamma * (e_pos^(gamma - 1)) + travel_p)

  m <- m_rev - m_cost
  m[!open_p] <- NA_real_
  m
}

# ----------------------------
# 3) FAST swap marginal u_swap[p] in O(P)
# ----------------------------

swap_marginal_from_m_fast <- function(e, m, open_p, eps = 1e-12) {
  # Budget-neutral "proportional donor" swap marginal:
  # donors are all other open patches; donor funding shares proportional to donor effort.
  #
  # For open patch p:
  #   donor_mean_m(p) = (sum_j e_j m_j - e_p m_p) / (sum_j e_j - e_p)
  #   u_swap[p] = m_p - donor_mean_m(p)
  #
  u_swap <- rep(NA_real_, length(e))

  idx <- which(open_p & is.finite(m))
  if (length(idx) <= 1) return(u_swap)

  e_i <- e[idx]
  m_i <- m[idx]

  E <- sum(e_i)
  S <- sum(e_i * m_i)

  denom <- pmax(E - e_i, eps)
  donor_mean <- (S - e_i * m_i) / denom

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
    e, E_tot, u_swap, open_p,
    beta = 2.5, rho = 0.1, norm = "iqr",
    cap_frac = 0
) {
  e[!open_p] <- 0

  nu <- normalize_util(u_swap, open_p, method = norm)
  u_star <- nu$u_scaled

  idx <- which(open_p & is.finite(u_swap))
  shares <- rep(0, length(e))
  shares[idx] <- softmax_stable(beta * u_star[idx])

  e_target <- E_tot * shares
  e_new <- (1 - rho) * e + rho * e_target

  if (cap_frac > 0) {
    cap <- cap_frac * E_tot
    e_new <- pmax(e - cap, pmin(e + cap, e_new))
  }

  e_new[!open_p] <- 0
  s <- sum(e_new)

  if (s > 0) {
    e_new <- e_new * (E_tot / s)
  } else {
    e_new <- rep(0, length(e))
    e_new[open_p] <- E_tot / sum(open_p)
  }

  list(e = e_new, shares = shares, util_scale = nu$scale, u_scaled = u_star)
}

# ----------------------------
# 5) Inner-loop allocator (UPDATED to use fast marginal + fast swap)
# ----------------------------

allocate_until_stable_patchwise <- function(
    e_init, E_tot,
    v_mats, other_mats, b_mats, price_s, dt,
    c0, gamma, travel_p, open_p,
    beta = 1, rho = 0.1,
    max_iter = 40, tol = 1e-5,
    norm = "iqr",
    cap_frac = 0,
    record = TRUE
) {
  P <- length(e_init)
  open_idx <- which(open_p)

  # init on simplex
  e <- pmax(e_init, 0)
  e[!open_p] <- 0
  if (sum(e) > 0) {
    e <- e * (E_tot / sum(e))
  } else {
    e[open_idx] <- E_tot / length(open_idx)
  }

  history <- if (record) vector("list", max_iter) else NULL

  for (k in seq_len(max_iter)) {
    e_old <- e

    # 1) fast analytic marginal
    m <- marginal_unconstrained_ragged_analytic_fast(
      e = e,
      v_mats = v_mats, other_mats = other_mats, b_mats = b_mats,
      price_s = price_s, dt = dt,
      c0 = c0, gamma = gamma, travel_p = travel_p,
      open_p = open_p
    )

    # 2) fast swap marginal
    u_swap <- swap_marginal_from_m_fast(e, m, open_p)

    # 3) softmax update + inertia + projection
    upd <- update_effort_softmax(
      e = e, E_tot = E_tot,
      u_swap = u_swap, open_p = open_p,
      beta = beta, rho = rho,
      norm = norm, cap_frac = cap_frac
    )
    e <- upd$e

    rel_change <- sum(abs(e - e_old)) / E_tot

    if (record) {
      history[[k]] <- tibble::tibble(
        iter = k,
        patch = seq_len(P),
        effort = e,
        m = m,
        u_swap = u_swap,
        open = open_p
      )
    }

    if (rel_change < tol) {
      if (record) history <- history[seq_len(k)]
      break
    }
  }

  list(e = e, history = history)
}

# ============================================================
# End of core suite
# ============================================================

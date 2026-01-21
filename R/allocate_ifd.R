#' Allocate Ideal Free Distribution
#'
#' @param Etot_target toatl effort
#' @param storage a place for biomass
#' @param fauna the fauna objects
#' @param fleets the fleet object
#' @param target_fleet integer indicating which fleet is being allocated
#' @param E_exo total exogenous effort in each patch
#' @param fleet_fishable list of fishing grounds available to each fleet
#' @param n_inner optim stufff
#' @param tol_outer optim stufff
#' @param active_eps optim stufff
#' @param lambda_init optim stufff
#'
#' @returns a list with the effort by patch of the target fleet and diagnostics
#' @export
#'
#' @examples
#
#'
# ----------------------------
# Baranov helper (single source of truth)
# ----------------------------
baranov_catch_and_dcatch <- function(fishing_mort, other_mort, biomass, time_step) {
  # fishing_mort: matrix [P x A] for the focal fleet/species
  # other_mort  : matrix [P x A] = M + other fleets' fishing mortality
  # biomass     : matrix [P x A]
  # time_step   : scalar > 0
  stopifnot(is.matrix(fishing_mort), is.matrix(other_mort), is.matrix(biomass))
  stopifnot(all(dim(fishing_mort) == dim(other_mort)), all(dim(biomass) == dim(other_mort)))
  stopifnot(is.numeric(time_step), length(time_step) == 1L, time_step > 0)

  total_mort <- other_mort + fishing_mort
  # guard against divide-by-zero if total_mort hits 0 exactly
  z_safe <- pmax(total_mort, 1e-15)

  surv_term <- exp(-time_step * z_safe)
  harvest_term <- 1 - surv_term

  catch <- (fishing_mort / z_safe) * harvest_term * biomass

  # dCatch/dF (F = fishing_mort), with other_mort held fixed
  # dC/dF = [other/z^2 * (1-exp(-dt z)) + (F/z) * (dt exp(-dt z))] * B
  dcatch_dF <- ((other_mort / (z_safe^2)) * harvest_term +
                  (fishing_mort / z_safe) * (time_step * surv_term)) * biomass

  list(
    catch = catch,
    dcatch_dF = dcatch_dF,
    total_mort = total_mort
  )
}

# ----------------------------
# Allocator
# ----------------------------
allocate_ifd <- function(
    Etot_target,
    storage,
    fauna,
    fleets,
    target_fleet = 1L,
    E_exo = NULL,
    fleet_fishable = NULL,
    time_step = 1,     # NEW: scales instantaneous rates inside Baranov
    n_inner = 40,
    tol_outer = 1e-10,
    active_eps = 1e-12,
    lambda_init = NULL
) {

  stopifnot(is.numeric(Etot_target), length(Etot_target) == 1L, Etot_target >= 0)
  stopifnot(is.numeric(time_step), length(time_step) == 1L, time_step > 0)
  stopifnot(is.list(storage), is.list(fauna), is.list(fleets))

  S <- length(storage)
  NF <- length(fleets)
  stopifnot(length(fauna) == S)
  stopifnot(target_fleet >= 1L, target_fleet <= NF)

  # ---- infer P and A_s from storage; validate fauna R6 ----
  stopifnot(!is.null(storage[[1]]$b_p_a), is.matrix(storage[[1]]$b_p_a))
  P <- nrow(storage[[1]]$b_p_a)

  A_s <- integer(S)
  for (s in seq_len(S)) {
    stopifnot(!is.null(storage[[s]]$b_p_a))
    biomass <- storage[[s]]$b_p_a
    stopifnot(is.matrix(biomass), nrow(biomass) == P, ncol(biomass) >= 1L)
    A_s[s] <- ncol(biomass)

    stopifnot(inherits(fauna[[s]], "R6"))
    nat_mort <- fauna[[s]]$m_at_age
    stopifnot(is.numeric(nat_mort), length(nat_mort) == A_s[s])
  }

  # ---- validate fleets / metiers ----
  for (fleet_idx in seq_len(NF)) {
    fl <- fleets[[fleet_idx]]
    stopifnot(is.numeric(fl$cost_per_patch), length(fl$cost_per_patch) == P)
    stopifnot(is.numeric(fl$cost_per_unit_effort),
              length(fl$cost_per_unit_effort) == 1L, fl$cost_per_unit_effort >= 0)
    stopifnot(is.numeric(fl$effort_cost_exponent),
              length(fl$effort_cost_exponent) == 1L, fl$effort_cost_exponent >= 1)

    stopifnot(is.list(fl$metiers), length(fl$metiers) == S)
    for (s in seq_len(S)) {
      met <- fl$metiers[[s]]
      stopifnot(is.numeric(met$sel_at_age), length(met$sel_at_age) == A_s[s])
      stopifnot(is.numeric(met$price), length(met$price) == 1L)
      stopifnot(is.numeric(met$spatial_catchability), length(met$spatial_catchability) == P)
    }
  }

  # ---- fleet fishable masks ----
  if (is.null(fleet_fishable)) {
    fleet_fishable <- replicate(NF, rep(1L, P), simplify = FALSE)
  } else {
    stopifnot(is.list(fleet_fishable), length(fleet_fishable) == NF)
    for (fleet_idx in seq_len(NF)) {
      v <- fleet_fishable[[fleet_idx]]
      stopifnot(length(v) == P)
      if (is.logical(v)) v <- as.integer(v)
      stopifnot(all(v %in% c(0L, 1L)))
      fleet_fishable[[fleet_idx]] <- as.integer(v)
    }
  }
  fishable_target <- fleet_fishable[[target_fleet]] == 1L

  if (!any(fishable_target)) {
    return(list(
      E_target = numeric(P),
      lambda = NA_real_,
      marginal_profit = rep(NA_real_, P),
      revenue_p = rep(0, P),
      cost_p_total = rep(0, P),
      profit_p = rep(0, P),
      diagnostic = list(note = "No fishable patches for target fleet; forced E=0", n_active = 0L)
    ))
  }

  # ---- exogenous effort ----
  if (is.null(E_exo)) {
    E_exo <- replicate(NF, numeric(P), simplify = FALSE)
  } else {
    stopifnot(is.list(E_exo), length(E_exo) == NF)
    for (fleet_idx in seq_len(NF)) {
      stopifnot(is.numeric(E_exo[[fleet_idx]]), length(E_exo[[fleet_idx]]) == P)
    }
  }

  # Enforce closures on exogenous effort (safe default)
  for (fleet_idx in seq_len(NF)) {
    closed <- fleet_fishable[[fleet_idx]] == 0L
    if (any(closed)) E_exo[[fleet_idx]][closed] <- 0
  }

  if (Etot_target == 0) {
    return(list(
      E_target = numeric(P),
      lambda = NA_real_,
      marginal_profit = rep(NA_real_, P),
      revenue_p = rep(0, P),
      cost_p_total = rep(0, P),
      profit_p = rep(0, P),
      diagnostic = list(note = "Etot_target = 0", n_active = 0L, n_fishable = sum(fishable_target))
    ))
  }

  # ---- target fleet cost params ----
  cost_patch <- fleets[[target_fleet]]$cost_per_patch
  c0 <- fleets[[target_fleet]]$cost_per_unit_effort
  gamma <- fleets[[target_fleet]]$effort_cost_exponent

  # ----------------------------
  # Precompute per-species pieces
  # alpha_t[p,a] = q_{t,s,p} * sel_{t,s,a}
  # other_mort[p,a] = M[a] + sum_{fleet!=t} q_{fleet,s,p}*sel_{fleet,s,a}*E_exo_fleet[p]
  # ----------------------------
  alpha_list <- vector("list", S)  # P x A_s
  other_mort_list <- vector("list", S)  # P x A_s
  biomass_list <- vector("list", S)  # P x A_s
  price_vec <- numeric(S)

  for (s in seq_len(S)) {
    biomass <- storage[[s]]$b_p_a
    nat_mort <- fauna[[s]]$m_at_age

    met_t <- fleets[[target_fleet]]$metiers[[s]]
    alpha_t <- met_t$spatial_catchability %o% met_t$sel_at_age
    price_vec[s] <- met_t$price

    # other fleets' instantaneous fishing mortality
    other_fleet_mort <- matrix(0, nrow = P, ncol = A_s[s])
    for (fleet_idx in seq_len(NF)) {
      if (fleet_idx == target_fleet) next
      met_f <- fleets[[fleet_idx]]$metiers[[s]]
      other_fleet_mort <- other_fleet_mort +
        (met_f$spatial_catchability %o% met_f$sel_at_age) * E_exo[[fleet_idx]]
    }

    other_mort <- matrix(nat_mort, nrow = P, ncol = A_s[s], byrow = TRUE) + other_fleet_mort

    alpha_list[[s]] <- alpha_t
    other_mort_list[[s]] <- other_mort
    biomass_list[[s]] <- biomass
  }

  # ---- marginal revenue wrt E_target[p] ----
  d_rev_dE <- function(E_p) {
    out <- numeric(P)

    for (s in seq_len(S)) {
      alpha <- alpha_list[[s]]
      other_mort <- other_mort_list[[s]]
      biomass <- biomass_list[[s]]
      price <- price_vec[s]

      fishing_mort <- alpha * E_p

      res <- baranov_catch_and_dcatch(
        fishing_mort = fishing_mort,
        other_mort = other_mort,
        biomass = biomass,
        time_step = time_step
      )

      # dF/dE = alpha; marginal revenue = price * sum_a alpha * dcatch_dF
      out <- out + rowSums(alpha * res$dcatch_dF) * price
    }

    out
  }

  d_cost_dE <- function(E_p) {
    c0 * (cost_patch + gamma * E_p^(gamma - 1))
  }

  d_pi_dE <- function(E_p) d_rev_dE(E_p) - d_cost_dE(E_p)

  # mp at zero; closed patches forced inactive
  mp0 <- d_pi_dE(rep(0, P))
  mp0[!fishable_target] <- -Inf

  solve_E_given_lambda <- function(lambda) {
    active <- fishable_target & (mp0 > lambda)
    E <- numeric(P)
    if (!any(active)) return(E)

    lo <- numeric(P)
    hi <- rep(Etot_target, P)
    hi[!fishable_target] <- 0

    for (k in seq_len(n_inner)) {
      mid <- 0.5 * (lo + hi)
      mp_mid <- d_pi_dE(mid)
      mp_mid[!fishable_target] <- -Inf

      lo[active & mp_mid > lambda] <- mid[active & mp_mid > lambda]
      hi[active & mp_mid <= lambda] <- mid[active & mp_mid <= lambda]
    }

    E[active] <- hi[active]
    E[!fishable_target] <- 0
    E
  }

  F_lambda <- function(lambda) sum(solve_E_given_lambda(lambda)) - Etot_target

  # ---- bracketing for uniroot ----
  mp0_open <- mp0[fishable_target]
  lambda_hi <- max(mp0_open) + 1e-9

  mp_at_hi <- d_pi_dE(rep(Etot_target, P))
  lambda_lo <- min(mp_at_hi[fishable_target]) - abs(min(mp_at_hi[fishable_target])) - 1

  if (!is.null(lambda_init) && is.finite(lambda_init)) {
    width <- 0.1 * max(1, abs(lambda_init))
    lo <- lambda_init - width
    hi <- lambda_init + width
    f_lo <- F_lambda(lo)
    f_hi <- F_lambda(hi)
    step <- 0L
    while (f_lo * f_hi > 0 && step < 60L) {
      width <- width * 2
      lo <- lambda_init - width
      hi <- lambda_init + width
      f_lo <- F_lambda(lo)
      f_hi <- F_lambda(hi)
      step <- step + 1L
    }
    if (f_lo * f_hi <= 0) {
      lambda_lo <- lo
      lambda_hi <- hi
    }
  }

  lambda_star <- stats::uniroot(F_lambda, lower = lambda_lo, upper = lambda_hi, tol = tol_outer)$root
  E_star <- solve_E_given_lambda(lambda_star)

  # ----------------------------
  # Accounting (target fleet only): revenue/cost/profit by patch
  # ----------------------------
  revenue_p <- numeric(P)
  for (s in seq_len(S)) {
    alpha <- alpha_list[[s]]
    other_mort <- other_mort_list[[s]]
    biomass <- biomass_list[[s]]
    price <- price_vec[s]

    fishing_mort <- alpha * E_star
    res <- baranov_catch_and_dcatch(
      fishing_mort = fishing_mort,
      other_mort = other_mort,
      biomass = biomass,
      time_step = time_step
    )

    revenue_p <- revenue_p + rowSums(res$catch) * price
  }
  revenue_p[!fishable_target] <- 0

  cost_p_total <- c0 * (cost_patch * E_star + E_star^gamma)
  profit_p <- revenue_p - cost_p_total

  # marginal profit at solution
  mp <- d_pi_dE(E_star)
  mp[!fishable_target] <- NA_real_
  active <- fishable_target & (E_star > active_eps)

  # extra diagnostics
  sum_E_closed_target <- sum(E_star[!fishable_target])
  max_E_closed_target <- if (any(!fishable_target)) max(E_star[!fishable_target]) else 0
  sum_Eexo_closed_by_fleet <- vapply(seq_len(NF), function(fleet_idx) {
    closed <- fleet_fishable[[fleet_idx]] == 0L
    sum(E_exo[[fleet_idx]][closed])
  }, numeric(1))

  diagnostic <- list(
    time_step = time_step,
    lambda = lambda_star,
    sum_E = sum(E_star),
    n_fishable = sum(fishable_target),
    n_closed = sum(!fishable_target),
    n_active = sum(active),
    mp_active_range = if (any(active)) range(mp[active]) else c(NA_real_, NA_real_),
    mp_minus_lambda_active_range = if (any(active)) range(mp[active] - lambda_star) else c(NA_real_, NA_real_),
    sum_E_closed_target = sum_E_closed_target,
    max_E_closed_target = max_E_closed_target,
    sum_Eexo_closed_by_fleet = sum_Eexo_closed_by_fleet
  )

  list(
    E_target = E_star,
    lambda = lambda_star,
    marginal_profit = mp,
    revenue_p = revenue_p,
    cost_p_total = cost_p_total,
    profit_p = profit_p,
    total_revenue = sum(revenue_p),
    total_cost = sum(cost_p_total),
    total_profit = sum(profit_p),
    diagnostic = diagnostic
  )
}


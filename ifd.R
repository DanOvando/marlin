# ============================================================
# WORKING VERSION (robust in tiny-biomass regime)
#
# Key fix vs your current version:
#   When include_costs == FALSE, we treat the KKT objective as potentially
#   numerically flat at small biomass. We detect flatness using BOTH:
#     1) sd(mp0_open_finite) <= flat_tol_sd   (your proven fix; default 1e-6)
#     2) maxabs(mp0_open_finite) <= flat_tol_abs (extra safety)
#   If flat, we SKIP uniroot entirely and allocate uniformly across fishable.
#
# Also:
# - C++ compiles on Apple clang 17: use as<NumericMatrix>(list[s])
# - Uses expm1 stability
# - Avoids using F as a variable name in R
# ============================================================

library(Rcpp)

Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <Rmath.h>      // expm1, R_FINITE
#include <algorithm>    // std::fill
#include <cmath>        // std::exp
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector cpp_d_rev_dE_baranov(
    const NumericVector& effort_p,
    const List& alpha_mats,
    const List& other_mort_mats,
    const List& biomass_mats,
    const NumericVector& price_s,
    const double time_step
) {
  const int n_species = alpha_mats.size();
  const int n_patch = effort_p.size();
  const double dt = time_step;

  NumericVector out(n_patch);
  std::fill(out.begin(), out.end(), 0.0);

  for (int s = 0; s < n_species; ++s) {
    NumericMatrix alpha = as<NumericMatrix>(alpha_mats[s]);
    NumericMatrix other = as<NumericMatrix>(other_mort_mats[s]);
    NumericMatrix biom  = as<NumericMatrix>(biomass_mats[s]);

    const int n_age = alpha.ncol();
    const double price = price_s[s];

    for (int a = 0; a < n_age; ++a) {
      for (int p = 0; p < n_patch; ++p) {
        const double alpha_pa  = alpha(p, a);
        const double fish_mort = alpha_pa * effort_p[p];
        const double other_m   = other(p, a);
        const double bio       = biom(p, a);

        double z = other_m + fish_mort;
        if (!R_FINITE(z)) z = 0.0;

        const double z_safe = (z > 1e-15) ? z : 1e-15;
        const double x = dt * z_safe;

        const double harvest = -expm1(-x);
        const double surv    = std::exp(-x);

        double dcatch_dmort =
          ((other_m / (z_safe * z_safe)) * harvest +
           (fish_mort / z_safe) * (dt * surv)) * bio;

        if (!R_FINITE(dcatch_dmort)) dcatch_dmort = 0.0;

        out[p] += alpha_pa * dcatch_dmort * price;
      }
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericVector cpp_revenue_baranov(
    const NumericVector& effort_p,
    const List& alpha_mats,
    const List& other_mort_mats,
    const List& biomass_mats,
    const NumericVector& price_s,
    const double time_step
) {
  const int n_species = alpha_mats.size();
  const int n_patch = effort_p.size();
  const double dt = time_step;

  NumericVector revenue_p(n_patch);
  std::fill(revenue_p.begin(), revenue_p.end(), 0.0);

  for (int s = 0; s < n_species; ++s) {
    NumericMatrix alpha = as<NumericMatrix>(alpha_mats[s]);
    NumericMatrix other = as<NumericMatrix>(other_mort_mats[s]);
    NumericMatrix biom  = as<NumericMatrix>(biomass_mats[s]);

    const int n_age = alpha.ncol();
    const double price = price_s[s];

    for (int a = 0; a < n_age; ++a) {
      for (int p = 0; p < n_patch; ++p) {
        const double alpha_pa  = alpha(p, a);
        const double fish_mort = alpha_pa * effort_p[p];
        const double other_m   = other(p, a);
        const double bio       = biom(p, a);

        double z = other_m + fish_mort;
        if (!R_FINITE(z)) z = 0.0;

        const double z_safe = (z > 1e-15) ? z : 1e-15;
        const double x = dt * z_safe;

        const double harvest = -expm1(-x);

        double catch_pa = (fish_mort / z_safe) * harvest * bio;
        if (!R_FINITE(catch_pa)) catch_pa = 0.0;

        revenue_p[p] += catch_pa * price;
      }
    }
  }

  return revenue_p;
}
', verbose = TRUE)

# ----------------------------
# Precompute alpha/other_mort/biomass mats
# ----------------------------
precompute_baranov_inputs <- function(storage, fauna, fleets, target_fleet, E_exo, P) {
  n_species <- length(storage)
  n_fleet <- length(fleets)

  alpha_mats <- vector("list", n_species)
  other_mort_mats <- vector("list", n_species)
  biomass_mats <- vector("list", n_species)
  price_s <- numeric(n_species)

  for (s in seq_len(n_species)) {
    biom <- storage[[s]]$b_p_a
    n_age <- ncol(biom)

    met_t <- fleets[[target_fleet]]$metiers[[s]]
    price_s[s] <- met_t$price

    alpha_mat <- met_t$spatial_catchability %o% met_t$sel_at_age

    other_fish <- matrix(0, nrow = P, ncol = n_age)
    if (n_fleet > 1) {
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

# ----------------------------
# Allocator (robust flat handling)
# ----------------------------
allocate_ifd_fixed_Etot_storage_fauna_fleets_rcpp <- function(
    Etot_target,
    storage,
    fauna,
    fleets,
    target_fleet = 1L,
    E_exo = NULL,
    fleet_fishable = NULL,
    time_step = 1,
    include_costs = TRUE,
    active_tol = 1e-14,
    # IMPORTANT: default now matches your empirical fix
    flat_tol_sd = 1e-6,
    flat_tol_abs = 1e-14,
    n_inner = 30,
    tol_outer = 1e-8,
    active_eps = 1e-12,
    lambda_init = NULL,
    validate_inputs = FALSE
) {

  if (validate_inputs) {
    stopifnot(is.numeric(Etot_target), length(Etot_target) == 1L, Etot_target >= 0)
    stopifnot(is.numeric(time_step), length(time_step) == 1L, time_step > 0)
    stopifnot(is.logical(include_costs), length(include_costs) == 1L)
    stopifnot(is.list(storage), is.list(fauna), is.list(fleets))
  }

  n_fleet <- length(fleets)
  P <- nrow(storage[[1]]$b_p_a)

  # fishable masks
  if (is.null(fleet_fishable)) {
    fleet_fishable <- replicate(n_fleet, rep(1L, P), simplify = FALSE)
  }
  for (fleet_idx in seq_len(n_fleet)) {
    v <- fleet_fishable[[fleet_idx]]
    if (is.logical(v)) v <- as.integer(v)
    fleet_fishable[[fleet_idx]] <- as.integer(v)
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

  # exogenous effort
  if (is.null(E_exo)) {
    E_exo <- replicate(n_fleet, numeric(P), simplify = FALSE)
  }

  # enforce closures on exogenous effort
  for (fleet_idx in seq_len(n_fleet)) {
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

  # target costs
  cost_patch <- fleets[[target_fleet]]$cost_per_patch
  c0 <- fleets[[target_fleet]]$cost_per_unit_effort
  gamma <- fleets[[target_fleet]]$effort_cost_exponent

  pre <- precompute_baranov_inputs(storage, fauna, fleets, target_fleet, E_exo, P)
  alpha_mats <- pre$alpha_mats
  other_mort_mats <- pre$other_mort_mats
  biomass_mats <- pre$biomass_mats
  price_s <- pre$price_s

  d_rev_dE <- function(E_p) {
    out <- cpp_d_rev_dE_baranov(E_p, alpha_mats, other_mort_mats, biomass_mats, price_s, time_step)
    out[!is.finite(out)] <- 0
    out
  }

  d_cost_dE <- function(E_p) {
    if (!include_costs) return(rep(0, length(E_p)))
    if (gamma == 1) {
      c0 * (cost_patch + gamma)
    } else {
      c0 * (cost_patch + gamma * (E_p^(gamma - 1)))
    }
  }

  d_obj_dE <- function(E_p) d_rev_dE(E_p) - d_cost_dE(E_p)

  mp0 <- d_obj_dE(rep(0, P))
  mp0[!is.finite(mp0)] <- -Inf
  mp0[!fishable_target] <- -Inf

  mp0_open <- mp0[fishable_target]
  mp0_open_finite <- mp0_open[is.finite(mp0_open)]
  mp_sd <- if (length(mp0_open_finite) >= 2) stats::sd(mp0_open_finite) else 0
  mp_maxabs <- if (length(mp0_open_finite) >= 1) max(abs(mp0_open_finite)) else 0

  # ---- Robust flat fallback (REVENUE ONLY) ----
  # This is the key behavior you observed: flat_tol_sd around 1e-6 fixes it.
  if (!include_costs && (mp_sd <= flat_tol_sd || mp_maxabs <= flat_tol_abs)) {
    E_star <- numeric(P)
    E_star[fishable_target] <- Etot_target / sum(fishable_target)

    revenue_p <- cpp_revenue_baranov(E_star, alpha_mats, other_mort_mats, biomass_mats, price_s, time_step)
    revenue_p[!fishable_target] <- 0

    cost_p_total <- c0 * (cost_patch * E_star + E_star^gamma)
    profit_p <- revenue_p - cost_p_total

    mp_sol <- d_obj_dE(E_star)
    mp_sol[!is.finite(mp_sol)] <- NA_real_
    mp_sol[!fishable_target] <- NA_real_

    diagnostic <- list(
      time_step = time_step,
      include_costs = include_costs,
      flat_fallback_used = TRUE,
      flat_fallback_mp_sd = mp_sd,
      flat_fallback_mp_maxabs = mp_maxabs,
      flat_tol_sd = flat_tol_sd,
      flat_tol_abs = flat_tol_abs,
      lambda = 0,
      sum_E = sum(E_star),
      n_fishable = sum(fishable_target),
      n_active = sum(fishable_target & (E_star > active_eps))
    )

    return(list(
      E_target = E_star,
      lambda = 0,
      marginal_profit = mp_sol,
      revenue_p = revenue_p,
      cost_p_total = cost_p_total,
      profit_p = profit_p,
      total_revenue = sum(revenue_p),
      total_cost = sum(cost_p_total),
      total_profit = sum(profit_p),
      diagnostic = diagnostic
    ))
  }

  # ---- KKT solver path ----
  solve_E_given_lambda <- function(lambda) {
    active <- fishable_target & (mp0 > (lambda + active_tol))
    effort <- numeric(P)
    if (!any(active)) return(effort)

    lo <- numeric(P)
    hi <- rep(Etot_target, P)
    hi[!fishable_target] <- 0

    for (k in seq_len(n_inner)) {
      mid <- 0.5 * (lo + hi)
      mp_mid <- d_obj_dE(mid)
      mp_mid[!is.finite(mp_mid)] <- -Inf
      mp_mid[!fishable_target] <- -Inf

      go_up <- mp_mid > (lambda + active_tol)
      lo[active & go_up] <- mid[active & go_up]
      hi[active & !go_up] <- mid[active & !go_up]
    }

    effort[active] <- hi[active]
    effort[!fishable_target] <- 0
    effort
  }

  F_lambda <- function(lambda) {
    effort_try <- solve_E_given_lambda(lambda)
    s <- sum(effort_try)
    if (!is.finite(s)) return(Inf)
    s - Etot_target
  }

  # scale-aware bracket bump (avoid hard 1e-9 when mp scale is tiny)
  mp0_max <- max(mp0_open_finite, 0)
  bump <- max(1e-12, 1e-6 * mp0_max)
  lambda_hi <- mp0_max + bump

  mp_at_hi <- d_obj_dE(rep(Etot_target, P))
  mp_at_hi[!is.finite(mp_at_hi)] <- -Inf
  lo_core <- mp_at_hi[fishable_target]
  lo_core <- lo_core[is.finite(lo_core)]
  if (length(lo_core) == 0) lo_core <- -1
  lambda_lo <- min(lo_core) - abs(min(lo_core)) - 1

  if (!is.null(lambda_init) && is.finite(lambda_init)) {
    width <- 0.1 * max(1, abs(lambda_init))
    lo <- lambda_init - width
    hi <- lambda_init + width
    f_lo <- F_lambda(lo)
    f_hi <- F_lambda(hi)
    step <- 0L
    while (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0 && step < 60L) {
      width <- width * 2
      lo <- lambda_init - width
      hi <- lambda_init + width
      f_lo <- F_lambda(lo)
      f_hi <- F_lambda(hi)
      step <- step + 1L
    }
    if (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi <= 0) {
      lambda_lo <- lo
      lambda_hi <- hi
    }
  }

  lambda_star <- stats::uniroot(F_lambda, lower = lambda_lo, upper = lambda_hi, tol = tol_outer)$root
  E_star <- solve_E_given_lambda(lambda_star)

  revenue_p <- cpp_revenue_baranov(E_star, alpha_mats, other_mort_mats, biomass_mats, price_s, time_step)
  revenue_p[!fishable_target] <- 0

  cost_p_total <- c0 * (cost_patch * E_star + E_star^gamma)
  profit_p <- revenue_p - cost_p_total

  mp_sol <- d_obj_dE(E_star)
  mp_sol[!is.finite(mp_sol)] <- NA_real_
  mp_sol[!fishable_target] <- NA_real_

  diagnostic <- list(
    time_step = time_step,
    include_costs = include_costs,
    flat_fallback_used = FALSE,
    mp0_open_sd = mp_sd,
    mp0_open_maxabs = mp_maxabs,
    n_inner = n_inner,
    tol_outer = tol_outer,
    lambda = lambda_star,
    sum_E = sum(E_star),
    n_fishable = sum(fishable_target),
    n_active = sum(fishable_target & (E_star > active_eps))
  )

  list(
    E_target = E_star,
    lambda = lambda_star,
    marginal_profit = mp_sol,
    revenue_p = revenue_p,
    cost_p_total = cost_p_total,
    profit_p = profit_p,
    total_revenue = sum(revenue_p),
    total_cost = sum(cost_p_total),
    total_profit = sum(profit_p),
    diagnostic = diagnostic
  )
}

# ============================================================
# WORKED EXAMPLE
# ============================================================

library(R6)

FaunaSpecies <- R6Class(
  "FaunaSpecies",
  public = list(
    m_at_age = NULL,
    initialize = function(m_at_age) {
      stopifnot(is.numeric(m_at_age), length(m_at_age) >= 1)
      self$m_at_age <- m_at_age
    }
  )
)

set.seed(123)

P <- 20^2
S <- 3
n_fleet <- 3
A_vec <- c(10, 42, 18)

storage <- vector("list", S)
for (s in seq_len(S)) {
  n_age <- A_vec[s]
  b_p_a <- matrix(rlnorm(P * n_age, meanlog = 0, sdlog = 0.7), nrow = P, ncol = n_age)
  b_p_a <- b_p_a * (seq(0.7, 1.3, length.out = P) %o% rep(1, n_age))
  storage[[s]] <- list(b_p_a = b_p_a)
}

fauna <- lapply(seq_len(S), function(s) {
  FaunaSpecies$new(m_at_age = runif(A_vec[s], 0.05, 0.35))
})

make_sel <- function(n_age) {
  peak <- sample(3:(n_age - 2), 1)
  sel <- exp(-0.5 * ((1:n_age) - peak)^2 / (2.2^2))
  sel / max(sel)
}

make_fleet <- function() {
  metiers <- vector("list", S)
  for (s in seq_len(S)) {
    n_age <- A_vec[s]
    metiers[[s]] <- list(
      sel_at_age = make_sel(n_age),
      price = runif(1, 0.5, 3.0),
      spatial_catchability = runif(P, 0.02, 0.25)
    )
  }
  list(
    cost_per_patch = runif(P, 0.2, 2.0),
    cost_per_unit_effort = runif(1, 0.03, 0.12),
    effort_cost_exponent = runif(1, 1.2, 1.7),
    metiers = metiers
  )
}

fleets <- lapply(seq_len(n_fleet), function(i) make_fleet())

fleet_fishable <- replicate(n_fleet, rep(1L, P), simplify = FALSE)
set.seed(99)
fleet_fishable[[1]][sample.int(P, round(0.30 * P))] <- 0L
fleet_fishable[[2]][sample.int(P, round(0.10 * P))] <- 0L
fleet_fishable[[3]] <- rep(1L, P)

E_exo <- replicate(n_fleet, numeric(P), simplify = FALSE)
E_exo[[1]] <- rep(0, P)
E_exo[[2]] <- rexp(P, rate = 1/10)
E_exo[[3]] <- rexp(P, rate = 1/6)

Etot_target <- 600
time_step <- 1
include_costs <- FALSE

cat("RUN 1 (normal biomass)\\n")
out1 <- allocate_ifd_fixed_Etot_storage_fauna_fleets_rcpp(
  Etot_target = Etot_target,
  storage = storage,
  fauna = fauna,
  fleets = fleets,
  target_fleet = 1L,
  E_exo = E_exo,
  fleet_fishable = fleet_fishable,
  time_step = time_step,
  include_costs = include_costs
)
cat("Sum E_target:", sum(out1$E_target), " target:", Etot_target, "\\n")
cat("Lambda:", out1$lambda, "\\n")
print(out1$diagnostic)

cat("\\nRUN 2 (tiny biomass)\\n")
storage_tiny <- storage
for (s in seq_len(S)) storage_tiny[[s]]$b_p_a <- storage_tiny[[s]]$b_p_a * 1e-11

out2 <- allocate_ifd_fixed_Etot_storage_fauna_fleets_rcpp(
  Etot_target = Etot_target,
  storage = storage_tiny,
  fauna = fauna,
  fleets = fleets,
  target_fleet = 1L,
  E_exo = E_exo,
  fleet_fishable = fleet_fishable,
  time_step = time_step,
  include_costs = include_costs
)
cat("Sum E_target:", sum(out2$E_target), " target:", Etot_target, "\\n")
cat("Lambda:", out2$lambda, "\\n")
print(out2$diagnostic)

cat("Any non-finite E_target?:", any(!is.finite(out2$E_target)), "\\n")

# ============================================================
# WORKED EXAMPLE: structured vs flat vs near-flat (not quite)
# - all patches open
# - E_exo = NULL
# - includes diagnostic plots for spatial allocation + objective surfaces
#
# NOTE: Requires:
#   - cpp_allocate_ifd_kkt_fullsolve_fast() already compiled/available
#   - precompute_baranov_inputs() available (your version with E_exo=NULL support)
# ============================================================

library(R6)
library(ggplot2)

stopifnot(exists("cpp_allocate_ifd_kkt_fullsolve_fast"))
stopifnot(exists("precompute_baranov_inputs"))

# ----------------------------
# minimal R6 fauna class for example
# ----------------------------
FaunaSpecies <- R6Class(
  "FaunaSpecies",
  public = list(
    m_at_age = NULL,
    initialize = function(m_at_age) self$m_at_age <- m_at_age
  )
)

# ----------------------------
# shared dimensions
# ----------------------------
set.seed(1)

grid_n <- 20
P <- grid_n^2
S <- 3
n_fleet <- 2
A_vec <- c(10, 42, 18)

patch_df <- expand.grid(x = seq_len(grid_n), y = seq_len(grid_n))
patch_df$p <- seq_len(P)

fishable_int <- rep(1L, P)  # all open

# helper: smooth spatial field (scaled to mean 1)
make_spatial_field <- function(grid_n, mode = c("hotspot", "gradient")) {
  mode <- match.arg(mode)
  xx <- rep(seq_len(grid_n), times = grid_n)
  yy <- rep(seq_len(grid_n), each  = grid_n)

  if (mode == "hotspot") {
    cx <- (grid_n + 1) / 2
    cy <- (grid_n + 1) / 2
    r2 <- (xx - cx)^2 + (yy - cy)^2
    field <- exp(-r2 / (2 * (0.22 * grid_n)^2))
  } else {
    field <- (xx + yy) / (2 * grid_n)
  }

  field / mean(field)
}

# helper: selectivity curve
make_sel <- function(n_age) {
  peak <- max(3, min(n_age - 2, round(0.6 * n_age)))
  a <- seq_len(n_age)
  sel <- exp(-0.5 * ((a - peak) / (0.2 * n_age))^2)
  sel / max(sel)
}

# fauna
fauna <- lapply(seq_len(S), function(s) {
  FaunaSpecies$new(m_at_age = rep(0.15, A_vec[s]))
})

# fleet builder (target fleet has spatial q if structured)
make_fleets <- function(mode = c("structured", "flat", "nearflat")) {
  mode <- match.arg(mode)

  spatial_signal <- switch(
    mode,
    structured = make_spatial_field(grid_n, "hotspot"),
    flat       = rep(1, P),
    nearflat   = rep(1, P) # keep q flat; nearflat perturbation will come from biomass
  )

  fleets <- vector("list", n_fleet)
  for (fleet_idx in seq_len(n_fleet)) {
    metiers <- vector("list", S)

    for (s in seq_len(S)) {
      n_age <- A_vec[s]
      q_patch <- if (fleet_idx == 1L) (0.10 * spatial_signal) else rep(0.08, P)

      metiers[[s]] <- list(
        sel_at_age = make_sel(n_age),
        price = c(2.0, 1.2, 0.8)[s],
        spatial_catchability = q_patch
      )
    }

    fleets[[fleet_idx]] <- list(
      cost_per_patch = rep(1.0, P),    # constant travel in this demo
      cost_per_unit_effort = 0.06,
      effort_cost_exponent = 1.4,
      metiers = metiers
    )
  }

  fleets
}

# storage builder with three modes:
# - structured: clear gradient field
# - flat: constant across patches + *extremely tiny* noise
# - nearflat: constant across patches + small lognormal noise (tunable)
make_storage <- function(mode = c("structured", "flat", "nearflat"),
                         nearflat_sdlog = 1e-4) {
  mode <- match.arg(mode)

  biom_field <- switch(
    mode,
    structured = make_spatial_field(grid_n, "gradient"),
    flat       = rep(1, P),
    nearflat   = rep(1, P)
  )

  storage <- vector("list", S)
  for (s in seq_len(S)) {
    n_age <- A_vec[s]
    age <- seq_len(n_age)

    age_shape <- exp(-0.05 * (age - 1))
    age_shape <- age_shape / mean(age_shape)

    b_p_a <- (biom_field %o% age_shape) * 100

    if (mode == "flat") {
      # essentially identical across patches; tiny noise to avoid exact ties
      b_p_a <- b_p_a * exp(matrix(rnorm(P * n_age, sd = 1e-10), nrow = P))
    }

    if (mode == "nearflat") {
      # small-but-real spatial noise -> near-flat objective surface
      # increasing nearflat_sdlog makes year-to-year "haphazard" effects more likely
      b_p_a <- b_p_a * exp(matrix(rnorm(P * n_age, sd = nearflat_sdlog), nrow = P))
    }

    storage[[s]] <- list(b_p_a = b_p_a)
  }

  storage
}

# ----------------------------
# scenario runner
# ----------------------------
run_scenario <- function(mode = c("structured", "flat", "nearflat"),
                         include_costs = FALSE,
                         flat_tol_sd = 1e-6,
                         flat_tol_abs = 1e-10,
                         nearflat_sdlog = 1e-4,
                         seed = 1) {

  mode <- match.arg(mode)
  set.seed(seed)

  storage <- make_storage(mode = mode, nearflat_sdlog = nearflat_sdlog)
  fleets  <- make_fleets(mode = mode)

  target_fleet <- 1L
  E_exo <- NULL

  pre <- precompute_baranov_inputs(
    storage = storage,
    fauna = fauna,
    fleets = fleets,
    target_fleet = target_fleet,
    E_exo = E_exo,
    P = P
  )

  alpha_mats <- pre$alpha_mats
  other_mort_mats <- pre$other_mort_mats
  biomass_mats <- pre$biomass_mats
  price_s <- pre$price_s

  cost_patch <- fleets[[target_fleet]]$cost_per_patch
  c0 <- fleets[[target_fleet]]$cost_per_unit_effort
  gamma <- fleets[[target_fleet]]$effort_cost_exponent

  Etot_target <- 600
  time_step <- 1

  out <- cpp_allocate_ifd_kkt_fullsolve_fast(
    Etot_target = Etot_target,
    alpha_mats = alpha_mats,
    other_mort_mats = other_mort_mats,
    biomass_mats = biomass_mats,
    price_s = price_s,
    cost_patch = cost_patch,
    c0 = c0,
    gamma = gamma,
    fishable_int = fishable_int,
    time_step = time_step,
    include_costs = include_costs,
    n_outer = 60,
    n_inner = 30,
    active_tol = 1e-14,
    flat_tol_sd = flat_tol_sd,
    flat_tol_abs = flat_tol_abs
  )

  list(out = out, storage = storage, fleets = fleets)
}

# ----------------------------
# plotting helpers
# ----------------------------
plot_maps <- function(out, title_prefix) {
  df <- patch_df
  df$E <- out$E_target
  df$rev_p <- out$revenue_p
  df$obj_p <- out$profit_p

  pE <- ggplot(df, aes(x = x, y = y, fill = E)) +
    geom_raster() +
    coord_equal() +
    labs(title = paste0(title_prefix, ": Effort allocation"), fill = "E") +
    theme_minimal()

  pR <- ggplot(df, aes(x = x, y = y, fill = rev_p)) +
    geom_raster() +
    coord_equal() +
    labs(title = paste0(title_prefix, ": Revenue surface at solution"), fill = "rev_p") +
    theme_minimal()

  pO <- ggplot(df, aes(x = x, y = y, fill = obj_p)) +
    geom_raster() +
    coord_equal() +
    labs(title = paste0(title_prefix, ": Accounting profit surface (rev - cost)"), fill = "profit_p") +
    theme_minimal()

  pH <- ggplot(df, aes(x = E)) +
    geom_histogram(bins = 50) +
    labs(title = paste0(title_prefix, ": Histogram of E"), x = "E", y = "Count") +
    theme_minimal()

  print(pE); print(pR); print(pO); print(pH)

  cat("\n", title_prefix, "summary(E):\n", sep = "")
  print(summary(out$E_target))
  cat("sd(E):", sd(out$E_target), "\n")
  cat("sum(E):", sum(out$E_target), "\n")
  if (!is.null(out$diagnostic)) print(out$diagnostic)
}

# A “haphazardness” metric across multiple runs of nearflat:
# shows how unstable the ranking/allocation is as you change the seed (year-to-year proxy)
nearflat_instability <- function(n_rep = 10, nearflat_sdlog = 1e-4, ...) {
  E_mat <- matrix(NA_real_, nrow = P, ncol = n_rep)
  for (i in seq_len(n_rep)) {
    rr <- run_scenario(mode = "nearflat", nearflat_sdlog = nearflat_sdlog, seed = 100 + i, ...)
    E_mat[, i] <- rr$out$E_target
  }
  # correlation of allocations between consecutive reps
  cors <- sapply(2:n_rep, function(i) cor(E_mat[, i - 1], E_mat[, i]))
  list(E_mat = E_mat, cor_consecutive = cors)
}

# ----------------------------
# Run all three scenarios
# ----------------------------
include_costs <- TRUE
flat_tol_sd <- 1e-6
flat_tol_abs <- 1e-10

res_struct <- run_scenario("structured", include_costs = include_costs,
                           flat_tol_sd = flat_tol_sd, flat_tol_abs = flat_tol_abs, seed = 1)

res_flat <- run_scenario("flat", include_costs = include_costs,
                         flat_tol_sd = flat_tol_sd, flat_tol_abs = flat_tol_abs, seed = 1)

# Near-flat: tune sdlog here:
nearflat_sdlog <- 1e-4
res_nearflat <- run_scenario("nearflat", include_costs = include_costs,
                             flat_tol_sd = flat_tol_sd, flat_tol_abs = flat_tol_abs,
                             nearflat_sdlog = nearflat_sdlog, seed = 1)

cat("Structured sum(E):", sum(res_struct$out$E_target), "\n")
cat("Flat sum(E):", sum(res_flat$out$E_target), "\n")
cat("Near-flat sum(E):", sum(res_nearflat$out$E_target), "\n")

# ----------------------------
# Plots
# ----------------------------
plot_maps(res_struct$out, "Scenario 1 (structured surface)")
plot_maps(res_flat$out,   "Scenario 2 (flat surface)")
plot_maps(res_nearflat$out, paste0("Scenario 3 (near-flat, sdlog=", nearflat_sdlog, ")"))

# ----------------------------
# Optional: instability check across repeated near-flat draws
# ----------------------------
inst <- nearflat_instability(n_rep = 8,
                             nearflat_sdlog = nearflat_sdlog,
                             include_costs = include_costs,
                             flat_tol_sd = flat_tol_sd,
                             flat_tol_abs = flat_tol_abs)

cat("\nNear-flat consecutive correlation of E across reps:\n")
print(inst$cor_consecutive)

inst_df <- data.frame(rep = 2:8, cor = inst$cor_consecutive)

pCor <- ggplot(inst_df, aes(x = rep, y = cor)) +
  geom_point() +
  geom_line() +
  labs(title = paste0("Near-flat instability: cor(E_t, E_{t-1})   sdlog=", nearflat_sdlog),
       x = "rep index", y = "correlation") +
  theme_minimal()

print(pCor)

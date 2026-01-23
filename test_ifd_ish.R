# ============================================================
# Worked example: ragged-age softmax/swap allocator + precompute
# + DIAGNOSTIC PLOTS (IFD-approx checks)
#
# - 20x20 patches (P = 400), numbered by expand_grid(x=1:20,y=1:20)
# - 2 species with DIFFERENT numbers of ages: A1=3, A2=5
# - 2 fleets: target fleet allocates effort; other fleet provides exogenous F
#
# Diagnostics added:
#   1) m vs effort (core-used highlighted)
#   2) u_swap vs effort (core-used highlighted)
#   3) summary metrics bars (smaller is better)
#   4) spatial maps of effort / u_swap / m / core-used
#
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

# ----------------------------
# 0) Grid helpers (20 x 20)
# ----------------------------

patch_to_xy <- function(patch) {
  x <- ((patch - 1) %% 20) + 1
  y <- ((patch - 1) %/% 20) + 1
  tibble(patch = patch, x = x, y = y)
}

plot_grid_field <- function(df, value_col, title, flip_y = TRUE) {
  df2 <- df |>
    left_join(patch_to_xy(df$patch), by = "patch") |>
    mutate(y_plot = if (flip_y) 20 - y + 1 else y)

  ggplot(df2, aes(x = x, y = y_plot, fill = .data[[value_col]])) +
    geom_tile() +
    coord_equal(expand = FALSE) +
    labs(title = title, x = "x", y = "y", fill = value_col) +
    theme_minimal()
}

# ----------------------------
# 1) Precompute (ragged by species)
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
# 2) Allocator (ragged)
# ----------------------------

patch_revenue_ragged <- function(e_p, v_mats, other_mats, b_mats, price_s, dt, p) {
  if (e_p <= 0) return(0)

  S <- length(v_mats)
  rev <- 0

  for (s in seq_len(S)) {
    v_a <- v_mats[[s]][p, ]
    other_a <- other_mats[[s]][p, ]
    b_a <- b_mats[[s]][p, ]

    F_a <- v_a * e_p
    Z_a <- F_a + other_a

    catch_frac_a <- ifelse(Z_a > 0, (F_a / Z_a) * (1 - exp(-dt * Z_a)), 0)
    rev <- rev + price_s[s] * sum(catch_frac_a * b_a)
  }

  rev
}

patch_obj_ragged <- function(e_p, v_mats, other_mats, b_mats, price_s, dt, p, c0, gamma, travel_p) {
  rev <- patch_revenue_ragged(e_p, v_mats, other_mats, b_mats, price_s, dt, p)
  cost <- c0 * (e_p^gamma + travel_p[p] * e_p)
  rev - cost
}

marginal_unconstrained_ragged <- function(
    e, delta,
    v_mats, other_mats, b_mats, price_s, dt,
    c0, gamma, travel_p,
    open_p
) {
  P <- length(e)
  m <- rep(NA_real_, P)

  for (p in seq_len(P)) {
    if (!open_p[p]) next

    base <- patch_obj_ragged(e[p], v_mats, other_mats, b_mats, price_s, dt, p, c0, gamma, travel_p)
    bump <- patch_obj_ragged(e[p] + delta, v_mats, other_mats, b_mats, price_s, dt, p, c0, gamma, travel_p)

    m[p] <- (bump - base) / delta
  }

  m
}

swap_marginal_from_m <- function(e, m, open_p) {
  P <- length(e)
  u_swap <- rep(NA_real_, P)

  for (p in seq_len(P)) {
    if (!open_p[p]) next

    donors <- which(open_p & seq_len(P) != p)
    pool <- sum(e[donors])
    if (pool <= 0) next

    w_donor <- e[donors] / pool
    donor_mean_m <- sum(w_donor * m[donors], na.rm = TRUE)

    u_swap[p] <- m[p] - donor_mean_m
  }

  u_swap
}

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

allocate_until_stable_patchwise <- function(
    e_init, E_tot,
    v_mats, other_mats, b_mats, price_s, dt,
    c0, gamma, travel_p, open_p,
    beta = 1, rho = 0.1,
    delta = NULL,
    max_iter = 40, tol = 1e-5,
    norm = "iqr",
    cap_frac = 0,
    record = TRUE
) {
  P <- length(e_init)
  open_idx <- which(open_p)

  e <- pmax(e_init, 0)
  e[!open_p] <- 0
  if (sum(e) > 0) {
    e <- e * (E_tot / sum(e))
  } else {
    e[open_idx] <- E_tot / length(open_idx)
  }

  if (is.null(delta)) {
    delta <- 0.02 * (E_tot / length(open_idx))
  }

  history <- if (record) vector("list", max_iter) else NULL

  for (k in seq_len(max_iter)) {
    e_old <- e

    m <- marginal_unconstrained_ragged(
      e = e, delta = delta,
      v_mats = v_mats, other_mats = other_mats, b_mats = b_mats,
      price_s = price_s, dt = dt,
      c0 = c0, gamma = gamma, travel_p = travel_p,
      open_p = open_p
    )

    u_swap <- swap_marginal_from_m(e, m, open_p)

    upd <- update_effort_softmax(
      e = e, E_tot = E_tot,
      u_swap = u_swap, open_p = open_p,
      beta = beta, rho = rho,
      norm = norm, cap_frac = cap_frac
    )
    e <- upd$e

    rel_change <- sum(abs(e - e_old)) / E_tot

    if (record) {
      history[[k]] <- tibble(
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

  list(e = e, history = history, delta = delta)
}

# ----------------------------
# 3) Diagnostics (IFD approximation checks)
# ----------------------------

ifd_diagnostic <- function(df_last, effort_floor = 0.1, top_frac = 0.9) {
  # "core-used": patches with effort > effort_floor and within top_frac of effort mass
  df <- df_last |>
    dplyr::filter(open) |>
    dplyr::arrange(dplyr::desc(effort)) |>
    dplyr::mutate(
      cum_share = cumsum(effort) / sum(effort),
      used = effort > effort_floor,
      core_used = used & cum_share <= top_frac
    )

  core_m <- df$m[df$core_used]
  core_u <- df$u_swap[df$core_used]

  metrics <- tibble(
    metric = c("IQR(m) on core-used",
               "SD(m) on core-used",
               "RMS(u_swap) on core-used",
               "MAD(u_swap) on core-used"),
    value  = c(stats::IQR(core_m, na.rm = TRUE),
               stats::sd(core_m, na.rm = TRUE),
               sqrt(mean(core_u^2, na.rm = TRUE)),
               stats::mad(core_u, constant = 1, na.rm = TRUE))
  )

  p_m <- ggplot(df, aes(x = effort, y = m)) +
    geom_point(aes(shape = core_used), size = 2.6, alpha = 0.9) +
    labs(
      title = "IFD approximation: marginal equalization check",
      subtitle = "Better when m is nearly constant across core-used patches",
      x = "Effort",
      y = "Unconstrained marginal m = d obj_p / d e_p",
      shape = "Core-used"
    ) +
    theme_minimal()

  p_u <- ggplot(df, aes(x = effort, y = u_swap)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(aes(shape = core_used), size = 2.6, alpha = 0.9) +
    labs(
      title = "IFD approximation: swap incentive residuals",
      subtitle = "Better when u_swap ~ 0 across core-used patches",
      x = "Effort",
      y = "Swap-based marginal u_swap",
      shape = "Core-used"
    ) +
    theme_minimal()

  p_metrics <- ggplot(metrics, aes(x = metric, y = value)) +
    geom_col() +
    labs(title = "IFD approximation summary metrics (smaller is better)", x = NULL, y = "Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  list(df = df, metrics = metrics, p_m = p_m, p_u = p_u, p_metrics = p_metrics)
}

# ============================================================
# 4) WORKED EXAMPLE DATA GENERATION (ragged ages)
# ============================================================

set.seed(1)

grid <- tidyr::expand_grid(x = 1:20, y = 1:20) |>
  mutate(patch = row_number())

P <- nrow(grid)

# Two species with different age counts
A1 <- 3
A2 <- 5
n_species <- 2
n_fleet <- 2

gauss_bump <- function(x, y, x0, y0, sigma) {
  exp(-((x - x0)^2 + (y - y0)^2) / (2 * sigma^2))
}

make_smooth_field <- function(bumps, base = 0.1, noise_sd = 0.0) {
  z <- rep(base, P)
  for (i in seq_len(nrow(bumps))) {
    z <- z + bumps$amp[i] * gauss_bump(grid$x, grid$y, bumps$x0[i], bumps$y0[i], bumps$sigma[i])
  }
  if (noise_sd > 0) z <- z + rnorm(P, 0, noise_sd)
  pmax(z, 0)
}

# Biomass surfaces
bio1 <- make_smooth_field(
  tibble(x0 = c(6, 15), y0 = c(7, 14), amp = c(1.2, 0.7), sigma = c(3.5, 5.0)),
  base = 0.15, noise_sd = 0.02
)
bio2 <- make_smooth_field(
  tibble(x0 = c(14, 5), y0 = c(6, 16), amp = c(1.0, 0.8), sigma = c(4.5, 4.0)),
  base = 0.12, noise_sd = 0.02
)

age_mult1 <- c(0.7, 1.0, 1.3)
age_mult2 <- c(0.5, 0.8, 1.0, 1.2, 1.4)

storage <- vector("list", n_species)
storage[[1]] <- list(b_p_a = 300 * bio1 %o% age_mult1) # [P x A1]
storage[[2]] <- list(b_p_a = 260 * bio2 %o% age_mult2) # [P x A2]

fauna <- vector("list", n_species)
fauna[[1]] <- list(m_at_age = c(0.15, 0.10, 0.08))
fauna[[2]] <- list(m_at_age = c(0.18, 0.14, 0.11, 0.09, 0.08))

make_vul <- function(base_field, age_mult, scale = 0.03) {
  # returns [P x A_s]
  scale * (base_field %o% (age_mult / max(age_mult)))
}

vul_base_f1 <- make_smooth_field(
  tibble(x0 = c(8, 16), y0 = c(14, 7), amp = c(1.0, 0.8), sigma = c(6, 5)),
  base = 0.8, noise_sd = 0.01
)
vul_base_f2 <- make_smooth_field(
  tibble(x0 = c(12, 4), y0 = c(5, 18), amp = c(1.0, 0.6), sigma = c(6, 6)),
  base = 0.7, noise_sd = 0.01
)

fleets <- vector("list", n_fleet)

fleets[[1]] <- list(
  metiers = list(
    list(price = 6,  vul_p_a = make_vul(vul_base_f1, age_mult1, scale = 0.030)),
    list(price = 9,  vul_p_a = make_vul(vul_base_f1, age_mult2, scale = 0.026))
  )
)

fleets[[2]] <- list(
  metiers = list(
    list(price = 6,  vul_p_a = make_vul(vul_base_f2, age_mult1, scale = 0.028)),
    list(price = 9,  vul_p_a = make_vul(vul_base_f2, age_mult2, scale = 0.024))
  )
)

# Exogenous effort (fleet 2 has fixed spatial effort)
E_exo <- vector("list", n_fleet)
E_exo[[1]] <- rep(0, P)
E_exo[[2]] <- {
  w <- 0.3 + gauss_bump(grid$x, grid$y, 16, 16, 6)
  200 * w / sum(w) * P
}

travel_p <- with(grid, sqrt((x - 10)^2 + (y - 10)^2))
travel_p <- travel_p / max(travel_p)

open_p <- rep(TRUE, P)
# open_p[with(grid, x >= 18 & y >= 18)] <- FALSE  # optional closure block

c0 <- 2
gamma <- 1.2
dt <- 1

# ============================================================
# 5) Run precompute + allocator
# ============================================================

target_fleet <- 1

pre <- precompute_baranov_inputs_softmax(
  storage = storage,
  fauna = fauna,
  fleets = fleets,
  target_fleet = target_fleet,
  E_exo = E_exo,
  P = P
)

E_tot <- 1000
e0 <- rep(0, P)
e0[open_p] <- E_tot / sum(open_p)

res <- allocate_until_stable_patchwise(
  e_init = e0, E_tot = E_tot,
  v_mats = pre$v_mats,
  other_mats = pre$other_mats,
  b_mats = pre$b_mats,
  price_s = pre$price_s,
  dt = dt,
  c0 = c0, gamma = gamma,
  travel_p = travel_p,
  open_p = open_p,
  beta = 7, rho = 0.1,
  max_iter = 15, tol = 1e-5,
  cap_frac = 0.03,
  record = TRUE
)

hist_df <- dplyr::bind_rows(res$history)
df_last <- hist_df |> dplyr::filter(iter == max(iter))

# ============================================================
# 6) Plots: convergence + IFD diagnostics + spatial maps
# ============================================================

# Convergence (effort trajectories)
ggplot(hist_df, aes(x = iter, y = effort, group = patch)) +
  geom_line(alpha = 0.15) +
  labs(title = "Convergence of patch effort (inner loop)", x = "Iteration", y = "Effort") +
  theme_minimal()

# ---- IFD approximation diagnostics ----
diag <- ifd_diagnostic(df_last, effort_floor = 0.1, top_frac = 0.9)
diag$p_m
diag$p_u
diag$p_metrics
diag$metrics

# ---- Spatial maps: effort, u_swap, m ----
plot_grid_field(df_last |> select(patch, effort), "effort", "Final effort allocation (target fleet)")
plot_grid_field(df_last |> select(patch, u_swap), "u_swap", "Final swap-based marginal u_swap (target fleet)")
plot_grid_field(df_last |> select(patch, m), "m", "Final unconstrained marginal m (target fleet)")

# ---- Spatial map: core-used patches (binary) ----
df_core <- diag$df |> select(patch, core_used)
plot_grid_field(df_core, "core_used", "Core-used patches (binary)")

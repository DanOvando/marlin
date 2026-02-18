#' Allocate fishing effort across patches
#'
#' @description
#' Updates patch-level effort for one or more fleets based on an objective signal
#' contained in the \code{buffet} output from \code{\link{go_fish}}.
#'
#' @details
#' This helper is called by \code{\link{simmar}} each time step. For each fleet,
#' it converts the selected objective (e.g. revenue, catch, profit, or a per-unit
#' variant) into a smooth multiplicative update so that higher-objective patches
#' receive more effort, subject to \code{open_patch}.
#'
#' @param effort_by_patch Numeric matrix with patches in rows and fleets in columns.
#' @param total_effort_by_fleet Optional numeric vector of total effort per fleet. If
#'   \code{NULL}, totals are computed from \code{effort_by_patch}.
#' @param buffet Named list from \code{\link{go_fish}} with matrix outputs (for example
#'   \code{r_p_fl}, \code{c_p_fl}, \code{prof_p_fl}, and per-unit variants).
#' @param fleets Named list of fleet objects (each with \code{$spatial_allocation}).
#' @param open_patch Logical matrix (patches x fleets) indicating where each fleet may fish.
#' @param clip_z Numeric scalar controlling how strongly extreme objective values affect the update.
#' @param scale Character; method used to scale the objective before updating effort.
#' @param eps_mix Numeric between 0 and 1; mix-in fraction for uniform exploration.
#' @param floor_frac Numeric; optional effort floor applied across open patches.
#' @param flatness_tol Numeric; threshold for treating the objective as flat across patches.
#' @param min_scale_abs Numeric; minimum scale to avoid division by near-zero values.
#' @param adaptive_floor_pct Numeric; fraction of the median objective used as an adaptive scale floor.
#'
#' @return
#' A named list containing the updated effort matrix (\code{effort_new}) and additional
#' diagnostic matrices/vectors used to track the update (for example standardised objectives
#' and flatness indicators).
#'
#' @seealso \code{\link{go_fish}}, \code{\link{calc_marginal_value}}, \code{\link{simmar}}

allocate_effort <- function(
    effort_by_patch,
    total_effort_by_fleet = NULL,
    buffet,
    fleets,
    open_patch,
    clip_z = 5,
    scale = c("mad", "iqr", "sd"),
    eps_mix = 0.0,
    floor_frac = 0.0,
    flatness_tol = 1e-3,
    min_scale_abs = 1e-10,
    adaptive_floor_pct = 0.01
) {
  scale <- match.arg(scale)

  if (is.null(total_effort_by_fleet)) {
    total_effort_by_fleet <- colSums(effort_by_patch)
  }
  # --- Coerce inputs to matrices -----------------------------------------------

  # incase total effort has changed, first normalize each fleets effort, then redistribute total effort
  effort_by_patch <- as.matrix(effort_by_patch)
  colnames(effort_by_patch) <- names(fleets)
  effort_by_patch <- (effort_by_patch / rep(colSums(effort_by_patch), each = nrow(effort_by_patch))) * rep(total_effort_by_fleet, each = nrow(effort_by_patch))

  n_patches <- nrow(effort_by_patch)
  n_fleets <- ncol(effort_by_patch)
  fleet_names <- (names(fleets))

  # open_patch: accept vector (shared) or matrix (per-fleet)
  if (is.null(dim(open_patch))) {
    open_patch <- matrix(open_patch, nrow = n_patches, ncol = n_fleets)
  } else {
    open_patch <- as.matrix(open_patch)
  }
  storage.mode(open_patch) <- "logical"

  # --- Map spatial_allocation to buffet matrices --------------------------------

  alloc_to_mat <- c(
    revenue          = "r_p_fl",
    catch            = "c_p_fl",
    profit           = "prof_p_fl",
    rpue             = "rpue_p_fl",
    cpue             = "cpue_p_fl",
    ppue             = "ppue_p_fl",
    marginal_revenue = "mr_p_fl",
    marginal_profit  = "mp_p_fl"
  )


  # --- Pre-allocate output -----------------------------------------------------
  effort_new <- matrix(0, nrow = n_patches, ncol = n_fleets)
  velocity_out <- matrix(0, nrow = n_patches, ncol = n_fleets)
  z_out <- matrix(0, nrow = n_patches, ncol = n_fleets)
  flat_out <- logical(n_fleets)
  cv_out <- numeric(n_fleets)
  range_out <- numeric(n_fleets)

  # --- Per-fleet allocation ----------------------------------------------------
  for (fl in seq_len(n_fleets)) {
    fl_name <- fleet_names[fl]
    # Look up which objective this fleet uses
    alloc_type <- fleets[[fl_name]]$spatial_allocation

    eta <- fleets[[fl_name]]$eta

    # --- Manual allocation: distribute effort by fishing_grounds weights --------
    if (alloc_type == "manual") {
      weights <- fleets[[fl_name]]$fishing_grounds$fishing_ground
      # Zero out weights in MPA-closed patches
      weights[!open_patch[, fl]] <- 0
      total_effort <- sum(effort_by_patch[, fl])

      weight_sum <- sum(weights)
      if (weight_sum <= 0) {
        warning(
          sprintf("Fleet '%s' has spatial_allocation = 'manual' but all fishing_grounds weights are zero (possibly due to MPA closures). Distributing effort uniformly across open patches.", fl_name),
          call. = FALSE
        )
        open_idx_manual <- which(open_patch[, fl])
        if (length(open_idx_manual) > 0) {
          effort_new[open_idx_manual, fl] <- total_effort / length(open_idx_manual)
        }
      } else {
        effort_new[, fl] <- total_effort * (weights / weight_sum)
      }

      # No velocity or z to report for manual fleets
      flat_out[fl] <- NA
      cv_out[fl] <- NA
      range_out[fl] <- NA
      next
    }

    # --- Uniform allocation: equal effort across all open patches ---------------
    if (alloc_type == "uniform") {
      total_effort <- sum(effort_by_patch[, fl])
      open_idx_uniform <- which(open_patch[, fl])

      if (length(open_idx_uniform) > 0) {
        effort_new[open_idx_uniform, fl] <- total_effort / length(open_idx_uniform)
      }

      # No velocity or z to report for uniform fleets
      flat_out[fl] <- NA
      cv_out[fl] <- NA
      range_out[fl] <- NA
      next
    }

    mat_name <- alloc_to_mat[alloc_type]

    if (is.na(mat_name)) {
      stop(sprintf(
        "Fleet '%s' has spatial_allocation = '%s', which is not one of: %s, manual, uniform",
        fl_name, alloc_type, paste(names(alloc_to_mat), collapse = ", ")
      ))
    }

    effort <- effort_by_patch[, fl]
    objective <- buffet[[mat_name]][, fl_name]
    open <- open_patch[, fl]

    open_idx <- which(open)
    total_effort <- sum(effort)

    # Step 1: FLATNESS DETECTION ------------------------------------------------
    obj_open <- objective[open_idx]
    obj_center <- stats::median(obj_open, na.rm = TRUE)
    obj_rng <- diff(range(obj_open, na.rm = TRUE))

    obj_scale <- switch(scale,
                        mad = stats::mad(obj_open, constant = 1, na.rm = TRUE),
                        iqr = stats::IQR(obj_open, na.rm = TRUE) / 1.349,
                        sd  = stats::sd(obj_open, na.rm = TRUE)
    )

    # Flatness detection uses the range, not obj_scale.
    # MAD/IQR can severely underestimate spread in bimodal distributions
    # (e.g. 80 patches at RPUE~1.44, 20 patches at RPUE~1.66 gives MAD-based
    # CV ~3e-5 even though the range-based CV is ~15%). The range captures
    # the full spread regardless of distribution shape.
    range_cv <- obj_rng / (abs(obj_center) + 1e-100)
    is_flat <- (range_cv < flatness_tol) || (obj_rng < min_scale_abs)

    cv_out[fl] <- range_cv
    range_out[fl] <- obj_rng

    if (is_flat) {
      flat_out[fl] <- TRUE
      # Distribute uniformly across open patches only; closed patches stay at 0
      effort_new[open_idx, fl] <- total_effort / length(open_idx)
      next
    }

    flat_out[fl] <- FALSE

    # Step 2: ADAPTIVE SCALE FLOOR ----------------------------------------------
    min_scale_adaptive <- pmax(adaptive_floor_pct * abs(obj_center), min_scale_abs)
    obj_scale <- pmax(obj_scale, min_scale_adaptive)

    # Step 3: STANDARDIZE -------------------------------------------------------
    z <- (objective - obj_center) / obj_scale
    z <- pmax(pmin(z, clip_z), -clip_z)
    z[!open] <- 0
    z_out[, fl] <- z

    # Step 4: VELOCITY FIELD ----------------------------------------------------
    z_open_mean <- mean(z[open_idx])
    v <- z - z_open_mean
    v[!open] <- 0
    velocity_out[, fl] <- v

    # Step 5: MULTIPLICATIVE UPDATE ---------------------------------------------
    e_open <- effort
    e_open[!open] <- 0

    if (floor_frac > 0) {
      floor_effort <- floor_frac * (total_effort / length(open_idx))
      e_open[open] <- pmax(e_open[open], floor_effort)
    }

    log_e <- rep(-Inf, n_patches)
    log_e[open] <- log(pmax(e_open[open], 1e-300))

    log_e_prop <- log_e + eta * v

    mx <- max(log_e_prop[open])
    w <- numeric(n_patches)
    w[open] <- exp(log_e_prop[open] - mx)
    w[!open] <- 0

    e_new <- w * (total_effort / sum(w))

    # Step 6: EXPLORATION MIX ---------------------------------------------------
    if (eps_mix > 0) {
      baseline <- numeric(n_patches)
      baseline[open] <- total_effort / length(open_idx)
      e_new <- (1 - eps_mix) * e_new + eps_mix * baseline
    }

    effort_new[, fl] <- e_new
  }

  # --- Label outputs -----------------------------------------------------------
  colnames(effort_new) <- fleet_names
  colnames(velocity_out) <- fleet_names
  colnames(z_out) <- fleet_names
  names(flat_out) <- fleet_names
  names(cv_out) <- fleet_names
  names(range_out) <- fleet_names

  list(
    effort_new = effort_new,
    velocity_by_patch = velocity_out,
    z_by_patch = z_out,
    flat_objective = flat_out,
    cv = cv_out,
    obj_range = range_out
  )
}

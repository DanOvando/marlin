#' Allocate Fishing Effort Using Multiplicative Velocity-Field Update
#'
#' @description
#' Updates patch-level effort for one or more fleets using a multiplicative
#' (replicator-style) gradient-flow heuristic derived from an objective signal.
#' Each fleet's objective is selected automatically from the \code{buffet} based
#' on its \code{spatial_allocation} setting (e.g., "rpue", "profit", "catch").
#' Fleets with \code{spatial_allocation = "manual"} bypass the optimizer
#' entirely: their effort is distributed proportionally to the continuous
#' weights in \code{fishing_grounds$fishing_ground}.
#' Fleets with \code{spatial_allocation = "uniform"} distribute effort
#' equally across all open patches, ignoring both the objective signal and
#' fishing grounds weights.
#'
#' The allocator interprets the objective as a **velocity field** rather than a
#' target distribution. Effort changes proportionally according to:
#'
#' \deqn{e_{new} \propto e_{old} \exp(\eta v)}
#'
#' where \eqn{v} is a centered, scaled, and clipped version of the objective.
#'
#' This is a **single-step response**: fleets observe the current objective
#' landscape and make one reallocation decision each. Convergence to equilibrium
#' emerges from the bioeconomic feedback loop in \code{\link{simmar}}, where
#' effort decisions play out through population dynamics and generate updated
#' objective signals at each timestep.
#'
#' This approach guarantees:
#' \itemize{
#'   \item Non-negative effort
#'   \item Conservation of total effort (per fleet)
#'   \item No change when incentives are equal across open patches
#'   \item High numerical stability via log-space updates
#'   \item Stable behavior when objective is flat (no spurious edge concentration)
#' }
#'
#' @param effort_by_patch Numeric matrix of effort by patch (rows) and fleet
#'   (columns). Column names should be fleet names.
#' @param buffet List returned by \code{\link{go_fish}} with
#'   \code{output_format = "matrix"}. Must contain: \code{r_p_fl}, \code{c_p_fl},
#'   \code{prof_p_fl}, \code{rpue_p_fl}, \code{cpue_p_fl}, \code{ppue_p_fl}.
#'   Each fleet picks its objective from this buffet based on its
#'   \code{spatial_allocation} setting.
#' @param fleets List of fleet objects. Each fleet must have a
#'   \code{spatial_allocation} element set to one of: \code{"revenue"},
#'   \code{"catch"}, \code{"profit"}, \code{"rpue"}, \code{"cpue"},
#'   \code{"ppue"}, \code{"manual"}, or \code{"uniform"}. When set to \code{"manual"}, the
#'   fleet's \code{fishing_grounds$fishing_ground} values (continuous,
#'   between 0 and 1) define the relative spatial distribution of effort.
#'   Total effort is preserved but distributed proportionally to these
#'   weights rather than optimized against an objective.
#'   When set to \code{"uniform"}, total effort is distributed equally
#'   across all open patches regardless of any objective or weights.
#' @param open_patch Logical matrix indicating which patches are open to fishing,
#'   by patch (rows) and fleet (columns). Same dimensions as
#'   \code{effort_by_patch}. Allows fleet-specific closures (e.g., different
#'   fishing grounds per fleet). Can also be a logical vector if all fleets share
#'   the same open patches.
#' @param eta Numeric responsiveness parameter controlling magnitude of effort
#'   movement (default: 0.2). Higher values = faster adjustment to incentives.
#' @param clip_z Numeric maximum absolute standardized objective value allowed
#'   (default: 5). Prevents extreme updates from outliers.
#' @param scale Character string specifying method to scale the objective:
#'   "mad" (median absolute deviation, default), "iqr" (interquartile range),
#'   or "sd" (standard deviation). MAD is most robust to outliers.
#' @param eps_mix Numeric optional exploration mixing fraction between 0 and 1
#'   (default: 0). When > 0, mixes optimal allocation with uniform distribution
#'   to allow re-entry into abandoned patches.
#' @param floor_frac Numeric optional fractional effort floor relative to mean
#'   open-patch effort (default: 0). Prevents complete abandonment of patches.
#' @param flatness_tol Numeric threshold for detecting flat objectives using
#'   range-based coefficient of variation (range / median; default: 1e-3 = 0.1%).
#'   Range-based CV is used instead of MAD/IQR-based CV because the latter
#'   can severely underestimate spread in bimodal distributions.
#' @param min_scale_abs Numeric absolute minimum scale value to prevent
#'   division by tiny numbers (default: 1e-10).
#' @param adaptive_floor_pct Numeric percentage of median for adaptive scale floor
#'   (default: 0.01 = 1%).
#'
#' @return A list with:
#' \describe{
#'   \item{effort_new}{Numeric matrix of updated effort (patches x fleets),
#'     with column names preserved from \code{effort_by_patch}}
#'   \item{velocity_by_patch}{Numeric matrix of velocity fields (patches x fleets)}
#'   \item{z_by_patch}{Numeric matrix of standardized objective values (patches x fleets)}
#'   \item{flat_objective}{Named logical vector (one per fleet): TRUE if objective was
#'     detected as flat for that fleet}
#'   \item{cv}{Named numeric vector (one per fleet): range-based coefficient of
#'     variation of objective (used for flatness detection)}
#'   \item{obj_range}{Named numeric vector (one per fleet): range of objective across open patches}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run go_fish in matrix mode, then allocate
#' buffet <- go_fish(e_p_fl, fauna, n_p_a, fleets, output_format = "matrix")
#'
#' result <- allocate_effort(
#'   effort_by_patch = e_p_fl,
#'   buffet = buffet,
#'   fleets = fleets,
#'   open_patch = fleet_fishable
#' )
#'
#' # Updated effort matrix, ready for next simmar timestep
#' e_p_fl_new <- result$effort_new
#'
#' # Check which fleets saw flat objectives
#' result$flat_objective
#' }
#'
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

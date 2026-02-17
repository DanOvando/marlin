#' Deep-clone all metiers in a fleet list
#'
#' Creates independent copies of all metier R6 objects so that
#' modifications during tuning don't mutate the caller's fleet objects.
#'
#' @param fleets a fleet object (list of fleets)
#'
#' @return cloned fleet list
#' @keywords internal
clone_fleet_metiers <- function(fleets) {
  tfleets <- fleets
  for (i in seq_along(tfleets)) {
    for (j in seq_along(tfleets[[i]]$metiers)) {
      tfleets[[i]]$metiers[[j]] <- fleets[[i]]$metiers[[j]]$clone(deep = TRUE)
    }
  }
  tfleets
}


#' Link function mapping raw catchability to (0, 1)
#'
#' Applies the transformation \code{q = raw_q / (1 + raw_q)}, which is
#' equivalent to \code{plogis(log(raw_q))}. This smoothly maps positive
#' real values to (0, 1), approximately preserving values much less than 1
#' (when \code{raw_q << 1}, \code{q ≈ raw_q}) while preventing q from
#' ever reaching or exceeding 1. The smooth gradient avoids the flat
#' boundary regions that hard clamping would create for the optimizer.
#'
#' @param raw_q positive numeric, the unconstrained catchability value
#'
#' @return numeric in (0, 1)
#' @keywords internal
q_link <- function(raw_q) {
  raw_q / (1 + raw_q)
}


#' Inverse link function mapping (0, 1) back to positive reals
#'
#' Applies the transformation \code{raw_q = q / (1 - q)}, the inverse of
#' \code{\link{q_link}}.
#'
#' @param q numeric in (0, 1)
#'
#' @return positive numeric
#' @keywords internal
q_link_inv <- function(q) {

  q / (1 - q)
}


#' Set catchability for a single fleet-species metier
#'
#' Assigns scalar catchability, rescales the spatial catchability vector
#' to have the correct mean, and rebuilds \code{vul_p_a}. Catchability
#' is enforced to lie in (0, 1).
#'
#' @param metier a metier R6 object (modified in place)
#' @param q_value the target scalar catchability. Must be in [0, 1] when
#'   \code{use_link = FALSE}. When \code{use_link = TRUE}, can be any
#'   non-negative value and will be mapped to (0, 1) via \code{\link{q_link}}.
#' @param use_link logical; if TRUE, applies \code{\link{q_link}} to map
#'   \code{q_value} from the positive reals into (0, 1). Use this during
#'   numerical optimization to ensure smooth gradients. If FALSE (default),
#'   \code{q_value} is used directly and must already be in [0, 1].
#'
#' @return the (modified) metier, invisibly
#' @keywords internal
set_metier_catchability <- function(metier, q_value, use_link = FALSE) {

  if (!is.finite(q_value) || q_value < 0) {
    q_value <- 0
  }

  if (use_link) {
    q <- q_link(q_value)
  } else {
    if (q_value > 1) {
      stop(
        sprintf(
          "Catchability = %.4g is > 1. Increase base_effort or reduce target exploitation.",
          q_value
        ),
        call. = FALSE
      )
    }
    q <- q_value
  }

  metier$catchability <- q

  # If spatial_catchability is all zeros (e.g. from a prior zeroed-out metier),
  # reset to a uniform vector so it can be rescaled
  if (all(metier$spatial_catchability == 0)) {
    metier$spatial_catchability <- rep(1, length(metier$spatial_catchability))
  }

  mean_sq <- mean(metier$spatial_catchability)
  mean_sq <- ifelse(mean_sq == 0, 1e-9, mean_sq)

  metier$spatial_catchability <- (metier$spatial_catchability / mean_sq) * q

  metier$vul_p_a <- outer(metier$spatial_catchability, metier$sel_at_age, `*`)

  invisible(metier)
}


#' Compute effective effort per fleet
#'
#' For each fleet, calculates the share of base effort that overlaps with
#' viable habitat (\code{b0_p > 0}) within the fleet's fishing grounds.
#' Effort is a fleet-level quantity (not species-specific): a fleet exerts
#' one pool of effort that affects each species differently through
#' catchability. Any single species' habitat can be used for the overlap
#' calculation.
#'
#' @param fauna_species a single fauna element (used only for \code{b0_p})
#' @param fleets a fleet object (list of fleets)
#'
#' @return named numeric vector of effective effort per fleet
#' @keywords internal
compute_effort_by_fleet <- function(fauna_species, fleets) {

  e_fl <- vapply(fleets, function(fl) {
    fg <- fl$fishing_grounds$fishing_ground
    n_active <- sum(fg > 0)
    if (n_active == 0) return(0)
    habitat_overlap <- mean(fauna_species$b0_p[fg > 0] > 0)
    (fl$base_effort / n_active) * habitat_overlap
  }, numeric(1))

  e_fl
}


#' Calibrate cost parameters from equilibrium simulation output
#'
#' Extracts revenue and effort from an equilibrium time step, then solves
#' for \code{cost_per_unit_effort} and \code{effort_reference} for each
#' fleet so that costs match the target cost-to-revenue ratio.
#'
#' The calibration formula is:
#' \deqn{c0 = \frac{\text{cr\_ratio} \times \text{revenue}}{E^{ref} \times P \times (1 + \theta)}}
#'
#' The effort cost exponent \eqn{\gamma} does not appear in the inversion
#' because tuning runs under constant effort, which distributes effort
#' roughly uniformly across open patches. With \eqn{E_l / E^{ref} \approx 1}
#' everywhere, \eqn{1^\gamma = 1} regardless of \eqn{\gamma}, so the
#' exponent drops out.
#'
#' @param eq a single equilibrium time step from \code{\link{simmar}} output
#'   (i.e. \code{sim[[length(sim)]]})
#' @param fleets a fleet object (list of fleets)
#'
#' @return a data.frame with columns: fleet, cost_per_unit_effort,
#'   effort_reference, revenue, total_effort, n_open_patches, cr_ratio,
#'   travel_weight, implied_cr
#' @keywords internal
calibrate_fleet_costs <- function(eq, fleets) {

  fleet_names <- names(fleets)

  # Revenue: sum across species, by fleet
  revenue <-
    purrr::map(
      eq,
      ~ tibble::rownames_to_column(
        data.frame(revenue = colSums(.x$r_p_fl, na.rm = TRUE)),
        "fleet"
      )
    ) |>
    purrr::list_rbind(names_to = "critter") |>
    dplyr::group_by(fleet) |>
    dplyr::summarise(revenue = sum(revenue), .groups = "drop")

  # Effort by patch and fleet (same for all critters; use first)
  effort <-
    purrr::map(
      eq[1],
      ~ data.frame(.x$e_p_fl) |> dplyr::mutate(patch = dplyr::row_number())
    ) |>
    purrr::list_rbind(names_to = "critter") |>
    tidyr::pivot_longer(
      -c(critter, patch),
      names_to = "fleet",
      values_to = "effort"
    )

  # Fleet-level parameters needed for calibration
  fleet_params <-
    purrr::map(
      fleets,
      ~ data.frame(
        cr_ratio = .x$cr_ratio,
        travel_weight = .x$travel_weight
      )
    ) |>
    purrr::list_rbind(names_to = "fleet")

  # Fishing grounds mask (to identify open patches)
  fishing_grounds <-
    purrr::map(
      fleets,
      ~ data.frame(
        patch = seq_along(.x$fishing_grounds$fishing_ground),
        fishing_ground = .x$fishing_grounds$fishing_ground
      )
    ) |>
    purrr::list_rbind(names_to = "fleet")

  # Effort statistics: mean effort per open patch (where fishing is allowed)
  effort_stats <- effort |>
    dplyr::left_join(fishing_grounds, by = c("fleet", "patch")) |>
    dplyr::filter(fishing_ground > 0) |>
    dplyr::group_by(fleet) |>
    dplyr::summarise(
      total_effort = sum(effort),
      n_open_patches = dplyr::n_distinct(patch),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      effort_reference = total_effort / n_open_patches
    )

  # Solve for c0 and validate
  cost_calibration <- effort_stats |>
    dplyr::left_join(revenue, by = "fleet") |>
    dplyr::left_join(fleet_params, by = "fleet") |>
    dplyr::mutate(
      cost_per_unit_effort = (cr_ratio * revenue) /
        (effort_reference * n_open_patches * (1 + travel_weight)),
      implied_cr = (cost_per_unit_effort * effort_reference *
                      n_open_patches * (1 + travel_weight)) / revenue
    )

  cost_calibration
}


#' Tune fleet catchability to achieve target initial conditions
#'
#' Adjusts the catchability coefficient of each fleet's metiers so that the
#' simulated fishery reaches either a target fishing mortality rate or a
#' target depletion level (B/B0). Optionally tunes cost parameters to match
#' a target cost-to-revenue ratio.
#'
#' Catchability is enforced to lie in (0, 1). In the \code{"f"} path,
#' infeasible values (q >= 1) produce an error. In the \code{"depletion"}
#' path, a logistic link function (\code{q = raw_q / (1 + raw_q)})
#' smoothly maps the optimizer's unconstrained catchability proposals
#' into (0, 1), preserving smooth gradients and avoiding hard boundary
#' issues during numerical solution.
#'
#' When \code{tune_type = "depletion"} and \code{tune_costs = TRUE}, cost
#' calibration uses a two-pass approach: (1) initial calibration using a
#' heuristic (F = 1 - depletion) before the depletion solver runs, then
#' (2) refined calibration using the actual equilibrium after depletion
#' tuning converges. This ensures costs accurately reflect the final
#' equilibrium conditions.
#'
#' Note that tuning is approximate: post-tuning values will not perfectly
#' match inputs since some tuning steps depend on prior steps.
#'
#' @param fauna a fauna object
#' @param fleets a fleet object
#' @param years the number of years to simulate during tuning
#' @param tune_type one of \code{"f"} (or \code{"explt"}) to tune to a target
#'   fishing mortality rate, or \code{"depletion"} to tune to a target B/B0
#' @param tune_costs logical; if TRUE, tune cost_per_unit_effort to match each
#'   fleet's target cost-to-revenue ratio
#'
#' @return tuned fleet object
#' @export
tune_fleets <- function(fauna,
                        fleets,
                        years = 50,
                        tune_type = "f",
                        tune_costs = TRUE,
                        depl_tol = 0.025) {

  # Treat "explt" as an alias for "f"
  if (tune_type == "explt") {
    tune_type <- "f"
  }

  tune_type <- match.arg(tune_type, c("f", "depletion"))

  tfleets <- clone_fleet_metiers(fleets)

  fleet_names <- names(tfleets)
  fauni <- names(fauna)

  # Store and override fleet models to constant effort for tuning
  og_fleet_model <- setNames(
    vapply(tfleets, function(fl) fl$fleet_model, character(1)),
    fleet_names
  )
  for (f in fleet_names) {
    tfleets[[f]]$fleet_model <- "constant_effort"
  }

  # --- Normalize p_explt to sum to 1 across fleets for each species ---
  for (s in fauni) {
    p_explts <- purrr::map_dbl(tfleets, c("metiers", s, "p_explt"))

    p_total <- sum(p_explts)
    p_explts <- p_explts / ifelse(p_total > 0, p_total, 1e-6)

    for (f in fleet_names) {
      tfleets[[f]]$metiers[[s]]$p_explt <- as.numeric(p_explts[f])
    }
  }

  # --- Compute effective effort per fleet ---
  # Effort is a fleet-level quantity (identical across species); use the
  # first species as the reference for habitat overlap
  e_fl <- compute_effort_by_fleet(fauna[[fauni[1]]], tfleets)

  # --- Set initial catchability estimates ---
  for (s in fauni) {
    p_explt <- purrr::map_dbl(tfleets, c("metiers", s, "p_explt"))

    if (tune_type == "f") {
      explt_by_fleet <- fauna[[s]]$init_explt * p_explt
    } else {
      # Depletion mode: initial guess uses (1 - depletion) as total exploitation
      explt_by_fleet <- (1 - fauna[[s]]$depletion) * p_explt
    }

    raw_q <- explt_by_fleet / e_fl
    # Replace NaN/Inf from zero effort with 0
    raw_q[!is.finite(raw_q)] <- 0

    for (f in fleet_names) {
      set_metier_catchability(
        tfleets[[f]]$metiers[[s]],
        raw_q[f],
        use_link = (tune_type == "depletion")
      )
    }
  }

  # --- Tune costs ---
  if (tune_costs) {
    init_sim <- simmar(
      fauna = fauna,
      fleets = tfleets,
      years = years
    )

    eq <- init_sim[[length(init_sim)]]

    cost_calibration <- calibrate_fleet_costs(eq, tfleets)

    for (f in fleet_names) {
      row <- cost_calibration$fleet == f
      tfleets[[f]]$cost_per_unit_effort <- cost_calibration$cost_per_unit_effort[row]
      tfleets[[f]]$effort_reference <- cost_calibration$effort_reference[row]
    }
  }

  # --- Depletion-based tuning via numerical solver ---
  if (tune_type == "depletion") {

    # Phase 1: coarse solve
    log_fs_phase_1 <- nleqslv::nleqslv(
      x = rep(log(0.1), length(fauna)),
      fn = fleet_tuner,
      method = "Broyden",
      global = "dbldog",
      control = list(
        maxit = 80, xtol = 1e-3, ftol = 1e-3,
        allowSingular = TRUE, stepmax = 1.5
      ),
      fleets = tfleets,
      e_fl = e_fl,
      fauna = fauna,
      years = years
    )

    # Phase 2: refine from phase 1 solution
    log_fs <- nleqslv::nleqslv(
      x = log_fs_phase_1$x,
      fn = fleet_tuner,
      method = "Broyden",
      global = "dbldog",
      control = list(
        maxit = 150, xtol = 1e-5, ftol = 1e-5,
        allowSingular = TRUE, stepmax = 0.8, trace = 0
      ),
      fleets = tfleets,
      e_fl = e_fl,
      fauna = fauna,
      years = years
    )

    # --- Solver convergence diagnostic ---
    if (!(log_fs$termcd %in% c(1, 2))) {
      warning(
        "nleqslv solver did not converge (termination code: ",
        log_fs$termcd, "). ",
        "Depletion targets may not have been achieved. ",
        "Check whether target depletion is plausible given supplied ",
        "selectivities, fishing grounds, and p_explt.",
        call. = FALSE
      )
    }

    # Apply solved catchabilities through the link function
    for (f in seq_along(tfleets)) {
      for (ff in seq_along(fauna)) {

        f_critter <- exp(log_fs$x[ff])
        f_metier <- tfleets[[f]]$metiers[[ff]]$p_explt * f_critter
        raw_q <- f_metier / e_fl[f]

        if (!is.finite(raw_q)) raw_q <- 0

        set_metier_catchability(
          tfleets[[f]]$metiers[[ff]],
          raw_q,
          use_link = TRUE
        )
      }
    }

    # --- Post-solve validation: check achieved vs target depletion ---
    validation_sim <- simmar(
      fauna = fauna,
      fleets = tfleets,
      years = years
    )

    val_eq <- validation_sim[[length(validation_sim)]]

    val_ssb <- purrr::map(
      val_eq,
      ~ as.data.frame(.x$ssb_p_a)
    ) |>
      purrr::list_rbind(names_to = "fauna")

    val_ssb <- data.frame(
      fauna = val_ssb$fauna,
      ssb = rowSums(val_ssb[, 2:ncol(val_ssb)], na.rm = TRUE)
    ) |>
      dplyr::group_by(fauna) |>
      dplyr::summarise(ssb = sum(ssb)) |>
      dplyr::ungroup() |>
      dplyr::arrange(fauna)

    ssb0s <- purrr::map_dbl(fauna, "ssb0")
    ssb0s <- ssb0s[sort(names(ssb0s))]

    target_depletion <- purrr::map_dbl(fauna, "depletion")
    target_depletion <- target_depletion[sort(names(target_depletion))]

    achieved_depletion <- val_ssb$ssb / ssb0s

    depl_tol <- depl_tol # relative tolerance: warn if off
    rel_error <- abs(achieved_depletion - target_depletion) / target_depletion

    species_off <- which(rel_error > depl_tol)

    if (length(species_off) > 0) {
      species_names <- sort(names(fauna))
      msg_lines <- vapply(species_off, function(i) {
        sprintf(
          "  %s: target = %.3f, achieved = %.3f (%.1f%% off)",
          species_names[i],
          target_depletion[i],
          achieved_depletion[i],
          rel_error[i] * 100
        )
      }, character(1))

      warning(
        "Achieved depletion does not match target (>",
        depl_tol * 100, "% relative error) for:\n",
        paste(msg_lines, collapse = "\n"), "\n",
        "Consider adjusting base_effort, selectivity, fishing grounds, ",
        "or p_explt.",
        call. = FALSE
      )
    }

    # --- Refine cost calibration using actual equilibrium ---
    # The initial cost tuning used a heuristic (F = 1 - depletion), which is
    # imprecise. Now that we have the actual equilibrium from depletion tuning,
    # recalibrate costs for greater accuracy.
    if (tune_costs) {
      cost_calibration_refined <- calibrate_fleet_costs(val_eq, tfleets)

      for (f in fleet_names) {
        row <- cost_calibration_refined$fleet == f
        tfleets[[f]]$cost_per_unit_effort <- cost_calibration_refined$cost_per_unit_effort[row]
        tfleets[[f]]$effort_reference <- cost_calibration_refined$effort_reference[row]
      }
    }
  }

  # Restore original fleet models
  for (f in fleet_names) {
    tfleets[[f]]$fleet_model <- og_fleet_model[f]
  }

  tfleets
}

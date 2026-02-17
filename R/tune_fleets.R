#' Deep-clone all metiers in a fleet list
#'
#' Creates independent copies of all metier R6 objects so that
#' modifications during tuning do not mutate the caller's fleet objects.
#'
#' @param fleets A fleet list (output of \code{\link{create_fleet}}).
#'
#' @return A deep-cloned copy of \code{fleets} with independent metier objects.
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


#' Logistic link function mapping positive reals to (0, 1)
#'
#' Applies the transformation \code{q = raw_q / (1 + raw_q)}, equivalent to
#' \code{plogis(log(raw_q))}. This smoothly maps positive real values to
#' (0, 1), approximately preserving values much less than 1 (when
#' \code{raw_q << 1}, \code{q ≈ raw_q}) while preventing q from ever
#' reaching or exceeding 1. Used during numerical optimisation of
#' catchability to maintain smooth gradients near the constraint boundary.
#'
#' @param raw_q Positive numeric; unconstrained catchability value.
#'
#' @return Numeric in (0, 1).
#' @keywords internal
q_link <- function(raw_q) {
  raw_q / (1 + raw_q)
}


#' Inverse link function mapping (0, 1) back to positive reals
#'
#' Applies the transformation \code{raw_q = q / (1 - q)}, the inverse of
#' \code{\link{q_link}}.
#'
#' @param q Numeric in (0, 1).
#'
#' @return Positive numeric.
#' @keywords internal
q_link_inv <- function(q) {
  q / (1 - q)
}


#' Set catchability for a single fleet-species metier
#'
#' Assigns scalar catchability to a metier, rescales its spatial catchability
#' vector to have the correct mean, and rebuilds \code{vul_p_a}. Catchability
#' is enforced to lie in (0, 1).
#'
#' @param metier A \code{\link{Metier}} R6 object (modified in place via R6
#'   reference semantics).
#' @param q_value Numeric. Target scalar catchability. Must be in \code{[0, 1]}
#'   when \code{use_link = FALSE}. When \code{use_link = TRUE}, any
#'   non-negative value is accepted and mapped to (0, 1) via
#'   \code{\link{q_link}}.
#' @param use_link Logical. If \code{TRUE}, applies \code{\link{q_link}} to map
#'   \code{q_value} from the positive reals into (0, 1). Use during numerical
#'   optimisation to ensure smooth gradients. Default \code{FALSE}.
#'
#' @return The modified metier, invisibly.
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
#' For each fleet, calculates the share of \code{base_effort} that overlaps
#' with viable habitat (\code{b0_p > 0}) within the fleet's fishing grounds.
#' Effort is a fleet-level quantity (not species-specific): a fleet exerts one
#' pool of effort that affects each species differently through catchability.
#'
#' @param fauna_species A single fauna element (used only for \code{b0_p}).
#' @param fleets A fleet list (output of \code{\link{create_fleet}}).
#'
#' @return Named numeric vector of effective effort per fleet.
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
#' Extracts revenue and effort from an equilibrium time step, then solves for
#' \code{cost_per_unit_effort} and \code{effort_reference} for each fleet so
#' that costs match the target cost-to-revenue ratio.
#'
#' The calibration formula is:
#' \deqn{c_0 = \frac{\text{cr\_ratio} \times R}{E^{ref} \times P \times (1 + \theta)}}
#'
#' The effort cost exponent \eqn{\gamma} drops out because tuning runs under
#' constant effort, distributing effort roughly uniformly (\eqn{E_l / E^{ref}
#' \approx 1}, so \eqn{1^\gamma = 1}).
#'
#' @param eq A single equilibrium time step from \code{\link{simmar}} output
#'   (i.e. \code{sim[[length(sim)]]}).
#' @param fleets A fleet list (output of \code{\link{create_fleet}}).
#'
#' @return A data frame with columns: \code{fleet},
#'   \code{cost_per_unit_effort}, \code{effort_reference}, \code{revenue},
#'   \code{total_effort}, \code{n_open_patches}, \code{cr_ratio},
#'   \code{travel_weight}, \code{implied_cr}.
#' @keywords internal
calibrate_fleet_costs <- function(eq, fleets) {

  fleet_names <- names(fleets)

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

  fleet_params <-
    purrr::map(
      fleets,
      ~ data.frame(
        cr_ratio = .x$cr_ratio,
        travel_weight = .x$travel_weight
      )
    ) |>
    purrr::list_rbind(names_to = "fleet")

  fishing_grounds <-
    purrr::map(
      fleets,
      ~ data.frame(
        patch = seq_along(.x$fishing_grounds$fishing_ground),
        fishing_ground = .x$fishing_grounds$fishing_ground
      )
    ) |>
    purrr::list_rbind(names_to = "fleet")

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


#' Tune Fleet Catchability and Costs to Target Initial Conditions
#'
#' @description
#' Adjusts the catchability coefficient of each fleet's metiers so that the
#' simulated fishery reaches either a target fishing mortality rate or a target
#' depletion level (B/B0). Optionally tunes \code{cost_per_unit_effort} for
#' each fleet to match a target cost-to-revenue ratio (\code{cr_ratio}).
#'
#' This function should be called after \code{\link{create_fleet}} and before
#' passing fleets to \code{\link{simmar}}; the tuned fleet list is returned
#' and should replace the original.
#'
#' @details
#' ## Tuning type
#' \describe{
#'   \item{\code{"f"} (or \code{"explt"})}{Sets catchability directly from each
#'     critter's \code{init_explt} and the metier's \code{p_explt} share.
#'     Analytical; no numerical solver required. Fast but does not guarantee a
#'     precise equilibrium depletion.}
#'   \item{\code{"depletion"}}{Uses a two-phase Broyden solver via
#'     \code{nleqslv} to find catchabilities that produce the target
#'     depletion specified in each critter's \code{depletion} field. A
#'     logistic link function keeps catchabilities in (0, 1) throughout
#'     optimisation. Post-solve validation warns if achieved depletion differs
#'     from the target by more than \code{depl_tol}.}
#' }
#'
#' When \code{tune_costs = TRUE} and \code{tune_type = "depletion"}, a
#' two-pass cost calibration is used: an initial heuristic pass before the
#' depletion solver, then a refined pass using the actual equilibrium. This
#' ensures costs accurately reflect final equilibrium conditions.
#'
#' Note that tuning is approximate: post-tuning values will not perfectly match
#' targets because some calibration steps depend on earlier steps.
#'
#' @param fauna A named list of fauna objects from \code{\link{create_critter}}.
#' @param fleets A named list of fleet objects from \code{\link{create_fleet}}.
#' @param years Integer. Number of years to simulate during tuning runs.
#'   Longer values ensure equilibrium is reached. Default \code{50}.
#' @param tune_type Character. One of \code{"f"} / \code{"explt"} (tune to
#'   fishing mortality rate) or \code{"depletion"} (tune to B/B0). See Details.
#' @param tune_costs Logical. If \code{TRUE} (default), calibrate
#'   \code{cost_per_unit_effort} for each fleet so that equilibrium costs
#'   match the fleet's \code{cr_ratio}.
#' @param depl_tol Numeric. Relative tolerance for the depletion-tuning
#'   convergence check. A warning is issued if any species' achieved depletion
#'   differs from the target by more than this fraction. Default \code{0.025}
#'   (2.5\%).
#'
#' @return A tuned copy of \code{fleets} with updated \code{catchability},
#'   \code{spatial_catchability}, \code{vul_p_a}, and (if
#'   \code{tune_costs = TRUE}) \code{cost_per_unit_effort} and
#'   \code{effort_reference} fields for each metier.
#'
#' @seealso \code{\link{create_fleet}}, \code{\link{simmar}},
#'   \code{\link{create_critter}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fauna  <- list(tuna = create_critter("Thunnus albacares",
#'                                      resolution = c(10, 10),
#'                                      fished_depletion = 0.4))
#' fleets <- list(fleet = create_fleet(metiers   = list(tuna = met),
#'                                     resolution = c(10, 10)))
#'
#' # Tune catchability to achieve target depletion
#' fleets <- tune_fleets(fauna, fleets, tune_type = "depletion")
#'
#' # Then run the simulation
#' sim <- simmar(fauna = fauna, fleets = fleets, years = 50)
#' }
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

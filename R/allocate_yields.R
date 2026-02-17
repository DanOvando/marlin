allocate_yields <- function(f_p_a_fl,p_p_a_fl, e_p_fl, critter, pop, patches, ages, fleets, fauna) {
  c_p_a_fl <- f_p_a_fl * array(
    pop$c_p_a,
    dim = c(patches, ages, length(fleets)),
    dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
  )


  fleet_names <- names(fleets)
  r_p_a_fl <- c_p_a_fl * p_p_a_fl

  tmp_e_p_fl <- e_p_fl


  c_p_fl <- matrix(nrow = patches, ncol = length(fleets), dimnames = list(NULL, fleet_names))
  r_p_fl <- matrix(nrow = patches, ncol = length(fleets), dimnames = list(NULL, fleet_names))
  prof_p_fl <- matrix(nrow = patches, ncol = length(fleets), dimnames = list(NULL, fleet_names))
  cost_p_fl <- matrix(nrow = patches, ncol = length(fleets), dimnames = list(NULL, fleet_names))


  for (fl in 1:length(fleets)) {

    c_p_fl[, fl] <- rowSums(c_p_a_fl[, , fl], na.rm = TRUE)

    r_p_fl[, fl] <- rowSums(r_p_a_fl[, , fl], na.rm = TRUE)

    # New normalized cost formula:
    # Cost_f = c0_f * E_ref_f * [sum_l (E_l,f / E_ref_f)^gamma_f + theta_f * sum_l d_tilde_l,f * (E_l,f / E_ref_f)]
    #
    # This is called once per critter in the loop, so we keep the /length(fauna) hack
    # to divide effort by number of species (it sums correctly across the loop)
    #
    # Breaking it down:
    # - e_normalized = (e_p / n_species) / effort_reference
    # - convex_term = sum(e_normalized^gamma)
    # - travel_term = theta * sum(cost_per_patch_normalized * e_normalized)
    # - total_cost = c0 * effort_reference * (convex_term + travel_term)

    e_ref <- fleets[[fl]]$effort_reference
    gamma <- fleets[[fl]]$effort_cost_exponent
    theta <- fleets[[fl]]$travel_weight
    c0 <- fleets[[fl]]$cost_per_unit_effort
    d_tilde <- fleets[[fl]]$cost_per_patch_normalized

    # Normalize effort by reference level (with /n_species hack)
    e_normalized <- (tmp_e_p_fl[, fl] / length(fauna)) / e_ref

    # Convex congestion term: (E_l / E_ref)^gamma
    convex_term <- e_normalized^gamma

    # Travel term: theta * d_tilde * (E_l / E_ref)
    travel_term <- theta * d_tilde * e_normalized

    # Total cost per patch: c0 * E_ref * (convex + travel)
    # Note: We multiply by e_ref here to convert the normalized costs back to actual units
    cost_p_fl[, fl] <- c0 * e_ref * (convex_term + travel_term)

    prof_p_fl[, fl] <-
      r_p_fl[, fl] - cost_p_fl[,fl]

  }

  out <- list(
    c_p_a_fl = c_p_a_fl,
    r_p_a_fl = r_p_a_fl,
    r_p_fl = r_p_fl,
    c_p_fl = c_p_fl,
    prof_p_fl = prof_p_fl
  )
} # close allocate_catch

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

    cost_p_fl[,fl] <- fleets[[fl]]$cost_per_unit_effort * ((
      as.matrix((tmp_e_p_fl[, fl]) / length(fauna))^fleets[[fl]]$effort_cost_exponent
    ) + as.matrix(tmp_e_p_fl[, fl] / length(fauna) * fleets[[fl]]$cost_per_patch)
    )

    prof_p_fl[, fl] <-
      r_p_fl[, fl] - cost_p_fl[,fl]

  }

  out <- list(
    c_p_a_fl = c_p_a_fl,
    r_p_a_fl = r_p_a_fl,
    r_p_fl = r_p_fl,
    c_p_fl = c_p_fl,
    cost_p_fl = cost_p_fl,
    prof_p_fl = prof_p_fl
  )
} # close allocate_catch

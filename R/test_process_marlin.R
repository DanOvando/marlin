fast_process_marlin <- function(sim, steps_to_keep = NA, time_step = NA, keep_age = TRUE) {
  stopifnot(length(sim) >= 1)
  nm <- sub("^step_", "", names(sim), perl = TRUE)
  names(sim) <- nm
  if (all(is.na(steps_to_keep))) steps_to_keep <- nm
  sim <- sim[as.character(steps_to_keep)]

  years   <- as.integer(sub("_.*$", "", nm))
  seasons <- as.integer(sub("^.*_", "", nm))

  resolution <- sim[[1]][[1]]$resolution
  n_patch <- resolution[1] * resolution[2]

  # coords once
  coords <- data.table::data.table(
    patch = seq_len(n_patch),
    x = rep(seq_len(resolution[1]), times = resolution[2]),
    y = rep(seq_len(resolution[2]), each  = resolution[1])
  )

  # infer time_step if needed
  if (is.na(time_step)) {
    if (length(sim) > 1) {
      time_step <- 1 / data.table::uniqueN(seasons)
    } else {
      time_step <- 1
      warning("Assuming time_step = 1; set `time_step` explicitly for single-step sims.")
    }
  }

  library(data.table)

  # --- fauna ---
  tidy_one_species <- function(y, critter_id) {
    # y$n_p_a etc are patch x age matrices
    mats <- list(
      n = y$n_p_a, b = y$b_p_a, ssb = y$ssb_p_a, c = y$c_p_a
    )

    # Melt all metrics, then join wide by (patch, age)
    long_list <- lapply(names(mats), function(m) {
      dt <- as.data.table(mats[[m]])
      dt[, patch := .I]
      mdt <- melt(dt, id.vars = "patch", variable.name = "age", value.name = m)
      mdt[, age := as.integer(sub("^V", "", as.character(age)))]
      setkey(mdt, patch, age)
      mdt
    })
    fauna <- Reduce(function(a, b) a[b], long_list)  # keyed join on (patch, age)

    if (!keep_age) {
      fauna <- fauna[, .(n = sum(n), b = sum(b), ssb = sum(ssb), c = sum(c)), by = .(patch)]
      fauna[, age := "all"]
    }

    fauna <- coords[fauna, on = "patch"]
    fauna[, critter := critter_id]
    if (keep_age) fauna[, age := y$ages[age]]  # map to actual ages if needed
    fauna[]
  }

  step_fauna <- function(step_list) {
    # step_list is a named list of species
    rbindlist(Map(tidy_one_species, step_list, names(step_list)), use.names = TRUE, fill = TRUE)
  }

  fauna_all <- rbindlist(
    Map(step_fauna, sim),
    use.names = TRUE, fill = TRUE, idcol = "step"
  )

  # add time columns without regex
  fauna_all[, year   := years[match(step, nm)]]
  fauna_all[, season := seasons[match(step, nm)]]
  fauna_all[, step   := year + (season * time_step - time_step)]

  # mean length lookup
  get_len <- function(x, nm) data.table(critter = nm, age = x$ages, mean_length = x$length_at_age)
  len_lookup <- rbindlist(Map(get_len, sim[[1]], names(sim[[1]])))
  if (keep_age) {
    fauna_all <- len_lookup[fauna_all, on = .(critter, age)]
  } else {
    fauna_all[, mean_length := NA_real_]
  }

  # --- fleets ---
  tidy_fleet_one <- function(sp, critter_id) {
    # c_p_a_fl and r_p_a_fl are patch x age x fleet arrays
    catch <- as.data.table(melt(sp$c_p_a_fl, value.name = "catch",
                                varnames = c("patch","age","fleet")))
    revenue <- as.data.table(melt(sp$r_p_a_fl, value.name = "revenue",
                                  varnames = c("patch","age","fleet")))
    setkey(catch, patch, age, fleet)
    setkey(revenue, patch, age, fleet)
    dt <- catch[revenue]

    if (!keep_age) {
      dt <- dt[, .(catch = sum(catch), revenue = sum(revenue)), by = .(patch, fleet)]
      dt[, age := "all"]
    }

    # effort: matrix patch x fleet (columns are fleets)
    eff <- as.data.table(sp$e_p_fl)
    eff[, patch := .I]
    eff <- melt(eff, id.vars = "patch", variable.name = "fleet", value.name = "effort")
    eff[, fleet := as.character(fleet)]  # keep as character to match above

    setkey(dt, patch, fleet)
    setkey(eff, patch, fleet)
    dt <- dt[eff]

    dt <- coords[dt, on = "patch"]
    dt[, critter := critter_id]
    dt[, cpue := fifelse(effort > 0, catch / effort, NA_real_)]
    dt[]
  }

  step_fleet <- function(step_list) {
    rbindlist(Map(tidy_fleet_one, step_list, names(step_list)), use.names = TRUE, fill = TRUE)
  }

  fleets_all <- rbindlist(
    Map(step_fleet, sim),
    use.names = TRUE, fill = TRUE, idcol = "step"
  )
  fleets_all[, year   := years[match(step, nm)]]
  fleets_all[, season := seasons[match(step, nm)]]
  fleets_all[, step   := year + (season * time_step - time_step)]

  if (keep_age) {
    fleets_all <- len_lookup[fleets_all, on = .(critter, age)]
  } else {
    fleets_all[, mean_length := NA_real_]
  }

  list(fauna = as_tibble(fauna_all),
       fleets = as_tibble(fleets_all))
}

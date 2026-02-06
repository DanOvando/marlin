#' simmar is the wrapper function for the marlin package
#'
#' when passed fauna and fleet objects, simmar will advance
#' the population for a number of steps
#'
#' @param fauna a list of fauna objects
#' @param fleets a list of fleet objects
#' @param habitat a list of habitat over time
#' @param years the number of years to run the simulation
#' @param initial_conditions initial conditions for the simulation, in the form simmar()[[final step]]
#' @param starting_step  the step to start the simulation from, used to keep track of steps across multiple runs of simmar
#' @param keep_starting_step should the starting step by kept (TRUE) or dropped (FALSE)
#' @param manager a list of management actions
#' @param steps the number of steps to run, as an alternative to years
#' @param starting_season the starting season for the simulation
#' @param log_rec_devs externally supplied log recruitment deviates
#' @param cor_rec correlation matrix in recruitment deviates across species
#'
#' @return a list containing the results of the simulation
#' @export
#'
simmar <- function(fauna = list(),
                   fleets = list(),
                   manager = list(),
                   habitat = list(),
                   years = 100,
                   steps = NA,
                   starting_season = NA,
                   initial_conditions = NA,
                   starting_step = NA,
                   keep_starting_step = TRUE,
                   log_rec_devs = NULL,
                   cor_rec = diag(length(fauna))) {

  init_cond_provided <-
    !all(is.na(initial_conditions)) # marker in case initial conditions were provided

  fauni <- names(fauna)

  # if initial_conditions is passed, enforce names/order too (cheap and catches a ton)
  if (!all(is.na(initial_conditions))) {
    initial_conditions <- rlang::set_names(initial_conditions, fauni)
    initial_conditions <- initial_conditions[fauni]
  }

  if (is.null(fauni) || anyNA(fauni) || any(fauni == "") || any(duplicated(fauni))) {
    stop("fauna must be a *named* list with unique, non-empty species names.")
  }

  n_critters <- length(fauna)

  fleet_names <- names(fleets)

  time_step <- unique(purrr::map_dbl(fauna, "time_step"))

  steps_per_year <- 1 / time_step

  patch_area <- unique(purrr::map_dbl(fauna, "patch_area"))

  sigma_recs <- purrr::map_dbl(fauna, "sigma_rec") # gather recruitment standard deviations

  ac_recs <- purrr::map_dbl(fauna, "ac_rec") # gather autocorrelation in recruitment standard deviations

  if (any(abs(ac_recs) >= 1)) {
    stop("all sigma_ac values must be between -1 and 1")
  }

  covariance_rec <- cor_rec * (sigma_recs %o% sigma_recs)

  # if (!is.null(rec_devs)){
  #   if (length(rec_devs) > 1 & length(rec_devs) != length(fauna)){
  #     stop("rev_devs must either be a vector of length")
  #   }
  # }

  if (length(time_step) > 1) {
    stop(
      paste(
        "All critters in fauna must have the same time step: current time steps are",
        paste(time_step, collapse = " ")
      )
    )
  }

  if (is.na(steps)) {
    steps <-
      (years) / time_step + 1 # tack on extra step for accounting
  } else {
    steps <- steps + 1 # to store initial conditions
  } # if steps are specified instead of years

  steps <- pmax(steps, 3) # need to be at least 1 initial + 2 running steps

  # generate recruitment deviates
  if (!is.na(starting_step)) {
    year_season <- marlin::clean_steps(starting_step)

    starting_year <- as.integer(sub("_.*$", "", year_season)) - 1
    starting_season <- as.integer(sub("^.*_", "", year_season)) - 1

    offset <- (starting_year * steps_per_year) + starting_season
  } else {
    offset <- 0
  }

  # create step names of form "year_season"
  step_names <- paste(
    rep(1:(steps + offset), each = steps_per_year),
    1:steps_per_year,
    sep = "_"
  )
  step_names <- step_names[(1:steps) + offset]

  # --- SPEED: parse year/season ONCE (avoid repeated stringr in the step loop) ---
  step_year   <- as.integer(sub("_.*$", "", step_names))
  step_season <- as.integer(sub("^.*_", "", step_names))
  # --- end step parsing ---

  # recruitment deviates
  if (is.null(log_rec_devs)) {
    log_rec_devs <- matrix(
      NA_real_,
      nrow = length(step_names),
      ncol = n_critters,
      dimnames = list(step_names, fauni)
    )

    # SPEED: draw all MVN innovations once, then do AR(1) recursion
    eps <- mvtnorm::rmvnorm(length(step_names), mean = rep(0, n_critters), sigma = covariance_rec)
    log_rec_devs[1, ] <- eps[1, ]
    if (length(step_names) > 1) {
      ar_scale <- sqrt(1 - ac_recs^2)
      for (i in 2:length(step_names)) {
        log_rec_devs[i, ] <- ac_recs * log_rec_devs[i - 1, ] + ar_scale * eps[i, ]
      }
    }
  } else {

    if (is.null(colnames(log_rec_devs)) ||
        !setequal(colnames(log_rec_devs), names(fauna))) {
      stop("column names of log_rec_devs must match those of fauna")
    }
    log_rec_devs <- log_rec_devs[, names(fauna), drop = FALSE]


  }

  # step_names <- step_names[1:steps] # chop back to the actual number of steps used; works this way in case there are multiple steps per year but an incomplete number of years
  patches <- unique(purrr::map_dbl(fauna, "patches"))

  # not building checks for same resolution here since almost redundant to patches unless product of resolution is identical
  resolution <-
    (purrr::map(fauna[1], ~ data.frame(t(.x$resolution)))) |>
    purrr::list_rbind() |>
    unlist()

  if (all(is.na(initial_conditions))) {
    initial_conditions <-
      purrr::map(fauna, c("unfished")) # pull out unfished conditions created by create_critter
  }

  if (length(patches) > 1) {
    stop(
      "fauna have different habitat resolutions: set resolution to same number for all species!"
    )
  }

  #  determine open fishing seasons

  season_foo <- function(fauni, name, manager) {
    open_seasons <- 1:fauni$seasons

    if (length(manager$closed_seasons[name]) > 0) {
      open_seasons <-
        open_seasons[!(open_seasons %in% manager$closed_seasons[[name]])]
    }

    return(open_seasons)
  }

  fishing_seasons <-
    purrr::imap(fauna, season_foo, manager = manager) # figure out which seasons are open for each critter

  if (!is.null(names(fishing_seasons)) && all(fauni %in% names(fishing_seasons))) {
    fishing_seasons <- fishing_seasons[fauni]
  }

  # --- Habitat handling (same logic, faster inner use) ---
  # We keep your expansion logic exactly, but ALSO precompute the per-step habitat vector
  # in the SAME patch order you were enforcing with:
  #   pivot_longer(as.data.frame(mat), everything(), names_to="x", names_transform=list(x=as.integer)) |> arrange(x)
  #
  # NOTE (spirit-of-comment): because pivot_longer over a data.frame stacks *columns* in order,
  # and arrange(x) doesn’t change that column order, the resulting $value is identical to as.vector(mat).
  habitat_vecs <- list()

  if (length(habitat) > 0) {
    if (!any(names(habitat) %in% fauni)) {
      stop(
        "names of critters in habitat must must match name of at least one critter in fauni list"
      )
    }

    for (f in seq_along(habitat)) {
      if (length(habitat[[f]]) != years &
          length(habitat[[f]]) != (steps - 1)) {
        stop(
          "supplied habitat vector must be same length as either number of years or number of years times number of seasons"
        )
      } else if (length(habitat[[f]]) != (steps - 1)) {
        new_habitat <- vector(mode = "list", length = steps - 1)

        for (i in 1:(steps - 1)) {
          supplied_years <- length(habitat[[f]])
          if (!is.na(starting_step)) {
            years_in <-
              as.integer(sub("_.*$", "", step_names[i])) - as.integer(sub("_.*$", "", starting_step)) + 1
            # put year in index form not named form
          } else {
            years_in <-
              as.integer(sub("_.*$", "", step_names[i]))
          }
          years_in <- min(years_in, supplied_years)

          # year <-  floor(step_names[i] - starting_step) # put year in index form not named form
          new_habitat[[i]] <- habitat[[f]][[years_in]]
        }

        habitat[[f]] <- new_habitat
      } # close expansion of habitat if needed
    } # close fauni loop

    # Precompute vectors once (speeds up the step loop)
    for (sp in intersect(names(habitat), fauni)) {
      habitat_vecs[[sp]] <- lapply(habitat[[sp]], function(hmat) {
        as.vector(hmat) # equivalent ordering to your pivot_longer(... ) |> arrange(x)
      })
    }
  } # close check on supplied habitat
  # --- end habitat handling ---

  storage <- replicate(
    steps,
    setNames(vector("list", length(fauna)), names(fauna)),
    simplify = FALSE
  )

  storage[[1]] <-
    initial_conditions # start populations at initial conditions


  fleets <- purrr::map(
    fleets,
    ~ purrr::list_modify(.x, e_p_s = matrix((.x$base_effort / patches), nrow = patches, ncol = steps))
  )

  if (init_cond_provided) {
    for (i in names(fleets)) {
      fleets[[i]]$e_p_s[, 1] <- getElement(initial_conditions[[1]]$e_p_fl, i)
    }
  }

  r_p_f <- matrix(0, patches, length(fauni))
  e_p_f <- matrix(0, patches, length(fauni))
  f_q <- rep(0, length(fauni))

  c_p_fl <- matrix(nrow = patches, ncol = length(fleets), dimnames = list(NULL, fleet_names))
  r_p_fl <- matrix(nrow = patches, ncol = length(fleets), dimnames = list(NULL, fleet_names))
  prof_p_fl <- matrix(nrow = patches, ncol = length(fleets), dimnames = list(NULL, fleet_names))

  fishable <- rep(1, patches)

  # loop over steps
  for (s in 2:steps) {

    last_season <- step_season[s - 1]
    year <- step_year[s]
    current_season <- step_season[s]

    # NOTE: starting_season is in signature but not used in your original script body here;
    # leaving behavior unchanged.

    if (length(manager$mpas) > 0) {
      manager$mpas$locations <- manager$mpas$locations |>
        dplyr::arrange(x, y)

      if (year == manager$mpas$mpa_year) {
        fishable <- manager$mpas$locations$mpa == 0
      }
    }

    fleet_fishable <- vector(mode = "list", length = length(fleet_names))
    fleet_concentrator <- vector(mode = "list", length = length(fleet_names))

    for (l in seq_along(fleet_names)) {
      tmp_fleet_fishable <- fishable

      if (!is.null(fleets[[l]]$fishing_grounds)) {
        tmp_fleet_fishable[fleets[[l]]$fishing_grounds$fishing_ground == FALSE] <- 0
      }

      fleet_fishable[[l]] <- tmp_fleet_fishable

      concentrator <- rep(1, patches)
      if (length(manager$mpas) > 0) {
        if (year >= manager$mpas$mpa_year && fleets[[l]]$mpa_response == "leave") {
          concentrator <- as.numeric(fishable)
        }
      }
      fleet_concentrator[[l]] <- concentrator
    }

    for (l in seq_along(fleet_names)) {
      # xx add ability to incorporate past revenues here. Idea. have a marker that lets you know if initial conditions were passed in. If s<= 2 or there were no initial conditions, pull this. If s<= 2 but there were initial conditions, pull the initial conditions for those steps

      if (s == 2 & init_cond_provided) {
        r_p_f <-
          (sapply(initial_conditions, function(x) {
            rowSums(x$r_p_a_fl[, , l], na.rm = TRUE)
          }))

        last_r_p <-
          rowSums(r_p_f, na.rm = TRUE) # pull out total revenue for fleet l
      } else if (s <= 2) {
        for (f in seq_along(fauni)) {

          critter <- fauni[f]

          last_b_p_a <- storage[[s - 1]][[critter]]$b_p_a

          last_e_p <- fleets[[l]]$e_p_s[, s - 1]

          # calculate fishable biomass in each patch for each species for that fleet

          # account for spatial catchability
          tmp <- 1 - exp(-(
            matrix(
              fleets[[l]]$metiers[[critter]]$spatial_catchability,
              nrow = nrow(last_b_p_a),
              ncol = ncol(last_b_p_a),
              byrow = FALSE
            ) *
              matrix(
                fleets[[l]]$metiers[[critter]]$sel_at_age,
                nrow = nrow(last_b_p_a),
                ncol = ncol(last_b_p_a),
                byrow = TRUE
              )
          ))

          last_b_p <-
            rowSums(last_b_p_a * tmp * (current_season %in% fishing_seasons[[critter]])) * fleet_fishable[[l]]

          r_p_f[, f] <-
            last_b_p * fleets[[l]]$metiers[[critter]]$price

          f_q[f] <- fleets[[l]]$metiers[[critter]]$catchability
        } # close fauni loop

        last_r_p <- rowSums(r_p_f, na.rm = TRUE)
      } else {
        r_p_f <-
          (sapply(storage[[s - 2]], function(x) {
            rowSums(x$r_p_a_fl[, , l], na.rm = TRUE)
          }))

        last_r_p <-
          rowSums(r_p_f, na.rm = TRUE) # pull out total revenue for fleet l
      }

      if (fleets[[l]]$fleet_model == "manual") {
        total_effort <- fleets[[l]]$effort[s - 1]
      } else {
        total_effort <- sum(fleets[[l]]$e_p_s[, s - 1] * fleet_concentrator[[l]])
      }

      if (sum(last_r_p, na.rm = TRUE) > 0) {
        # the only way for revenues in the last step to be literally zero is 100% mpas or all target species seasons closed or last effort = zero
        last_revenue <-
          sum(last_r_p, na.rm = TRUE) # pull out total revenue for fleet l

        last_cost <-
          fleets[[l]]$cost_per_unit_effort * (
            sum((fleets[[l]]$e_p_s[, s - 1])^
                  fleets[[l]]$effort_cost_exponent) + sum(fleets[[l]]$cost_per_patch * fleets[[l]]$e_p_s[, s - 1])
          )
      }

      if (fleets[[l]]$fleet_model == "open_access") {
        if (is.na(fleets[[l]]$cost_per_unit_effort) || is.na(fleets[[l]]$responsiveness)) {
          stop("open access fleet model requires both cost_per_unit_effort and responsiveness parameters")
        }

        if (exists("last_revenue")) {
          effort_cap <- Inf
          if (length(manager$effort_cap[[l]]) > 0) effort_cap <- manager$effort_cap[[l]]

          total_effort <- pmin(
            effort_cap,
            total_effort * pmin(1.5, exp(
              fleets[[l]]$responsiveness * log(pmax(last_revenue, 1e-6) / pmax(1e-6, last_cost))
            ))
          )
        }
      }

      e_p <- fleets[[l]]$e_p_s[, s - 1]


      e_floor <- 1e-2 * mean(e_p)

      # allocate fleets in space
      if (fleets[[l]]$spatial_allocation == "revenue") {
        if (sum(fleet_fishable[[l]]) == 0) {
          alloc <- 0
        } else if (sum(last_r_p, na.rm = TRUE) == 0) {
          alloc <- fleet_fishable[[l]] / sum(fleet_fishable[[l]])
        } else {
          alloc <- ((last_r_p * fleet_fishable[[l]]) / sum(last_r_p * fleet_fishable[[l]], na.rm = TRUE))
          alloc <- alloc - min(alloc, na.rm = TRUE) + 1
          alloc <- alloc / sum(alloc)

          if (s > 2) {
            alloc <- rowSums(sapply(storage[[s - 2]], function(x) x$r_p_fl[, l]), na.rm = TRUE)
            alloc <- alloc - min(alloc, na.rm = TRUE) + 1
            alloc <- alloc * fleet_fishable[[l]]
            alloc <- alloc / sum(alloc)
          }
        }

        fleets[[l]]$e_p_s[, s] <- total_effort * alloc
      } else if (fleets[[l]]$spatial_allocation == "uniform"){

        if (sum(fleet_fishable[[l]]) == 0) {
          alloc <- 0
        } else {
          alloc <- fleet_fishable[[l]] / sum(fleet_fishable[[l]])
        }

      } else if (fleets[[l]]$spatial_allocation == "rpue") {
        if (sum(fleet_fishable[[l]]) == 0) {
          alloc <- 0
        } else if (sum(last_r_p, na.rm = TRUE) == 0) {
          alloc <- fleet_fishable[[l]] / sum(fleet_fishable[[l]])
        } else {

          if (s<=2){
          if (init_cond_provided) {
            rpue <-  rpue_bar <- (last_r_p / pmax(e_p, e_floor))
          } else {
            rpue <- rpue_bar <-  fleet_fishable[[l]] / sum(fleet_fishable[[l]])
          }
          } else {
            rpue <- rowSums(sapply(storage[[s - 2]], function(x) x$r_p_fl[, l]), na.rm = TRUE) / pmax(e_p, e_floor)
          }
        }
        rpue[!is.finite(rpue)] <- 0
        rpue <- rpue * fleet_fishable[[l]]

        smoother <- 0.4

        rpue_bar <- (1-smoother)*rpue_bar + smoother*rpue

        beta <- 0.1

        z <- rpue_bar - median(rpue_bar)
        z <- z / (mad(rpue_bar) + 1e-12)   # robust scale
        w <- exp(beta * (z - max(z)))

        # w <- exp(beta * (rpue_bar - max(rpue_bar)))
        w <- w * fleet_fishable[[l]]
        if (sum(w) == 0) w <- fleet_fishable[[l]]
        alloc <- w / sum(w)

        fleets[[l]]$e_p_s[, s] <- total_effort * alloc
      } else if (fleets[[l]]$spatial_allocation == "ppue" && !is.na(fleets[[l]]$cost_per_unit_effort)) {
        if (sum(fleet_fishable[[l]]) == 0) {
          alloc <- 0
        } else if (sum(last_r_p, na.rm = TRUE) == 0) {
          alloc <- fleet_fishable[[l]] / sum(fleet_fishable[[l]])
        } else {
          alloc <- ((
            last_r_p - fleets[[l]]$cost_per_unit_effort * ((e_p)^fleets[[l]]$effort_cost_exponent + fleets[[l]]$cost_per_patch * e_p) / (e_p + 1)
          )) * fleet_fishable[[l]]

          alloc[!is.finite(alloc)] <- 0
          alloc <- alloc - min(alloc, na.rm = TRUE) + 1
          alloc <- alloc / sum(alloc, na.rm = TRUE)

          if (s > 2) {
            alloc <- rowSums(sapply(storage[[s - 2]], function(x) x$prof_p_fl[, l]), na.rm = TRUE) / (e_p)
            alloc[!is.finite(alloc)] <- 0
            alloc <- fleet_fishable[[l]] * (alloc - min(alloc, na.rm = TRUE) + 1)
            alloc <- alloc / sum(alloc)
          }
        }

        fleets[[l]]$e_p_s[, s] <- total_effort * alloc
      } else if (fleets[[l]]$spatial_allocation == "profit" && !is.na(fleets[[l]]$cost_per_unit_effort)) {
        if (sum(fleet_fishable[[l]]) == 0) {
          alloc <- 0
        } else if (sum(last_r_p, na.rm = TRUE) == 0) {
          alloc <- fleet_fishable[[l]] / sum(fleet_fishable[[l]])
        } else {
          alloc <- (last_r_p - fleets[[l]]$cost_per_unit_effort * ((e_p)^fleets[[l]]$effort_cost_exponent + fleets[[l]]$cost_per_patch * e_p))
          alloc[!is.finite(alloc)] <- 0
          alloc <- fleet_fishable[[l]] * (alloc - min(alloc, na.rm = TRUE) + 1)
          alloc <- alloc / sum(alloc, na.rm = TRUE)

          if (s > 2) {
            alloc <- rowSums(sapply(storage[[s - 2]], function(x) x$prof_p_fl[, l]), na.rm = TRUE)
            alloc <- fleet_fishable[[l]] * (alloc - min(alloc, na.rm = TRUE) + 1)
            alloc <- alloc / sum(alloc)
          }
        }

        fleets[[l]]$e_p_s[, s] <- total_effort * alloc
      } else if (fleets[[l]]$spatial_allocation == "manual") {
        alloc <- fleet_fishable[[l]] * fleets[[l]]$fishing_grounds$fishing_ground
        alloc <- alloc / sum(alloc)
        fleets[[l]]$e_p_s[, s] <- total_effort * alloc
      } else if (fleets[[l]]$spatial_allocation == "ifd") {
        # leaving your IFD branch intact; optimization here is mostly upstream (habitat + step parsing)
        cost_patch <- fleets[[l]]$cost_per_patch
        c0 <- 2
        gamma <- 1.2
        fishable_int <- as.integer(fleet_fishable[[l]])

        pre <- precompute_baranov_inputs(
          storage = storage[[s - 1]],
          fauna = fauna,
          fleets = fleets,
          target_fleet = l,
          E_exo = NULL,
          P = length(patches) # IMPORTANT: preserve your original behavior (even though patches is scalar)
        )

        out <- cpp_allocate_ifd_kkt_fullsolve_fast(
          Etot_target = sum(e_p, na.rm = TRUE),
          alpha_mats = pre$alpha_mats,
          other_mort_mats = pre$other_mort_mats,
          biomass_mats = pre$biomass_mats,
          price_s = pre$price_s,
          cost_patch = cost_patch,
          c0 = c0,
          gamma = gamma,
          fishable_int = fishable_int,
          time_step = time_step,
          include_costs = TRUE,
          n_outer = 24L,
          n_inner = 18L,
          active_tol = 1e-14,
          flat_tol_sd = 1e-6,
          flat_tol_abs = 1e-10
        )

        alloc <- out$E_target / sum(out$E_target)
        fleets[[l]]$e_p_s[, s] <- total_effort * alloc
      } else if (fleets[[l]]$spatial_allocation == "marginal_profits") {


        pre <- precompute_baranov_inputs_softmax(
          storage[[s - 1]],
          fauna,
          fleets,
          l,
          E_exo = purrr::map(fleets, \(x,s) x$e_p_s[,s-1],s = s),
          (patches)
        )

        res <- allocate_until_stable_patchwise(
          eff_init_p = e_p, E_tot = sum(e_p),
          v_mats = pre$v_mats, other_mats = pre$other_mats, b_mats = pre$b_mats, price_s = pre$price_s, time_step,
          c0 = fleets[[l]]$cost_per_unit_effort, gamma = fleets[[l]]$effort_cost_exponent, travel_p = fleets[[l]]$cost_per_patch, open_p = (fleet_fishable[[l]] == 1),
          beta = 6, rho = 0.1,
          max_iter = 40, tol = 0.01,
          norm = "iqr",
          cap_frac = 0.03,
          record = FALSE
        )

        fleets[[l]]$e_p_s[, s] <- res$eff_p

        } else {
        stop("spatial effort allocation strategy not properly defined, check spatial_allocation and cost_per_unit_effort in fleet object")
      }
    } # close loop over fleets

    for (f in seq_along(fauni)) {

      critter <- fauni[f]

      ages <- length(fauna[[critter]]$length_at_age)

      last_n_p_a <- storage[[s - 1]][[critter]]$n_p_a

      f_p_a <-
        matrix(0, nrow = patches, ncol = ages) # total fishing mortality by patch and age

      f_p_a_fl <-
        array(
          0,
          dim = c(patches, ages, length(fleets)),
          dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
        ) # storage for proportion of fishing mortality by patch, age, and fleet

      p_p_a_fl <-
        array(
          0,
          dim = c(patches, ages, length(fleets)),
          dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
        ) # storage for price by patch, age, and fleet

      for (l in seq_along(fleet_names)) {
        tmp <-
          matrix(
            fleets[[l]]$metiers[[critter]]$spatial_catchability,
            nrow = nrow(last_n_p_a),
            ncol = ncol(last_n_p_a),
            byrow = FALSE
          ) *
          matrix(
            fleets[[l]]$metiers[[critter]]$sel_at_age,
            nrow = nrow(last_n_p_a),
            ncol = ncol(last_n_p_a),
            byrow = TRUE
          )
        ## could add in the effective discard factor here, where that would be a multipliier as a function of 1 - (discard_rate * discard_survival)

        f_p_a <-
          f_p_a + fleets[[l]]$e_p_s[, s] * tmp

        f_p_a_fl[, , l] <-
          fleets[[l]]$e_p_s[, s] * tmp

        p_p_a_fl[, , l] <- fleets[[l]]$metiers[[critter]]$price
      } # calculate cumulative f at age by patch

      f_p_a_fl <-
        f_p_a_fl / array(
          f_p_a,
          dim = c(patches, ages, length(fleets)),
          dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
        ) # f by patch, age, and fleet

      movement <- fauna[[critter]]$movement_matrix

      # if there is updated habitat for the critter in question in current time step, update habitat
      if ((length(habitat) > 0) &
          (critter %in% names(habitat))) {
        season_block <-
          which(sapply(fauna[[critter]]$movement_seasons, function(x, y) {
            any(y %in% x)
          }, x = current_season)) # figure out which season block you are in

        # update habitat in this time step
        current_habitat <- habitat[[critter]][[s - 1]]

        # need to use pivot_longer to match patch order from expand_grid
        # SPEED (spirit-of-comment): we precomputed the equivalent patch-ordered vector with as.vector()
        hab_vals <- habitat_vecs[[critter]][[s - 1]]

        neighbors <- find_neighbors(resolution = resolution, water = is.finite(hab_vals))

        edges <- Matrix::summary(neighbors)

        to <- edges$i

        from   <- edges$j

        delta_h <- hab_vals[to] - hab_vals[from]

        mult <- exp((time_step * delta_h) / sqrt(patch_area))

        mult <- pmin(mult, fauna[[critter]]$max_hab_mult)

        P <- length(hab_vals)

        current_habitat <-  Matrix::sparseMatrix(
          i = to,
          j = from,
          x = mult,
          dims = c(P, P),
          dimnames = dimnames(neighbors)
        )


        diffusion_and_taxis <-
          fauna[[critter]]$diffusion_foundation[[season_block]] * current_habitat

        inst_movement_matrix <-
          prep_movement(diffusion_and_taxis, resolution = resolution)

        movement[[season_block]] <- sparsify_transition(as.matrix(expm::expm(as.matrix(inst_movement_matrix))))

        if (any(!is.finite(movement[[season_block]]))) {
          stop("scale of supplied habitat differences are too extreme, try rescaling so that the exponent of the differences are less extreme in magnitude")
        }
      }

      if (!(current_season %in% fishing_seasons[[critter]])) {
        f_p_a <- f_p_a * 0
      }

      fauna_rec_devs <- rep(exp(log_rec_devs[s, critter] - fauna[[critter]]$sigma_rec^2 / 2), patches)

      pop <- fauna[[critter]]$swim(
        season = current_season,
        adult_movement = movement,
        f_p_a = f_p_a,
        last_n_p_a = last_n_p_a,
        rec_devs = fauna_rec_devs
      )

      c_p_a_fl <- f_p_a_fl * array(
        pop$c_p_a,
        dim = c(patches, ages, length(fleets)),
        dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
      )

      fmult <- 1
      if (length(manager$quotas[critter]) > 0) {
        if (manager$quotas[[critter]] < sum(c_p_a_fl, na.rm = TRUE)) {
          quota <- manager$quotas[[critter]]

          fmulter <- optim(
            par = 0.9,
            fn = marlin::quota_finder,
            quota = quota,
            fauna = fauna,
            current_season = current_season,
            movement = movement,
            f_p_a = f_p_a,
            last_n_p_a = last_n_p_a,
            f_p_a_fl = f_p_a_fl,
            critter = critter,
            patches = patches,
            ages = ages,
            fleets = fleets,
            rec_devs = fauna_rec_devs,
            lower = 0,
            upper = 1,
            method = "L-BFGS-B"
          )

          fmult <- fmulter$par
          f_p_a <- f_p_a * fmult

          pop <- fauna[[critter]]$swim(
            season = current_season,
            adult_movement = movement,
            f_p_a = f_p_a,
            last_n_p_a = last_n_p_a,
            rec_devs = fauna_rec_devs
          )

          c_p_a_fl <- f_p_a_fl * array(
            pop$c_p_a,
            dim = c(patches, ages, length(fleets)),
            dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
          )
        }
      }

      r_p_a_fl <- c_p_a_fl * p_p_a_fl

      tmp_e_p_fl <-
        purrr::list_cbind(unname(purrr::map(
          fleets, ~ data.frame(x = as.numeric(.x$e_p_s[, s] * fmult))
        )), name_repair = "unique_quiet")
      colnames(tmp_e_p_fl) <- names(fleets)

      for (fl in 1:length(fleets)) {
        c_p_fl[, fl] <- rowSums(c_p_a_fl[, , fl], na.rm = TRUE)

        r_p_fl[, fl] <- rowSums(r_p_a_fl[, , fl], na.rm = TRUE)

        prof_p_fl[, fl] <-
          r_p_fl[, fl] - fleets[[fl]]$cost_per_unit_effort * ((
            as.matrix(tmp_e_p_fl[, fl]) / length(fauna)^fleets[[fl]]$effort_cost_exponent
          ) + as.matrix(tmp_e_p_fl[, fl] / length(fauna) * fleets[[fl]]$cost_per_patch)
          )
      }

      storage[[s - 1]][[critter]]$c_p_fl <-
        c_p_fl # store catch by patch  by fleet

      storage[[s - 1]][[critter]]$r_p_fl <-
        r_p_fl # store revenue by patch  by fleet

      storage[[s - 1]][[critter]]$prof_p_fl <-
        prof_p_fl # store profits by patch by fleet

      storage[[s - 1]][[critter]]$c_p_a_fl <-
        c_p_a_fl # catch stored in each model is the catch that came from the last time step, so put in the right place here

      storage[[s - 1]][[critter]]$r_p_a_fl <-
        r_p_a_fl # revenue stored in each model is the revenue that came from the last time step, so put in the right place here

      storage[[s - 1]][[critter]]$c_p_a <-
        pop$c_p_a # catch stored in each model is the catch that came from the last time step, so put in the right place here

      storage[[s - 1]][[critter]]$e_p_fl <-
        tmp_e_p_fl # store effort by patch by fleet (note that this is the same across species)

      storage[[s - 1]][[critter]]$f_p_a_fl <-
        f_p_a_fl # store effort by patch by fleet (note that this is the same across species)

      if (any(tmp_e_p_fl < 0)) {
        stop("something hase gone very wrong, effort is negative")
      }

      storage[[s]][[critter]] <- pop

      storage[[s]][[critter]]$b0 <- fauna[[critter]]$b0
    } # close fauni, much faster this way than dopar, who knew
  } # close steps

  trimmed_names <- step_names[ifelse(keep_starting_step, 1, 2):steps]

  trimmed_names <-
    trimmed_names[1:(length(trimmed_names) - 1)] # a seriously annoying edge case in case you are only running one step

  storage <-
    storage[ifelse(keep_starting_step, 1, 2):pmax(2, steps - 1)] # since catch is retrospective, chop off last time step to ensure that every step has a catch history, and drop starting step is specified

  storage <-
    rlang::set_names(storage, nm = paste0(trimmed_names))

  storage <- purrr::map(storage, ~ rlang::set_names(.x, fauni))

  return(storage)
} # close function

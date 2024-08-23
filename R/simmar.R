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

  n_critters <- length(fauna)

  fleet_names <- names(fleets)

  time_step <- unique(purrr::map_dbl(fauna, "time_step"))

  steps_per_year <- 1 / time_step

  patch_area <- unique(purrr::map_dbl(fauna, "patch_area"))

  sigma_recs <- purrr::map_dbl(fauna, "sigma_rec") # gather recruitment standard deviations

  ac_recs <- purrr::map_dbl(fauna, "ac_rec") # gather autocorrelation in recruitment standard deviations

  if (any(abs(ac_recs) >= 1)){
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
      (years) / time_step + 1  # tack on extra step for accounting
  } else {
    steps <- steps + 1 # to store initial conditions
  } # if steps are specified instead of years

  steps <- pmax(steps, 3) # need to be at least 1 initial + 2 running steps

  # generate recruitment deviates


  if (!is.na(starting_step)) {
    year_season <- marlin::clean_steps(starting_step)

    starting_year <-
      as.integer(gsub("_.*$", "", year_season)) - 1 # -1 to account for years being 1 indexed

    starting_season <-
      as.integer(gsub("^.*_", "", year_season)) - 1 # -1 to account for it would be starting in the third season

    offset <- (starting_year * steps_per_year) + starting_season

    # tack on a number of steps equal to the number of

  } else {
    offset <- 0
  }

  step_names <-
    paste(rep(1:(steps + offset), each = steps_per_year), 1:steps_per_year, sep = "_") # generate at least as many steps as you're going to need

  step_names <-
    step_names[1:steps + offset] # chop back to the actual number of steps used; works this way in case there are multiple steps per year but an incomplete number of years

  if (is.null(log_rec_devs)) {
  log_rec_devs <- matrix(NA, nrow = length(step_names), ncol = n_critters, dimnames = list(step_names, fauni))

  log_rec_devs[1,] <- mvtnorm::rmvnorm(1,rep(0, n_critters),sigma = covariance_rec)

  for (i in 2:length(step_names)){

    log_rec_devs[i, ] <- ac_recs *  log_rec_devs[i - 1, ] + sqrt(1 - ac_recs ^ 2) * mvtnorm::rmvnorm(1, rep(0,n_critters), sigma = covariance_rec)

  }


}
  # step_names <-  paste(rep((1:steps) + offset, each = steps_per_year), 1:steps_per_year, sep = "_") # generate at least as many steps as you're going to need

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
              as.integer(gsub("_.*$", "", step_names[i])) -  as.integer(gsub("_.*$", "", starting_step)) + 1
            # put year in index form not named form

          } else {
            years_in <-
              as.integer(gsub("_.*$", "", step_names[i]))
          }
          years_in <- min(years_in, supplied_years)


          # year <-  floor(step_names[i] - starting_step) # put year in index form not named form
          new_habitat[[i]] <- habitat[[f]][[years_in]]
        }

        habitat[[f]] <- new_habitat

      } # close expansion of habitat if needed

    } # close fauni loop

  } # close check on supplied habitat


  storage <- vector("list", steps)

  storage[[1]] <-
    initial_conditions # start populations at initial conditions

  ## XX need to modify this to allow for fishable patches
  fleets <-
    purrr::map(fleets, ~ purrr::list_modify(.x, e_p_s = matrix(((.x$base_effort / patches)
    ), nrow = patches, ncol = steps))) # create blank for effort by fleet, space, and time

  if (init_cond_provided) {
    # fill in first step of effort from initial conditions
    for (i in names(fleets)) {
      fleets[[i]]$e_p_s[, 1] <-
        getElement(initial_conditions[[1]]$e_p_fl, i)

    }

  }

  r_p_f <-
    matrix(0, patches, length(fauni)) # fishable revenue by patch and fauna

  e_p_f <-
    matrix(0, patches, length(fauni)) # total fishing effort by patch and fauna

  f_q <- rep(0, length(fauni)) # storage for q by fauna

  c_p_fl <-
    matrix(
      nrow = patches,
      ncol = length(fleets),
      dimnames = list(NULL, fleet_names)
    )

  r_p_fl <-
    matrix(
      nrow = patches,
      ncol = length(fleets),
      dimnames = list(NULL, fleet_names)
    )

  prof_p_fl <-
    matrix(
      nrow = patches,
      ncol = length(fleets),
      dimnames = list(NULL, fleet_names)
    )


  fishable <- rep(1, patches)

  # loop over steps
  for (s in 2:steps) {
    last_season <- as.integer(gsub("^.*_", "", step_names[s - 1]))

    year <-  as.integer(gsub("_.*$", "", step_names[s]))

    current_season <- as.integer(gsub("^.*_", "", step_names[s]))

    if (length(manager$mpas) > 0) {
      manager$mpas$locations <- manager$mpas$locations |>
        dplyr::arrange(x, y)

      # assign MPAs if needed
      if (year == manager$mpas$mpa_year) {
        fishable <- manager$mpas$locations$mpa == 0
      }

    } # close MPA if statement


    for (l in seq_along(fleet_names)) {
      # distribute fleets in space based on revenues

      fleet_fishable <- fishable # keep base fishable, but modify by fleet-specific fishing grounds as needed

      if (!is.null(fleets[[l]]$fishing_grounds)) {
        fleets[[l]]$fishing_grounds <- fleets[[l]]$fishing_grounds |>
          dplyr::arrange(x, y) # make sure patches are in correct order

        fleet_fishable[fleets[[l]]$fishing_grounds$fishing_ground == FALSE] <- 0


      }

      concentrator <-
        rep(1, patches) # reset fishing effort concentrator by fleet

      if (length(manager$mpas) > 0) {
        if (year >= manager$mpas$mpa_year &
            fleets[[l]]$mpa_response == "leave") {
          concentrator <- as.numeric(fishable)
        }
      }

      # xx add ability to incorporate past revenues here. Idea. have a marker that lets you know if initial conditions were passed in. If s<= 2 or there were no initial conditions, pull this. If s<= 2 but there were initial conditions, pull the initial conditions for those steps

      if (s == 2 & init_cond_provided) {
        r_p_f <-
          (sapply(initial_conditions, function(x)
            rowSums(x$r_p_a_fl[, , l], na.rm = TRUE)))

        last_r_p <-
          rowSums(r_p_f, na.rm = TRUE) # pull out total revenue for fleet l
      } else if (s <= 2) {
        for (f in seq_along(fauni)) {
          last_b_p_a <- storage[[s - 1]][[f]]$b_p_a

          last_e_p <- fleets[[l]]$e_p_s[, s - 1]

          # calculate fishable biomass in each patch for each species for that fleet

          # account for spatial catchability
          tmp = 1 - exp(-(
            matrix(
              fleets[[l]]$metiers[[fauni[f]]]$spatial_catchability,
              nrow = nrow(last_b_p_a),
              ncol = ncol(last_b_p_a),
              byrow = FALSE
            ) *
              matrix(
                fleets[[l]]$metiers[[fauni[f]]]$sel_at_age,
                nrow = nrow(last_b_p_a),
                ncol = ncol(last_b_p_a),
                byrow = TRUE
              )
          ))

          last_b_p <-
            rowSums(last_b_p_a * tmp * (current_season %in% fishing_seasons[[f]])) * fleet_fishable

          r_p_f[, f] <-
            last_b_p * fleets[[l]]$metiers[[fauni[f]]]$price

          f_q[f] <- fleets[[l]]$metiers[[fauni[f]]]$catchability

        } # close fauni loop

        last_r_p <- rowSums(r_p_f, na.rm = TRUE)

      } else {
        r_p_f <-
          (sapply(storage[[s - 2]], function(x)
            rowSums(x$r_p_a_fl[, , l], na.rm = TRUE)))

        last_r_p <-
          rowSums(r_p_f, na.rm = TRUE) # pull out total revenue for fleet l


      }

      total_effort <- sum(fleets[[l]]$e_p_s[, s - 1] * concentrator)

      if (sum(last_r_p, na.rm = TRUE) > 0) {
        # the only way for revenues in the last step to be literally zero is 100% mpas or all target species seasons closed or last effort = zero
        last_revenue <-
          sum(last_r_p, na.rm = TRUE) # pull out total revenue for fleet l

        last_cost <-
          fleets[[l]]$cost_per_unit_effort * (
            sum((fleets[[l]]$e_p_s[, s - 1]) ^
                  fleets[[l]]$effort_cost_exponent) + sum(fleets[[l]]$cost_per_patch * fleets[[l]]$e_p_s[, s -
                                                                                                           1])
          )

        last_profits <-
          last_revenue - last_cost # calculate profits in the last time step.

        last_r_to_c <- last_revenue / last_cost
      }

      if (fleets[[l]]$fleet_model == "open_access") {
        if (is.na(fleets[[l]]$cost_per_unit_effort) |
            is.na(fleets[[l]]$responsiveness)) {
          stop(
            "open access fleet model requires both cost_per_unit_effort and responsiveness parameters"
          )
        }


        if (exists("last_revenue")) {
          effort_cap <- Inf

          if (length(manager$effort_cap[[l]]) > 0) {
            effort_cap <- manager$effort_cap[[l]]

          }

          total_effort <-
            pmin(effort_cap, total_effort * pmin(1.5, exp(
              fleets[[l]]$responsiveness * log(pmax(last_revenue, 1e-6) / pmax(1e-6, last_cost))
            ))) # adjust effort per an open access dynamics model
        } # in edge case where the fishery is closed for the first few seasons of the simulation stick with last value

      } # close open access

      e_p <- fleets[[l]]$e_p_s[, s - 1]

      # allocate fleets in space

      if (fleets[[l]]$spatial_allocation == "revenue") {
        if (sum(fleet_fishable) == 0) {
          alloc <- 0
        } else if (sum(last_r_p, na.rm = TRUE) == 0) {
          # if there is no revenue anywhere just distribute fleet evenly as an edge case for extreme overfishing

          alloc <- fleet_fishable / sum(fleet_fishable)

          #1 / nrow(r_p_f)
        } else {
          alloc <-
            ((last_r_p * fleet_fishable) / sum(last_r_p * fleet_fishable, na.rm = TRUE)) # just extra cautios.

          alloc <-
            alloc - min(alloc, na.rm = TRUE) + 1 # rescale to be greater than or equal to 1

          alloc <- alloc / sum(alloc)

          if (s > 2) {
            alloc <-
              rowSums(sapply(storage[[s - 2]], function(x)
                x$r_p_fl[, l]), na.rm = TRUE)

            alloc <-
              alloc - min(alloc, na.rm = TRUE) + 1 # rescale to be greater than or equal to 1

            alloc <- alloc * fleet_fishable

            alloc <- alloc / sum(alloc)

          }


        }

        fleets[[l]]$e_p_s[, s] <-
          total_effort * alloc # distribute fishing effort by fishable biomass


      } else if (fleets[[l]]$spatial_allocation == "rpue") {
        if (sum(fleet_fishable) == 0) {
          alloc <- 0
        } else if (sum(last_r_p, na.rm = TRUE) == 0) {
          # if there is no revenue anywhere just distribute fleet evenly as an edge case for extreme overfishing

          alloc <- fleet_fishable / sum(fleet_fishable)

          #1 / nrow(r_p_f)
        } else {
          alloc = (last_r_p / (e_p)) * fleet_fishable

          alloc[!is.finite(alloc)] <-  0

          alloc <-
            alloc - min(alloc, na.rm = TRUE) + 1 # rescale to be greater than or equal to 1, in case prices are negative

          alloc <-  alloc / sum(alloc, na.rm = TRUE)

          if (s > 2) {
            alloc <-
              rowSums(sapply(storage[[s - 2]], function(x)
                x$r_p_fl[, l]), na.rm = TRUE) / (e_p)

            alloc[!is.finite(alloc)] <-  0

            alloc <-
              alloc - min(alloc, na.rm = TRUE) + 1 # rescale to be greater than or equal to 1

            alloc <- alloc * fleet_fishable

            alloc <- alloc / sum(alloc)

            # if (s == (steps-10)  & fleet_names[l] != "artisanal"){
            #   browser()
            # }

          }

        }

        fleets[[l]]$e_p_s[, s] <-
          total_effort * alloc # distribute fishing effort by fishable biomass


      } else if (fleets[[l]]$spatial_allocation == "ppue" &&
                 !is.na(fleets[[l]]$cost_per_unit_effort)) {
        if (sum(fleet_fishable) == 0) {
          alloc <- 0
        } else if (sum(last_r_p, na.rm = TRUE) == 0) {
          # if there is no revenue anywhere just distribute fleet evenly as an edge case for extreme overfishing

          alloc <- fleet_fishable / sum(fleet_fishable)

        } else {
          alloc = ((
            last_r_p - fleets[[l]]$cost_per_unit_effort * ((e_p) ^ fleets[[l]]$effort_cost_exponent
                                                           + fleets[[l]]$cost_per_patch * e_p
            ) / (e_p + 1)
          )) * fleet_fishable

          alloc[!is.finite(alloc)] <-  0

          alloc <-
            alloc - min(alloc, na.rm = TRUE) + 1 # rescale to be greater than or equal to 1

          alloc <-  alloc / sum(alloc, na.rm = TRUE)

          if (s > 2) {
            alloc <-
              rowSums(sapply(storage[[s - 2]], function(x)
                x$prof_p_fl[, l]), na.rm = TRUE) / (e_p)

            alloc[!is.finite(alloc)] <-  0

            alloc <-
              fleet_fishable * (alloc - min(alloc, na.rm = TRUE) + 1) # rescale to be greater than or equal to 1

            alloc <- alloc / sum(alloc)


          }

        }

        fleets[[l]]$e_p_s[, s] <-
          total_effort * alloc # distribute fishing effort by fishable biomass



      } else if (fleets[[l]]$spatial_allocation == "profit" &&
                 !is.na(fleets[[l]]$cost_per_unit_effort)) {
        if (sum(fleet_fishable) == 0) {
          alloc <- 0
        } else if (sum(last_r_p, na.rm = TRUE) == 0) {
          # if there is no revenue anywhere just distribute fleet evenly as an edge case for extreme overfishing

          alloc <- fleet_fishable / sum(fleet_fishable)

          #1 / nrow(r_p_f)
        } else {
          alloc = (
            last_r_p - fleets[[l]]$cost_per_unit_effort * ((e_p) ^ fleets[[l]]$effort_cost_exponent +  fleets[[l]]$cost_per_patch * e_p
            )
          )

          alloc[!is.finite(alloc)] <-  0

          alloc <-
            fleet_fishable * (alloc - min(alloc, na.rm = TRUE) + 1) # rescale to be greater than or equal to 1

          alloc <-  alloc / sum(alloc, na.rm = TRUE)

          if (s > 2) {
            alloc <-
              rowSums(sapply(storage[[s - 2]], function(x)
                x$prof_p_fl[, l]), na.rm = TRUE)

            alloc <-
              fleet_fishable * (alloc - min(alloc, na.rm = TRUE) + 1) # rescale to be greater than or equal to 1

            alloc <- alloc / sum(alloc)

          }


        }

        fleets[[l]]$e_p_s[, s] <-
          total_effort * alloc # distribute fishing effort



      } else if (fleets[[l]]$spatial_allocation == "ideal_free" &&
                 !is.na(fleets[[l]]$cost_per_unit_effort)) {
        stop(
          "ideal free distribution not yet supported. Set spatial_allocation = 'revenue' in create_fleet"
        )

        # calculate expected marginal revenue when effort = 0 in each patch
        # marginal revenue is fishable revenue (r_p_f) - marginal cost per unit effort

        worth_fishing <-
          ((last_r_p - fleets[[l]]$cost_per_unit_effort) * fleet_fishable) > 0 # check whether any effort could be positive, factoring in potential for negative price,

        # for patches that will support any fishing, solve for effort such that marginal profits are equal in space

        # optfoo <- function(log_effort, mp = 0, revs, qs, cost) {
        #   mp_hat <- sum(revs * exp(-qs * exp(log_effort))) - cost
        #
        #   ss <- (mp_hat - mp) ^ 2
        # }

        id_e_p <- rep(0, patches)

        fishable_patches <-
          (1:patches)[(fleet_fishable == 1) & worth_fishing]

        r_p_f <-
          matrix(r_p_f[fishable_patches, ], nrow  = length(fishable_patches)) # seriously annoying step to preserve matrix structure when there is only one species
        tmp <-
          matrix(rep(f_q, nrow(r_p_f)),
                 nrow = nrow(r_p_f),
                 byrow = TRUE)

        id_e_p[fishable_patches] <-
          ((rowSums(log(r_p_f * tmp))) - log(fleets[[l]]$cost_per_unit_effort)) / sum(tmp) # see physical paper TNC notebook for derivation of this, asumes that profits = p * b * exp(-(q * E)) - c


        # for (pp in seq_along(fishable_patches)) {
        #   id_e_p[fishable_patches[pp]] <-
        #     exp(
        #       nlminb(
        #         log(total_effort / length(fishable_patches + 1e-3)),
        #         optfoo,
        #         revs = r_p_f[fishable_patches[pp], ],
        #         qs = f_q,
        #         cost = fleets[[l]]$cost_per_unit_effort
        #       )$par
        #     )
        #
        # }

        # allocate available effort

        choices <-
          dplyr::arrange(data.frame(patch = 1:patches, effort = id_e_p),
                         desc(effort))

        choices$cumu_effort <- cumsum(choices$effort)

        choices <- choices[choices$cumu_effort < total_effort, ]

        choices$alloc_effort <-
          total_effort * (choices$effort / sum(choices$effort))

        fleets[[l]]$e_p_s[, s] <- 0

        fleets[[l]]$e_p_s[choices$patch, s] <- choices$alloc_effort


      } # close ifd else
      else if ((fleets[[l]]$spatial_allocation == "manual")) {
        alloc <- fleet_fishable * fleets[[l]]$fishing_grounds$fishing_ground

        alloc <- alloc / sum(alloc)

        fleets[[l]]$e_p_s[, s] <-
          total_effort * alloc # distribute fishing effort by fishable biomass

      }
      else {
        stop(
          "spatial effort allocation strategy not properly defined, check spatial_allocation and cost_per_unit_effort in fleet object"
        )
      }

    } # close loop over fleets

    for (f in seq_along(fauni)) {
      # run population model for each species

      ages <-  length(fauna[[f]]$length_at_age)

      last_n_p_a <-
        storage[[s - 1]][[f]]$n_p_a # numbers by patch and age in the last time step

      f_p_a <-
        matrix(0, nrow = patches, ncol = ages) # total fishing mortality by patch and age

      f_p_a_fl <-
        array(
          0,
          dim = c(patches, ages, length(fleets)),
          dimnames = list(1:patches, fauna[[f]]$ages, names(fleets))
        ) # storage for proportion of fishing mortality by patch, age, and fleet

      p_p_a_fl <-
        array(
          0,
          dim = c(patches, ages, length(fleets)),
          dimnames = list(1:patches, fauna[[f]]$ages, names(fleets))
        ) # storage for price by patch, age, and fleet

      for (l in seq_along(fleet_names)) {
        tmp =
          matrix(
            fleets[[l]]$metiers[[fauni[f]]]$spatial_catchability,
            nrow = nrow(last_n_p_a),
            ncol = ncol(last_n_p_a),
            byrow = FALSE
          ) *
          matrix(
            fleets[[l]]$metiers[[fauni[f]]]$sel_at_age,
            nrow = nrow(last_n_p_a),
            ncol = ncol(last_n_p_a),
            byrow = TRUE
          )
        ## could add in the effective discard factor here, where that would be a multipliier as a function of 1 - (discard_rate * discard_survival)

        f_p_a <-
          f_p_a + fleets[[l]]$e_p_s[, s] * tmp

        f_p_a_fl[, , l] <-
          fleets[[l]]$e_p_s[, s] * tmp

        p_p_a_fl[, , l] <- fleets[[l]]$metiers[[fauni[f]]]$price

      } # calculate cumulative f at age by patch

      f_p_a_fl <-
        f_p_a_fl / array(
          f_p_a,
          dim = c(patches, ages, length(fleets)),
          dimnames = list(1:patches, fauna[[f]]$ages, names(fleets))
        ) # f by patch, age, and fleet




      movement <- fauna[[f]]$movement_matrix


      # if there is updated habitat for the critter in question in current time step, update habitat
      if ((length(habitat) > 0) &
          (names(fauna)[[f]] %in% names(habitat))) {
        season_block <-
          which(sapply(fauna[[f]]$movement_seasons, function(x, y)
            any(y %in% x), x = current_season)) # figure out which season block you are in

        # update habitat in this time step

        current_habitat <- habitat[[names(fauna)[[f]]]][[s - 1]]

        current_habitat <-
          tidyr::pivot_longer(as.data.frame(current_habitat), tidyr::everything()) # need to use pivot_longer to match patch order from expand_grid

        current_habitat <-
          pmin(exp((
            time_step * outer(current_habitat$value, current_habitat$value, "-")
          ) / sqrt(patch_area)), fauna[[f]]$max_hab_mult) # convert habitat gradient into diffusion multiplier


        diffusion_and_taxis <-
          fauna[[f]]$diffusion_foundation[[season_block]] * current_habitat

        inst_movement_matrix <-
          prep_movement(diffusion_and_taxis, resolution = resolution)


        # update movement matrix with current habitat
        movement[[season_block]] <-
          as.matrix(expm::expm(inst_movement_matrix))


        if (any(!is.finite(movement[[season_block]]))) {
          stop(
            "scale of supplied habitat differences are too extreme, try rescaling so that the exponent of the differences are less extreme in magnitude"
          )
        }

      } # close habitat update




      if (!(current_season %in% fishing_seasons[[names(fauna)[f]]])) {
        f_p_a <- f_p_a * 0 # turn off fishing if the season is closed

      }

      fauna_rec_devs <- rep(exp(log_rec_devs[s,f] - fauna[[f]]$sigma_rec ^ 2 / 2), patches)

      pop <-
        fauna[[f]]$swim(
          season = current_season,
          adult_movement = movement,
          f_p_a = f_p_a,
          last_n_p_a = last_n_p_a,
          rec_devs = fauna_rec_devs
        )
      # process catch data
      c_p_a_fl <-
        f_p_a_fl * array(
          pop$c_p_a,
          dim = c(patches, ages, length(fleets)),
          dimnames = list(1:patches, fauna[[f]]$ages, names(fleets))
        )

      # if there are any quotas to evaluate
      fmult <- 1

      if (length(manager$quotas[names(fauna)[f]]) > 0) {
        if (manager$quotas[[names(fauna)[f]]] < sum(c_p_a_fl, na.rm = TRUE)) {
          quota <- manager$quotas[[names(fauna)[f]]]

          fmulter <-
            nlminb(
              0.9,
              marlin::quota_finder,
              quota = quota,
              fauna = fauna,
              current_season = current_season,
              movement = movement,
              f_p_a = f_p_a,
              last_n_p_a = last_n_p_a,
              f_p_a_fl = f_p_a_fl,
              f = f,
              patches = patches,
              ages = ages,
              fleets = fleets,
              rec_devs = rec_devs,
              lower = 0,
              upper = 1
            )

          fmult <- fmulter$par

          f_p_a <- f_p_a * fmult

          pop <-
            fauna[[f]]$swim(
              season = current_season,
              adult_movement = movement,
              f_p_a = f_p_a,
              last_n_p_a = last_n_p_a,
              rec_devs = rec_devs
            )
          # process catch data
          c_p_a_fl <-
            f_p_a_fl * array(
              pop$c_p_a,
              dim = c(patches, ages, length(fleets)),
              dimnames = list(1:patches, fauna[[f]]$ages, names(fleets))
            )

        } # if the quota is less than the catch, enforce the quota


      } # close quota evaluation

      r_p_a_fl <- c_p_a_fl * p_p_a_fl

      tmp_e_p_fl <-
        purrr::list_cbind(unname(purrr::map(
          fleets, ~ data.frame(x = as.numeric(.x$e_p_s[, s] * fmult))
        )), name_repair = "unique_quiet")


      colnames(tmp_e_p_fl) <-  names(fleets)


      for (fl in 1:length(fleets)) {
        c_p_fl[, fl] <- rowSums(c_p_a_fl[, , fl], na.rm = TRUE)

        r_p_fl[, fl] <-  rowSums(r_p_a_fl[, , fl], na.rm = TRUE)

        prof_p_fl[, fl] <-
          r_p_fl[, fl] - fleets[[fl]]$cost_per_unit_effort * ((
            as.matrix(tmp_e_p_fl[, fl]) / length(fauna) ^ fleets[[fl]]$effort_cost_exponent
          )  + as.matrix(tmp_e_p_fl[, fl] / length(fauna) * fleets[[fl]]$cost_per_patch)
          )

      }

      storage[[s - 1]][[f]]$c_p_fl <-
        c_p_fl # store catch by patch  by fleet

      storage[[s - 1]][[f]]$r_p_fl <-
        r_p_fl # store revenue by patch  by fleet

      storage[[s - 1]][[f]]$prof_p_fl <-
        prof_p_fl # store profits by patch by fleet

      storage[[s - 1]][[f]]$c_p_a_fl <-
        c_p_a_fl # catch stored in each model is the catch that came from the last time step, so put in the right place here

      storage[[s - 1]][[f]]$r_p_a_fl <-
        r_p_a_fl # revenue stored in each model is the revenue that came from the last time step, so put in the right place here

      storage[[s - 1]][[f]]$c_p_a <-
        pop$c_p_a # catch stored in each model is the catch that came from the last time step, so put in the right place here

      storage[[s - 1]][[f]]$e_p_fl <-
        tmp_e_p_fl # store effort by patch by fleet (note that this is the same across species)


      storage[[s - 1]][[f]]$f_p_a_fl <-
        f_p_a_fl # store effort by patch by fleet (note that this is the same across species)

      if (any(tmp_e_p_fl < 0)) {
        stop("something hase gone very wrong, effort is negative")
      }

      storage[[s]][[f]] <- pop
    } # close fauni, much faster this way than dopar, who knew



  } #close steps
  trimmed_names <- step_names[ifelse(keep_starting_step, 1, 2):steps]

  trimmed_names <-
    trimmed_names[1:(length(trimmed_names) - 1)] # a seriously annoying edge case in case you are only running one step

  storage <-
    storage[ifelse(keep_starting_step, 1, 2):pmax(2, steps - 1)] # since catch is retrospective, chop off last time step to ensure that every step has a catch history, and drop starting step is specified
  storage <-
    rlang::set_names(storage, nm = paste0(trimmed_names))

  storage <- purrr::map(storage, ~ rlang::set_names(.x, fauni))
} # close function

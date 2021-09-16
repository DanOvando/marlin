#' simmar is the wrapper function for the marlin package
#'
#' when passed fauna and fleet objects, simmar will advance
#' the population for a number of steps
#'
#' @param fauna
#' @param fleets
#' @param mpas
#' @param steps
#' @param tune_unfished
#'
#' @return
#' @export
#'
simmar <- function(fauna = list(),
                   fleets = list(),
                   mpas = list(),
                   years = 100,
                   tune_unfished = 0,
                   initial_conditions = NA) {
  fauni <- names(fauna)
  
  fleet_names <- names(fleets)
  
  time_step <- unique(purrr::map_dbl(fauna, "time_step"))
  
  if (length(time_step) > 1) {
    stop(
      paste(
        "All critters in fauna must have the same time step: current time steps are",
        paste(time_step, collapse = " ")
      )
    )
  }
  
  steps <-
    (years) / time_step  + 1 # tack on extra year for accounting
  
  
  patches <- unique(purrr::map_dbl(fauna, "patches"))
  if (all(is.na(initial_conditions))){
    
    initial_conditions <-
      purrr::map(fauna, c("unfished")) # pull out unfished conditions created by create_critter
    
  }

  if (length(patches) > 1) {
    stop(
      "fauna have different habitat resolutions: set resolution to same number for all species!"
    )
  }
  
  storage <- vector("list", steps)
  
  storage[[1]] <-
    initial_conditions # start populations at initial conditions
  
  fleets <-
    purrr::map(fleets, ~ purrr::list_modify(.x,
                                            e_p_s = matrix(((.x$base_effort / patches)
                                            ), nrow = patches, ncol = steps))) # create blank for effort by fleet, space, and time
  
  r_p_f <-
    matrix(0, patches, length(fauni)) # fishable revenue by patch and fauna
  
  e_p_f <-
    matrix(0, patches, length(fauni)) # total fishing mortality by patch and fauna
  
  f_q <- rep(0, length(fauni)) # storage for q by fauna
  
  fishable <- rep(1, patches)
  
  # step_seq <-
  #   seq(0, steps, by = time_step) # create sequence of seasonal time steps
  #
  step_names <- seq(0, years + 1, by = time_step)
  
  # browser()
  
  # loop over steps
  for (s in 2:steps) {
    # season <-
    #   step_seq[s - 1] - floor(step_seq[s - 1]) # determine what season the last time step was
    #
    season <-
      step_names[s - 1] - floor(step_names[s - 1]) # determine what season the last time step was
    
    year <-  floor(step_names[s])
    
    if (length(mpas) > 0) {
      # assign MPAs if needed
      if (year == mpas$mpa_year) {
        fishable <- mpas$locations$mpa == 0
      }
      
    } # close MPA if statement
    
    
    for (l in seq_along(fleet_names)) {
      # distribute fleets in space based on revenues
      
      concentrator <-
        rep(1, patches) # reset fishing effort concentrator by fleet
      
      if (length(mpas) > 0) {
        if (year >= mpas$mpa_year & fleets[[l]]$mpa_response == "leave") {
          concentrator <- as.numeric(fishable)
        }
      }
      
      
      if (s <= 2){
      for (f in seq_along(fauni)) {
        last_b_p_a <- storage[[s - 1]][[f]]$b_p_a
        
        last_e_p <- fleets[[l]]$e_p_s[, s - 1]

        # calculate fishable biomass in each patch for each species for that fleet
        
        # browser()
        # last_revenue <-  sum((sapply(storage[[s-2]], function(x) sum(x$r_p_a_fl[,,l], na.rm = TRUE))), na.rm = TRUE) # pull out total revenue for fleet l
        # 
        # account for spatial catchability
        tmp = 1 - exp(-(
          matrix(
            fleets[[l]]$metiers[[fauni[f]]]$spatial_catchability,
            nrow = nrow(last_b_p_a),
            ncol = ncol(last_b_p_a),
            byrow = TRUE
          ) *
            matrix(
              fleets[[l]]$metiers[[fauni[f]]]$sel_at_age,
              nrow = nrow(last_b_p_a),
              ncol = ncol(last_b_p_a),
              byrow = TRUE
            )
        ))
        
        last_b_p <- rowSums(last_b_p_a * tmp) * fishable
        
        r_p_f[, f] <-
          last_b_p * fleets[[l]]$metiers[[fauni[f]]]$price
        
        f_q[f] <- fleets[[l]]$metiers[[fauni[f]]]$catchability
        
      } # close fauni loop
        
        last_r_p <- rowSums(r_p_f, na.rm = TRUE)
      
      } else {
        
        last_r_p <-  rowSums((sapply(storage[[s-2]], function(x) rowSums(x$r_p_a_fl[,,l], na.rm = TRUE))), na.rm = TRUE) * fishable # pull out total revenue for fleet l
        
      }
      
      total_effort <- sum(fleets[[l]]$e_p_s[, s - 1] * concentrator)
      
      if (fleets[[l]]$fleet_model == "open access"){
        
        if (is.na(fleets[[l]]$cost_per_unit_effort) | is.na(fleets[[l]]$profit_sensitivity)){
          stop("open access fleet model requires both cost_per_unit_effort and profit_sensitivity parameters")
        }

        if (s > 2){ # no past revenue available in first two time steps for accoutning
        

          
          last_revenue <-  sum(last_r_p, na.rm = TRUE) # pull out total revenue for fleet l
          
          
        # last_revenue <-  sum((sapply(storage[[s-2]], function(x) sum(x$r_p_a_fl[,,l], na.rm = TRUE))), na.rm = TRUE) # pull out total revenue for fleet l
        # 
        last_profits <- last_revenue - fleets[[l]]$cost_per_unit_effort * sum(fleets[[l]]$e_p_s[, s - 1])^2 # calculate profits in the last time step
        # 
        total_effort <- total_effort * exp(fleets[[l]]$profit_sensitivity * last_profits) # adjust effort per open access

        }
        
      }
      
      e_p <- fleets[[l]]$e_p_s[, s - 1]
      
      if (fleets[[l]]$spatial_allocation == "revenue") {
        
        if (sum(fishable) == 0){
          alloc <- 0
        } else if (sum(last_r_p, na.rm = TRUE) == 0){
          # if there is no revenue anywhere just distribute fleet evenly as an edge case for extreme overfishing
          
          alloc <- fishable / sum(fishable)
            
            #1 / nrow(r_p_f)
        } else {
          
          # alloc <- (last_r_p / sum(last_r_p, na.rm = TRUE)) / e_p # just extra cautios.
          # 
          # alloc[!is.finite(alloc)] <-  0
          # 
          # alloc <-  alloc / sum(alloc)
          
          alloc <-
            (last_r_p/ sum(last_r_p, na.rm = TRUE)) # just extra cautios.

          
        }
        
        fleets[[l]]$e_p_s[, s] <-
          total_effort * alloc # distribute fishing effort by fishable biomass
        
        
      } else if (fleets[[l]]$spatial_allocation == "rpue") {
        
        
        if (sum(fishable) == 0){
          alloc <- 0
        } else if (sum(last_r_p, na.rm = TRUE) == 0){
          # if there is no revenue anywhere just distribute fleet evenly as an edge case for extreme overfishing
          
          alloc <- fishable / sum(fishable)
          
          #1 / nrow(r_p_f)
        } else {
          
          alloc <- (last_r_p / sum(last_r_p, na.rm = TRUE)) / e_p # just extra cautios.
          
        }
        
        fleets[[l]]$e_p_s[, s] <-
          total_effort * alloc # distribute fishing effort by fishable biomass
        
        
        
      } else if (fleets[[l]]$spatial_allocation == "ideal_free" &&
                 !is.na(fleets[[l]]$cost_per_unit_effort)) {
        
        
        
        stop("ideal free distribution not yet supported. Set spatial_allocation = 'revenue' in create_fleet")

        # calculate expected marginal revenue when effort = 0 in each patch
        # marginal revenue is fishable revenue (r_p_f) - marginal cost per unit effort
        
        worth_fishing <-
          (last_r_p - fleets[[l]]$cost_per_unit_effort) > 0 # check whether any effort could be positive, factoring in potential for negative price,
        
        # for patches that will support any fishing, solve for effort such that marginal profits are equal in space
        
        # optfoo <- function(log_effort, mp = 0, revs, qs, cost) {
        #   mp_hat <- sum(revs * exp(-qs * exp(log_effort))) - cost
        #   
        #   ss <- (mp_hat - mp) ^ 2
        # }
        
        id_e_p <- rep(0, patches)
        
        fishable_patches <- (1:patches)[(fishable == 1) & worth_fishing]
        
        stop("broken here")
        r_p_f <-
          matrix(r_p_f[fishable_patches, ], nrow  = length(fishable_patches)) # seriously annoying step to preserve matrix structure when there is only one species
        
        id_e_p[fishable_patches] <- ((rowSums(log(r_p_f + 1))) - log(fleets[[l]]$cost_per_unit_effort)) / sum(f_q) # see TNC notebook for derivation of this
        
        
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
        
        choices <- dplyr::arrange(data.frame(patch = 1:patches, effort = id_e_p), desc(effort))
        
        choices$cumu_effort <- cumsum(choices$effort)
        
        choices <- choices[choices$cumu_effort < total_effort,]
        
        choices$alloc_effort <- total_effort * (choices$effort / sum(choices$effort))
        
        # choices <- data.frame(patch = 1:patches, effort = id_e_p) %>% 
        #   dplyr::arrange(desc(effort)) %>% 
        #   dplyr::mutate(cumu_effort = cumsum(effort)) %>% 
        #   dplyr::filter(cumu_effort <= total_effort) %>% 
        #   dplyr::mutate(alloc_effort = total_effort * (effort / sum(effort))) # for now, force model to use the same total effort no matter what, so marginal profits will not be zero if the models says not to fish at all, etc. 
        # 
        fleets[[l]]$e_p_s[,s] <- 0
        
        fleets[[l]]$e_p_s[choices$patch, s] <- choices$alloc_effort
        
        
      } # close ifd else
      else {
        stop(
          "spatial effort allocation strategy not properly defined, check spatial_allocation and cost_per_unit_effort in fleet object"
        )
      }
      
    } # close allocate fleets in space based on fishable revenue
    
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
            byrow = TRUE
          ) *
          matrix(
            fleets[[l]]$metiers[[fauni[f]]]$sel_at_age,
            nrow = nrow(last_n_p_a),
            ncol = ncol(last_n_p_a),
            byrow = TRUE
          )
        f_p_a <-
          f_p_a + fleets[[l]]$e_p_s[, s] * tmp
        
        f_p_a_fl[, , l] <-
          fleets[[l]]$e_p_s[, s] * tmp
        
        p_p_a_fl[, , l] <- fleets[[l]]$metiers[[fauni[f]]]$price
        
      } # calculate cumulative f at age by patch
      # you can build a series of if statements here to sub in the correct species module
      f_p_a_fl <-
        f_p_a_fl / array(
          f_p_a,
          dim = c(patches, ages, length(fleets)),
          dimnames = list(1:patches, fauna[[f]]$ages, names(fleets))
        ) # f by patch, age, and fleet
      # season <- (season / self$seasons) - self$time_step
      
      pop <-
        fauna[[f]]$swim(
          season = (season + fauna[[f]]$time_step) * fauna[[f]]$seasons,
          # annoying step to get seasons back to which season is it, not decimal season
          f_p_a = f_p_a,
          last_n_p_a = last_n_p_a
        )
      
      
      
      # process catch data
      c_p_a_fl <-
        f_p_a_fl * array(
          pop$c_p_a,
          dim = c(patches, ages, length(fleets)),
          dimnames = list(1:patches, fauna[[f]]$ages, names(fleets))
        )
      
      r_p_a_fl <- c_p_a_fl * p_p_a_fl # revenue
      
      # storage[[s - 1]][[f]]$rpue_p_a_fl <- r_p_a_fl / 
        
      
      storage[[s - 1]][[f]]$c_p_a_fl <-
        c_p_a_fl # catch stored in each model is the catch that came from the last time step, so put in the right place here
      
      storage[[s - 1]][[f]]$r_p_a_fl <-
        r_p_a_fl # revenue stored in each model is the revenue that came from the last time step, so put in the right place here
      
      storage[[s - 1]][[f]]$c_p_a <-
        pop$c_p_a # catch stored in each model is the catch that came from the last time step, so put in the right place here
      
      # storage[[s - 1]][[f]]$f_p_a_fl <-
      #   f_p_a_fl # proportion of f by fleet
      
      tmp_e_p_fl = purrr::map_dfc(fleets, ~ .x$e_p_s[, s])
      
      storage[[s - 1]][[f]]$e_p_fl <-
        tmp_e_p_fl # store effort by patch by fleet (note that this is the same across species)
      storage[[s]][[f]] <- pop
      
    } # close fauni, much faster this way than dopar, who knew
    
  } #close steps
  
  # Sys.time() - a
  storage <-
    storage[1:(steps - 1)] # since catch is retrospective, chop off last time step to ensure that every step has a catch history
  storage <-
    rlang::set_names(storage, nm = step_names[1:(steps - 1)])
  
  storage <- purrr::map(storage, ~ rlang::set_names(.x, fauni))
  
  
} # close function

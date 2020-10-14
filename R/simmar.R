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
                   tune_unfished = 0) {
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
    (years + 1) / time_step #tack on extra year for accounting
  
  patches <- unique(purrr::map_dbl(fauna, "patches"))
  
  initial_conditions <-
    purrr::map(fauna, c("unfished")) # pull out unfished conditions created by create_critter
  
  if (length(patches) > 1) {
    stop(
      "fauna have different habitat resolutions: set resolution to same number for all species!"
    )
  }
  
  storage <- vector("list", steps)
  
  storage[[1]] <-
    initial_conditions # start populations at initial conditions
  
  fleets <-
    purrr::map(fleets, ~ purrr::list_modify(
      .x,
      e_p_s = matrix(.x$base_effort / patches, nrow = patches, ncol = steps)
    )) # create blank for effort by fleet, space, and time
  
  r_p_f <-
    matrix(0, patches, length(fauni)) # fishable revenue by patch and fauna
  
  e_p_f <-
    matrix(0, patches, length(fauni)) # total fishing mortality by patch and fauna
  
  f_q <- rep(0, length(fauni)) # storage for q by fauna
  
  fishable <- rep(1, patches)
  
  step_seq <-
    seq(1, steps, by = time_step) # create sequence of seasonal time steps
  
  # loop over steps
  for (s in 2:steps) {
    season <-
      step_seq[s - 1] - floor(step_seq[s - 1]) # determine what season the last time step was
    
    year <-  floor(step_seq[s])
    
    if (length(mpas) > 0) {
      # assign MPAs if needed
      if (year == mpas$mpa_year) {
        fishable <- mpas$locations$mpa == 0
        
      }
      
    } # close MPA if statement
    for (l in seq_along(fleet_names)) {
      # distribute fleets in space based on revenues
      
      concentrator <- rep(1, patches) # reset fishing effort concentrator by fleet
      
      for (f in seq_along(fauni)) {
        last_b_p_a <- storage[[s - 1]][[f]]$b_p_a
        
        if (length(mpas) > 0) {

          if (year >= mpas$mpa_year & fleets[[l]]$mpa_response == "leave") {
            concentrator <- as.numeric(fishable)
          }
        }
        
        # calculate fishable biomass in each patch for each species for that fleet
        
        # tmp <-
        #   matrix((1 - exp(
        #     -(fleets[[l]]$metiers[[fauni[f]]]$catchability * fleets[[l]]$metiers[[fauni[f]]]$sel_at_age)
        #   )),
        #   nrow = nrow(last_b_p_a),
        #   ncol = ncol(last_b_p_a),
        #   byrow = TRUE
        #   )
        
        
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
        
        r_p_f[, f] <- last_b_p * fleets[[l]]$metiers[[fauni[f]]]$price
        
        f_q[f] <- fleets[[l]]$metiers[[fauni[f]]]$catchability
        
      } # close fauni loop

      fleets[[l]]$e_p_s[, s] <-
        sum(fleets[[l]]$e_p_s[, s - 1] * concentrator) * pmax(rowSums(r_p_f), 0) / max(sum(pmax(rowSums(r_p_f),0)), 1e-6) # distribute fishing effort by fishable biomass

    } # close allocate fleets in space based on fishable revenue
    
    for (f in seq_along(fauni)) {
      # run population model for each species
      
      ages <-  length(fauna[[f]]$length_at_age)
      
      last_n_p_a <-
        storage[[s - 1]][[f]]$n_p_a # numbers by patch and age in the last time step
      
      f_p_a <-
        matrix(0, nrow = patches, ncol = ages) # total fishing mortality by patch and age
      
      f_p_a_fl <-
        array(0, dim = c(patches, ages, length(fleets)),
              dimnames = list(1:patches, fauna[[f]]$ages, names(fleets))) # storage for proportion of fishing mortality by patch, age, and fleet

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
      } # calculate cumulative f at age by patch
      # you can build a series of if statements here to sub in the correct species module
      f_p_a_fl <-
        f_p_a_fl / array(f_p_a, dim = c(patches, ages, length(fleets)),
                         dimnames = list(1:patches, fauna[[f]]$ages, names(fleets))) # f by patch, age, and fleet
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
        f_p_a_fl * array(pop$c_p_a, dim = c(patches, ages, length(fleets)),
                        dimnames = list(1:patches, fauna[[f]]$ages, names(fleets)))
      storage[[s - 1]][[f]]$c_p_a_fl <-
        c_p_a_fl # catch stored in each model is the catch that came from the last time step, so put in the right place here
      
      storage[[s - 1]][[f]]$c_p_a <-
        pop$c_p_a # catch stored in each model is the catch that came from the last time step, so put in the right place here
      
      # storage[[s - 1]][[f]]$f_p_a_fl <-
      #   f_p_a_fl # proportion of f by fleet

      tmp_e_p_fl = purrr::map_dfc(fleets, ~.x$e_p_s[,s])

      storage[[s - 1]][[f]]$e_p_fl <- tmp_e_p_fl # store effort by patch by fleet (note that this is the same across species)
      storage[[s]][[f]] <- pop
      
    } # close fauni, much faster this way than dopar, who knew
    
  } #close steps
  
  # Sys.time() - a
  
  storage <-
    storage[1:(steps - 1)] # since catch is retrospective, chop off last year to ensure that every step has a catch history
  
  storage <- purrr::map(storage, ~ rlang::set_names(.x, fauni))
  
} # close function

simmar <- function(fauna = list(),
                   fleets = list(),
                   mpas = list(),
                   steps = 100,
                   n_cores = 1) {
  fauni <- names(fauna)
  
  fleet_names <- names(fleets)
  
  patches <- unique(purrr::map_dbl(fauna, "patches"))
  
  initial_conditions <- purrr::map(fauna, c("unfished"))
  
  # a <- matrix(0.01, nrow = patches, ncol = length(fauna[[1]]$length_at_age))
  #
  # a[,1] <- fauna[[1]]$r0 / (patches)
  #
  # initial_conditions$`white seabass`$n_p_a <- a
  
  if (length(patches) > 1) {
    stop(
      "fauna have different habitat resolutions: set resolution to same number for all species!"
    )
  }
  
  # you'll need to fish down to joint EQ here
  
  
  # forecast from joint EQ, basically start MPA here
  
  
  storage <- vector("list", steps)
  
  storage[[1]] <- initial_conditions
  
  
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
  
  # a <- Sys.time()
  
  for (s in 2:steps) {
    if (length(mpas) > 0) {
      if (s == mpas$mpa_step) {
        fishable <- mpas$locations$mpa == 0
      }
      
    } # close MPA if statement
    for (l in seq_along(fleet_names)) {
      for (f in seq_along(fauni)) {
        last_b_p_a <- storage[[s - 1]][[f]]$b_p_a
        
        tmp <-
          matrix((1 - exp(
            -(fleets[[l]][[fauni[f]]]$catchability * fleets[[l]][[fauni[f]]]$sel_at_age)
          )),
          nrow = nrow(last_b_p_a),
          ncol = ncol(last_b_p_a),
          byrow = TRUE
          )
        
        last_b_p <- rowSums(last_b_p_a * tmp) * fishable
        
        r_p_f[, f] <- last_b_p * fleets[[l]][[fauni[f]]]$price
        
        f_q[f] <- fleets[[l]][[fauni[f]]]$catchability
        
      } # close fauni loop
  
      fleets[[l]]$e_p_s[, s] <-
        sum(fleets[[l]]$e_p_s[, s - 1]) * rowSums(r_p_f) / pmax(sum(rowSums(r_p_f)), 1e-6) # distribute fishing effort by fishable biomass
      
      
    } # close allocate fleets in space based on fishable revenue
    
    # eff <- fleets$longline$e_p_s[,2]
    #
    # ssb <- rowSums(last_b_p_a)
    #
    # plot(eff, ssb)
    
    # calculate f at age by patch here?
    
    for (f in seq_along(fauni)) {
      # run population model for each species
      ages <-  length(fauna[[f]]$length_at_age)
      
      last_n_p_a <- storage[[s - 1]][[f]]$n_p_a
      
      f_p_a <- matrix(0, nrow = patches, ncol = ages)
      
      for (l in seq_along(fleet_names)) {
        f_p_a <-
          f_p_a + fleets[[l]]$e_p_s[, s] * matrix(rep(fleets[[l]][[fauni[f]]]$catchability * fleets[[l]][[fauni[f]]]$sel_at_age),
                                                  patches,
                                                  ages,
                                                  byrow = TRUE)
        
        # if (s >= (mpas$mpa_step + 10)){
        #   browser()
        # }

        
      } # calculate cumulative f at age by patch
      # you can build a series of if statements here to sub in the correct species module
      pop <- marlin::sim_fish_pop(
        length_at_age = fauna[[f]]$length_at_age,
        weight_at_age = fauna[[f]]$weight_at_age,
        maturity_at_age = fauna[[f]]$maturity_at_age,
        steepness = fauna[[f]]$steepness,
        m = fauna[[f]]$m,
        patches = patches,
        burn_steps = 0,
        r0 = fauna[[f]]$r0,
        ssb0 = fauna[[f]]$ssb0,
        movement = fauna[[f]]$move_mat,
        f_p_a = f_p_a,
        last_n_p_a = last_n_p_a,
        tune_unfished = 0
      )
      
      storage[[s]][[f]] <- pop
      
    } # close fauni, much faster this way than dopar
    
  } #close steps
  
  # Sys.time() - a
  
  storage <- purrr::map(storage, ~ rlang::set_names(.x, fauni))
  
  # rec <- map(storage, ~.x[[2]]$n_p_a) %>%
  #   map_df(~tibble(rec = .x[,1]),.id = "i") %>%
  #   mutate(i = as.numeric(i)) %>%
  #   filter(i > 1) %>%
  #   group_by(i) %>%
  #   summarise(recs = sum(rec))
  #
  # ssb <- map(storage, ~.x[[2]]$ssb_p_a)%>%
  #   map_df(~tibble(ssb = rowSums(.x)),.id = "i") %>%
  #   mutate(i = as.numeric(i)) %>%
  #   filter(i > 1) %>%
  #   group_by(i) %>%
  #   summarise(ssb = sum(ssb))
  #
  # plot(ssb$ssb, rec$recs)
  #
  # plot(ssb$ssb)
  #
  # last(ssb$ssb) / fauna[[1]]$ssb0
  #
  # ya <- rowSums(storage[[s]][[1]]$ssb_p_a)
  #
  # check <- expand_grid(x = 1:sqrt(patches), y = 1:sqrt(patches)) %>%
  #   mutate(ssb = ya)
  #
  #  ggplot(check, aes(x,y, fill = ssb)) +
  #   geom_tile() +
  #    scale_fill_viridis_c()
  
  
  
  
} # close function

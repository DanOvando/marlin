#' R6 Class Representing a Fish
#'
#' @description
#' A fish object has all the life history data needed to model a finfish population
#'
#' @details
#' Fish like to swim

Fishold <- R6::R6Class(
  "Fishold",
  lock_objects = FALSE,
  public = list(
    #' @description
    #' fill in fish object
    #'
    #' @param common_name common name
    #' @param scientific_name
    #' @param linf
    #' @param vbk
    #' @param t0
    #' @param cv_len
    #' @param length_units
    #' @param min_age
    #' @param max_age
    #' @param weight_a
    #' @param weight_b
    #' @param weight_units
    #' @param length_50_mature
    #' @param length_95_mature
    #' @param delta_mature
    #' @param age_50_mature
    #' @param age_95_mature
    #' @param age_mature
    #' @param length_mature
    #' @param m
    #' @param steepness
    #' @param r0
    #' @param ssb0
    #' @param density_dependence_form
    #' @param adult_diffusion
    #' @param recruit_diffusion
    #' @param query_fishlife
    #' @param sigma_r
    #' @param rec_ac
    #' @param cores
    #' @param mat_mode
    #' @param default_wb
    #' @param tune_weight
    #' @param density_movement_modifier
    #' @param linf_buffer
    #' @param resolution
    #' @param season_blocks
    #' @param recruit_habitat
    #' @param fished_depletion
    #' @param rec_form
    #' @param burn_years
    #' @param seasonal_hab
    #' @param seasons
    #' @param init_explt
    #' @param pups 
    #' @param fec_form 
    #' @param fec_expo exponent for fecundity at weight relationship, 1 = isometric > 1 hyperallometric
    #' @param get_common_name 
    #' @param explt_type
    initialize = function(common_name = NA,
                          scientific_name = NA,
                          linf = NA,
                          vbk = NA,
                          t0 = -0.1,
                          cv_len = 0.1,
                          length_units = 'cm',
                          min_age = 0,
                          max_age = NA,
                          weight_a = NA,
                          weight_b = NA,
                          pups = 10,
                          weight_units = 'kg',
                          fec_form = "weight",
                          fec_expo = 1,
                          length_50_mature = NA,
                          length_95_mature = NA,
                          delta_mature = .1,
                          age_50_mature = NA,
                          age_95_mature = NA,
                          age_mature = NA,
                          length_mature = NA,
                          m = NA,
                          steepness = 0.8,
                          r0 = 10000,
                          ssb0 = NA,
                          density_dependence_form = 1,
                          adult_diffusion = 4,
                          recruit_diffusion = 10,
                          query_fishlife = T,
                          sigma_r = 0,
                          rec_ac = 0,
                          cores = 4,
                          mat_mode = "age",
                          default_wb = 2.8,
                          tune_weight = FALSE,
                          density_movement_modifier = 1,
                          linf_buffer = 1.2,
                          resolution = NA,
                          seasonal_habitat = list(),
                          season_blocks = list(),
                          recruit_habitat = NA,
                          fished_depletion = 1,
                          rec_form = 1,
                          burn_years = 50,
                          seasonal_hab = NA,
                          seasons = 1,
                          explt_type = "f",
                          init_explt = .1,
                          get_common_name = FALSE) {
      
      seasons <- as.integer(seasons)
      
      if (seasons < 1) {
        stop("seasons must be greater than or equal to 1")
      }
      #
      
      if (length(seasonal_habitat) > 1){
        resolution <- nrow(seasonal_habitat[[1]])
      }
      
      if (length(seasonal_habitat) == 0) {
        seasonal_habitat <-
          purrr::map(1:seq_along(seasons), function(x, res)
            matrix(1, nrow = res, ncol = res), res = resolution)
        
      }
      
      if (length(season_blocks) == 0) {
        
        seasons_per_habitat <- seasons / length(seasonal_habitat) # determine how many seasons to assign to each habitat block
        
        tmp <- data.frame(block =  rep(seq_along(seasonal_habitat), each = seasons_per_habitat), season =  1:seasons)
        
        season_blocks <- vector(mode = "list", length = length(seasonal_habitat))
        
        for (i in seq_along(seasonal_habitat)){
          
          season_blocks[[i]] <- tmp$season[tmp$block == i]
          
        }
        
      }
      
      if (length(seasonal_habitat) != length(season_blocks)) {
        stop("length of seasonal_habitat and movement must equal length of season_blocks")
      }
      
      
      if (length(adult_diffusion) > 1 & length(adult_diffusion) != length(seasonal_habitat)){
        stop("As of now adult movement and seasonal habitat lists but be same length")
      }
      
      # if adult movement is constant, make it same shape as seasonal habitat
      if (length(season_blocks) != length(adult_diffusion) & length(adult_diffusion) == 1) {
        
        
        adult_diffusion <- as.list(rep(adult_diffusion, length(seasonal_habitat)))
        
      }
      
      
      
      time_step = 1 / seasons
      
      burn_steps <-
        burn_years * seasons # convert years to time steps
      
      season_foo <-
        function(x, seasons, time_step) {
          if (any(x > seasons)) {
            stop("habitat seasons must be less than or equal to the number of seasons")
            
          }
          
          y = (x / seasons) - time_step
          
        }
      
      all_habitat_seaons <- purrr::map_dfr(season_blocks, ~ data.frame(season = .x))
      
      if (!all(1:seasons %in% all_habitat_seaons$season)){
        stop("all seasons must be represented in season_blocks")
      }
      
      season_blocks <-
        purrr::map(season_blocks,
                   season_foo,
                   time_step = time_step,
                   seasons = seasons)
      
      if (!is.null(dim(seasonal_habitat[[1]]))) {
        resolution <- nrow(seasonal_habitat[[1]])
        
      }
      
      patches <- resolution ^ 2
      
      if (!is.na(scientific_name) & get_common_name == TRUE) {
        common_name <-
          taxize::sci2comm(scientific_name, db = "ncbi")[[1]][1]
        
      }
      
      if (is.na(scientific_name) &
          !is.na(common_name)) {
        scientific_name <-
          taxize::comm2sci(common_name, db = "ncbi")[[1]][1]
        
      }
      
      # check fishbase -------------
      if (is.na(scientific_name) == F &
          query_fishlife == T) {
        sq <- purrr::safely(purrr::quietly(get_traits))
        
        
        genus_species <-
          stringr::str_split(scientific_name, " ", simplify = T) %>%
          as.data.frame() %>%
          rlang::set_names(c("genus", "species"))
        
        fish_life <- genus_species %>%
          dplyr::mutate(life_traits = purrr::pmap(list(
            Genus = genus, Species = species
          ),
          sq))
        
        if (!is.null(fish_life$life_traits[[1]]$error)) {
          stop("No match in FishLife: check spelling or supply your own life history values")
        }
        
        fish_life <- fish_life %>%
          dplyr::mutate(fish_life_worked = purrr::map(life_traits, 'error') %>% purrr::map_lgl(is.null)) %>%
          dplyr::filter(fish_life_worked) %>%
          dplyr::mutate(life_traits = purrr::map(life_traits, c('result', "result"))) %>%
          tidyr::unnest(cols = life_traits) %>%
          dplyr::mutate(taxa = glue::glue('{genus} {species}')) %>%
          rlang::set_names(tolower)
        
        
        if (weight_units == "kg") {
          fish_life$winfinity <- fish_life$winfinity / 1000
        }
        
        if (tune_weight == T) {
          weight_stan <- "
   data {
    real winf;
    real linf;
  }
    parameters {
    real<lower = 0> wa;
    real<lower = 2.7, upper = 3.2> wb;
    real<lower = 0, upper = 1> sigma;
    }
    transformed parameters{
    real w_hat;
    w_hat = wa*linf^wb;

    }
    model {
    winf ~ normal(w_hat, sigma);
    wb ~ normal(3,.1);
    }
    "
          
          weight_fit <-
            rstan::stan(
              model_code = weight_stan,
              data = list(winf = fish_life$winfinity * 2, linf = fish_life$loo),
              verbose = F,
              cores = cores
            )
          
          weight_fit <-
            broom::tidy(weight_fit) %>%
            dplyr::select(term, estimate) %>%
            tidyr::spread(term, estimate)
        } else{
          weight_fit <-
            dplyr::tibble(wa = fish_life$winfinity / (fish_life$loo ^ default_wb),
                          wb = default_wb)
        }
        # process lengths ---------------------------------------------------------
        
        if (is.na(linf)) {
          linf <- fish_life$loo
        }
        
        if (is.na(vbk)) {
          vbk <- fish_life$k
          
        }
        
        if (is.na(weight_a)) {
          weight_a <- weight_fit$wa
          
        }
        if (is.na(weight_b)) {
          weight_b <- weight_fit$wb
          
        }
        
        if (is.na(max_age)) {
          max_age <- ceiling(fish_life$tmax)
          
        }
        
        if (is.na(age_mature)) {
          age_mature <- fish_life$tm
          
        }
        
        if (is.na(length_mature)) {
          length_mature <- fish_life$lm
          
        }
        
        if (is.na(m)) {
          m <- fish_life$m
          
        }
        
      } #close fishlife query
      
      ages <- seq(min_age, max_age, by = time_step)
      
      self$ages <- ages
      
      length_at_age <-
        linf * (1 - exp(-vbk * (ages - t0)))
      
      # process weight
      
      weight_at_age <-
        weight_a * length_at_age ^ weight_b
      
      
      
      lmat_to_linf_ratio <- length_mature / linf
      
      m_at_age <-
        rep(m, length(weight_at_age)) # place holder to allow for different m at age
      
      #
      
      length_at_age_key <-
        generate_length_at_age_key(
          min_age = min_age,
          max_age = max_age,
          cv = cv_len,
          linf = linf,
          k = vbk,
          t0 = t0,
          time_step = time_step,
          linf_buffer = linf_buffer
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(age, length_bin, p_bin) %>%
        tidyr::spread(length_bin, p_bin) %>%
        dplyr::select(-age)
      
      # process maturity
      if ((is.na(age_50_mature) |
           is.na(age_95_mature)) &
          is.na(age_mature) == F) {
        age_50_mature <- age_mature
        
        age_95_mature <-
          age_50_mature + delta_mature
        
        maturity_at_age <-
          ((1 / (1 + exp(-log(
            19
          ) * ((ages - age_50_mature) / (age_95_mature - age_50_mature)
          )))))
        
      } else if (is.na(age_mature) |
                 mat_mode == "length") {
        if (is.na(length_mature)) {
          length_mature <-  linf * lmat_to_linf_ratio
        }
        
        length_bins <-
          as.numeric(colnames(length_at_age_key))
        
        mat_at_bin <- ((1 / (1 + exp(-log(
          19
        ) * ((length_bins - length_mature) / (delta_mature)
        )))))
        
        p_mat_at_age <-
          (as.matrix(length_at_age_key) %*% mat_at_bin)
        
        mat_at_age <-
          dplyr::tibble(age = ages,
                        mean_mat_at_age = p_mat_at_age)
        
        age_mature <-
          mat_at_age$age[mat_at_age$mean_mat_at_age >= 0.5][1]
        
        age_50_mature <- age_mature
        
        age_95_mature <-
          mat_at_age$age[mat_at_age$mean_mat_at_age >= 0.95][1]
        
        maturity_at_age <-
          mat_at_age$mean_mat_at_age
      }
      
      
      if (is.na(length_50_mature)) {
        length_50_mature <- length_mature
        
        length_95_mature <-
          length_50_mature + delta_mature
        
      }
      
      if (fec_form == "pups"){
        
        fec_at_age <- rep(pups, length(maturity_at_age))
        
        self$fec_at_age <- fec_at_age
        
      } else {
        
        
        fec_at_age <- weight_at_age^fec_expo
        
        self$fec_at_age <- fec_at_age
      }
      
      self$maturity_at_age <- maturity_at_age
      
      self$weight_at_age <- weight_at_age
      
      self$ssb_at_age <-
        maturity_at_age * weight_at_age
      
      # create habitat and movement matrices
      
      
      
      
      calc_diffusion_mat <-
        function(diffusion_rate,
                 time_step,
                 resolution) {
          
          # reminder, 
          # 
          # .2^2 == exp(2 * log(.2)), hence strange notation in Thorson et al. going back and forth
          
          
          adjacent <-
            tidyr::expand_grid(x = 1:resolution, y = 1:resolution) %>%
            dist() %>%
            as.matrix()
          
          
          adjacent[adjacent != 1] <- 0
          
          diffusion <- adjacent * diffusion_rate
          
          diag(diffusion) <-  -1 * colSums(diffusion)
          
          diffusion <- as.matrix(diffusion)
          
          # as.matrix(Matrix::expm(diffusion)) %>% image()
          return(diffusion)
          
        } # close calc_move_mat
      self$seasonal_diffusion <-
        purrr::pmap(
          list(diffusion_rate = adult_diffusion),
          calc_diffusion_mat,
          resolution = resolution,
          time_step = time_step
        )
      
      # ideally, you would set things up with mean environmental conditions, but for now, set up a placeholder for movement ignoring taxis for unfished conditions... 
      
      tmp_movement <- lapply(self$seasonal_diffusion, function(x, seasons) as.matrix(Matrix::expm(x * 1 / seasons)), seasons = seasons)
      
      self$movement_seasons <- season_blocks
      
      self$recruit_movement <-
        calc_diffusion_mat(
          diffusion_rate = recruit_diffusion,
          resolution = resolution,
          time_step = time_step
        )
      
      self$recruit_movement <-as.matrix(Matrix::expm(self$recruit_movement * 1 / seasons))
      
      
      # set up unfished recruitment by patch
      if (is.null(dim(recruit_habitat))) {
        recruit_habitat <- matrix(1, nrow = resolution, ncol = resolution)
        
      }
      
      r0s <- recruit_habitat %>%
        as.data.frame() %>%
        dplyr::mutate(x = 1:nrow(.)) %>%
        tidyr::pivot_longer(
          -x,
          names_to = "y",
          values_to = "rec_habitat",
          names_prefix = "V",
          names_ptypes = list(rec_habitat = integer())
        ) %>%
        dplyr::mutate(rec_habitat = rec_habitat / sum(rec_habitat))
      
      local_r0s <- r0 * r0s$rec_habitat
      
      
      
      if (!is.na(ssb0)){
        
        
        tune_ssb0 <-
          function(r0,
                   ssb0_target,
                   rec_habitat,
                   m,
                   max_age,
                   time_step,
                   patches,
                   length_at_age,
                   weight_at_age,
                   maturity_at_age) {
            tmp_r0s <- r0 * rec_habitat
            
            init_pop <-
              tmp_r0s * matrix(
                rep(exp(-m * seq(0, max_age, by = time_step)), patches),
                nrow = patches,
                ncol = length(length_at_age),
                byrow = TRUE
              )
            
            ssb0 <-
              sum(colSums(init_pop) * fec_at_age * maturity_at_age)
            
            delta <- ((ssb0) - (ssb0_target)) ^ 2
            
          }
        tuned_r0 <- nlminb(
          r0,
          tune_ssb0,
          lower = 1e-3,
          ssb0_target = ssb0,
          rec_habitat = r0s$rec_habitat,
          m = m,
          max_age = max_age,
          time_step = time_step,
          patches = patches,
          length_at_age = length_at_age,
          weight_at_age = weight_at_age,
          maturity_at_age = maturity_at_age
        )
        
        local_r0s <- tuned_r0$par * r0s$rec_habitat
        
      }
      
      
      # tune SSB0 and unfished equilibrium
      init_pop <-
        local_r0s * matrix(
          rep(exp(-m * seq(0, max_age, by = time_step)), patches),
          nrow = patches,
          ncol = length(length_at_age),
          byrow = TRUE
        )
      
      self$r0s <- local_r0s
      
      
      f_p_a <-
        matrix(0, nrow = patches, ncol = length(length_at_age))
      
      self$max_age <- max_age # in years
      
      self$min_age <- min_age # starting age for the model
      
      self$linf <- linf
      
      self$m <- m
      
      self$density_dependence_form <-
        density_dependence_form
      
      self$time_step <- time_step
      
      self$seasons <- seasons
      
      self$rec_form <- rec_form
      
      self$common_name <- common_name
      
      self$scientific_name <- scientific_name
      
      self$length_at_age <- length_at_age
      
      self$length_at_age_key <- length_at_age_key
      
      self$length_50_mature <- length_50_mature
      
      self$patches <- patches
      
      self$steepness <- steepness
      
      self$m_at_age <- m_at_age
      
      self$fished_depletion <- fished_depletion
      
      self$explt_type <- explt_type
      
      self$init_explt <- init_explt
      
      unfished <- marlin::sim_fish(
        length_at_age = length_at_age,
        weight_at_age = weight_at_age,
        fec_at_age = fec_at_age,
        maturity_at_age = maturity_at_age,
        steepness = steepness,
        m_at_age = m_at_age,
        patches = resolution ^ 2,
        burn_steps = burn_steps,
        time_step = time_step,
        season = 0,
        r0s = self$r0s,
        ssb0 = NA,
        ssb0_p = rep(-999, patches),
        f_p_a = f_p_a,
        seasonal_movement = tmp_movement,
        movement_seasons = self$movement_seasons,
        recruit_movement = self$recruit_movement,
        last_n_p_a = init_pop,
        tune_unfished = 1,
        rec_form = rec_form
      )
      
      unfished$tmppop$ages <- ages
      
      self$b0 <- sum(unfished$tmppop$b_p_a)
      
      self$ssb0 <- unfished$ssb0
      
      self$ssb0_p <- unfished$ssb0_p
      
      self$n_p_a_0 <- unfished$tmppop$n_p_a
      
      self$unfished <- unfished$tmppop
      
      self$ref_points <- NA
      
    }, # close initialize
  plot = function(type = 2) {
    
    tmp <- as.list(self)
    
    # ogives <- tmp[which(str_detect(names(tmp), "at_age$"))]
    
    ogives <- tmp[which(grepl("at_age$",names(tmp)))]
    
    
    tidy_ogives <-
      purrr::map_df(ogives, ~ data.frame(
        val_at_age = .x,
        age = seq(self$min_age, self$max_age, by = self$time_step)
      ), .id = "trait")
    
    tidy_ogives %>%
      ggplot2::ggplot(aes(age, val_at_age, color = trait)) +
      ggplot2::geom_line(show.legend = FALSE) +
      ggplot2::facet_wrap( ~ trait, scales = "free_y") + 
      marlin::theme_marlin() + 
      ggplot2::scale_fill_manual(values = marlin::marlin_pal("diverging_fish")(length(unique(tidy_ogives$trait))))    
    
  },
  #' Swim
  #'
  #' Swim advances the population one time step
  #'
  #' @param burn_steps number of steps for burn in period if applicable
  #' @param season the current season
  #' @param f_p_a matrix of fishing mortality by patch and age
  #' @param last_n_p_a matrix of initial numbers by patch and age
  #' @param tune_unfished boolean indicating whether to tune unfished
  #'
  #' @return the population in the next time step
  swim = function(burn_steps = 0,
                  season = 1,
                  f_p_a = NULL,
                  last_n_p_a = NULL,
                  adult_movement = NULL,
                  tune_unfished = 0) {
    
    season <- (season / self$seasons) - self$time_step
    
    if (is.null(f_p_a)) {
      f_p_a <-
        matrix(0,
               nrow = self$patches,
               ncol = length(self$length_at_age))
      
    }
    
    if (is.null(last_n_p_a)) {
      last_n_p_a <- self$n_p_a_0
    }
    
    pop <- marlin::sim_fish(
      length_at_age = self$length_at_age,
      weight_at_age = self$weight_at_age,
      fec_at_age = self$fec_at_age,
      maturity_at_age = self$maturity_at_age,
      steepness = self$steepness,
      m_at_age = self$m_at_age,
      patches = self$patches,
      burn_steps = burn_steps,
      time_step = self$time_step,
      season = season,
      r0s = self$r0s,
      ssb0 = self$ssb0,
      ssb0_p = self$ssb0_p,
      seasonal_movement = adult_movement,
      movement_seasons = self$movement_seasons,
      recruit_movement = self$recruit_movement,
      f_p_a = f_p_a,
      last_n_p_a = last_n_p_a,
      tune_unfished = tune_unfished,
      rec_form = self$rec_form
    )
    pop$ages <- self$ages
    
    return(pop)
    
  } # close swim
  
  ) # close public
) # close object
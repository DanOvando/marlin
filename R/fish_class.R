#' R6 Class Representing a Fish
#'
#' @description
#' A fish object has all the life history data needed to model a finfish population
#'
#' @details
#' Fish like to swim

Fish <- R6::R6Class(
  "Fish",
  lock_objects = FALSE,
  public = list(
    #' @description
    #' Used to store parameters describing the behavior of an object of class Fish
    #'
    #' @param common_name common name of the species (can be used to lookup information using `taxize` and `fishlife`)
    #' @param scientific_name common name of the species (can be used to lookup information in `fishlife`)
    #' @param linf asymptotic length of the species in a von Bertalanffy growth function
    #' @param vbk  growth parameter *k* of the species in a von Bertalanffy growth function
    #' @param t0  hypothetical age at which the fish would have length 0 (e.g. -0.5 years)
    #' @param cv_len coefficient of variation around the length at age relationship (in log space)
    #' @param length_units  units of the length at age function (arbitrary)
    #' @param min_age  minimum age tracked in the model. Best to leave at 0, as the model does not explicitly track delays for recruitment
    #' @param max_age  maximum age tracked by the model (individuals this age or older are tracked in the plus group)
    #' @param weight_a  alpha parameter in the allometric weight function alpha x length ^ beta
    #' @param weight_b  beta parameter in the allometric weight function alpha x length ^ beta
    #' @param weight_units  units of the allometric weight function (defaults to kg)
    #' @param length_50_mature  length at 50% maturity in a logistic maturity ogive
    #' @param length_95_mature  length at 95% maturity in a logistic maturity ogive
    #' @param delta_mature as an alternative, the different in units of length between length_50_mature and length_95_mature
    #' @param age_50_mature  age at 50% maturity in a logistic maturity ogive if maturity is age based
    #' @param age_95_mature  age at 95% maturity in a logistic maturity ogive if maturity is age based
    #' @param age_mature an alternative option to just set one age mature, which ends up as the age_50_mature
    #' @param m instantaneous natural mortality rate. When `lorenzen_m = TRUE`, this is the average natural mortality across all ages
    #' @param steepness  steepness parameter (h) in a Beverton-Holt spawner-recruit function
    #' @param r0  asymptotic number of recruits under unfished conditions
    #' @param ssb0  asymptotic spawning stock biomass of recruits under unfished conditions. Tunnes r0 to achieve
    #' @param density_dependence  timing and nature of density dependence in the Beverton-Holt spawner recruit function, one of one of 'global_habitat','local_habitat','pre_dispersal','post_dispersal','global_ssb'
    #' @param adult_diffusion  diffusion parameter *D* in the CTMC movement function for "adults" (not recruits)
    #' @param recruit_diffusion  diffusion parameter *D* in the CTMC movement function for recruits
    #' @param query_fishlife TRUE or FALSE to query `Fishlife` for missing life history values. When set to FALSE all required life history values must be supplied by the user
    #' @param sigma_r the standard deviation of recruitment deviates in log-normal space
    #' @param rec_ac the autocorrelation of recruitment deviates
    #' @param cores the number of cores used to tun the weight relationship if used (deprecated)
    #' @param mat_mode specifies whether maturity is a function of age (default) or length
    #' @param default_wb deprecated
    #' @param tune_weight deprecated
    #' @param linf_buffer multiplier around linf to create length at age key, taking into account that some fish will be larger than Linf
    #' @param resolution a vector of length two with number of patches in X and Y dimensions
    #' @param season_blocks list with elements indicating blocks of seasons. For example, if there are four seasons, setting `season_blocks = list(c(1,2),c(3,4))` indicates that seasons 1 and 2 are one block, 3 and 4 another. Allows for the model to be run at fine time scales while allowing some processes like movement or spawning to operate at coarser scales
    #' @param recruit_habitat a matrix with dimensions X and Y with quality of recruit habitat (scales r0 in space as well as recruit diffusion under applicable forms of density dependence)
    #' @param fished_depletion  depletion (SSB/SSB0) under initial fished conditions
    #' @param burn_years  number of years used to burn in the population to tune parameters without analytical solutions like SSB0
    #' @param seasons  number of seasons per year. 4 would indicate quarterly time steps, 12 monthly, 365 daily.
    #' @param init_explt  instantaneous fishing mortality rate under initial fished conditions
    #' @param pups  number of pups per individual for animals with pups rather than larvae
    #' @param fec_form one of of  "weight" (default) or "pups". When "weight", fecundity is a function of weight. When "pups", constant number of pups per individual produced
    #' @param fec_expo exponent for fecundity at weight relationship, 1 = isometric > 1 hyperallometric
    #' @param get_common_name TRUE or FALSE to lookup common name from scientific name. Requires internet connection
    #' @param habitat a matrix with dimensions X and Y specifying quality of adult (non-recruit) habitat. Determines taxis matrix
    #' @param spawning_seasons  seasons in which spawning occurs
    #' @param patch_area  area of each patch (e.g. KM2)
    #' @param max_hab_mult maximum value of that habitat matrix multiplier (to prevent some habitats from being >>> good, defaults to 2)
    #' @param lorenzen_m TRUE or FALSE to use the Lorenzen function to calculate natural mortality at age
    #' @param explt_type deprecated
    initialize = function(common_name = NA,
                          scientific_name = NA,
                          linf = NA,
                          vbk = NA,
                          t0 = -0.5,
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
                          m = NA,
                          steepness = 0.8,
                          r0 = 10000,
                          ssb0 = NA,
                          density_dependence = "global_habitat",
                          adult_diffusion = 4,
                          recruit_diffusion = 10,
                          query_fishlife = T,
                          sigma_r = 0,
                          rec_ac = 0,
                          cores = 4,
                          mat_mode = "age",
                          default_wb = 2.8,
                          tune_weight = FALSE,
                          linf_buffer = 1.2,
                          resolution = NA,
                          patch_area = 1,
                          habitat = list(),
                          season_blocks = list(),
                          recruit_habitat = NA,
                          fished_depletion = 1,
                          burn_years = 50,
                          seasons = 1,
                          explt_type = "f",
                          init_explt = .1,
                          get_common_name = FALSE,
                          spawning_seasons = NA,
                          max_hab_mult = 2,
                          lorenzen_m = TRUE) {
      seasons <- as.integer(seasons)

      # if only one resolution dimension provided, assume square


      if (seasons < 1) {
        stop("seasons must be greater than or equal to 1")
      }

      if (all(is.na(spawning_seasons))) {
        spawning_seasons <- 1:seasons
      }

      # spawning_seasons <- (spawning_seasons / seasons) - 1 / seasons
      if (length(habitat)>0 && class(habitat) == "list") {
        resolution <- c(nrow(habitat[[1]]), ncol(habitat[[1]]))
      } else if (any(class(habitat) == "matrix")){
        resolution <- c(nrow(habitat), ncol(habitat))
      }

      # if habitat is an empty list
      if (length(resolution) == 1 & all(!is.na(resolution))){
        resolution <- rep(resolution,2)
      }

      self$resolution <- resolution

      patches <- prod(resolution)

      self$patches <- patches

      self$closest_taxa_match <- NA

      if (purrr::is_empty(habitat)) {
        habitat <-
          purrr::map(1:seq_along(seasons), function(x, res)
            matrix(0, nrow = res[1], ncol = res[2]), res = self$resolution)

      }

      if (length(season_blocks) == 0) {
        seasons_per_habitat <-
          seasons / length(habitat) # determine how many seasons to assign to each habitat block

        tmp <-
          data.frame(
            block =  rep(seq_along(habitat), each = seasons_per_habitat),
            season =  1:seasons
          )

        season_blocks <-
          vector(mode = "list", length = length(habitat))

        for (i in seq_along(habitat)) {
          season_blocks[[i]] <- tmp$season[tmp$block == i]

        }

      }

      if (length(habitat) != length(season_blocks)) {
        stop("length of habitat and movement must equal length of season_blocks")
      }


      if (length(adult_diffusion) > 1 &
          length(adult_diffusion) != length(habitat)) {
        stop("As of now adult movement and seasonal habitat lists but be same length")
      }

      # if adult movement is constant, make it same shape as seasonal habitat
      if (length(season_blocks) != length(adult_diffusion) &
          length(adult_diffusion) == 1) {
        adult_diffusion <-
          as.list(rep(adult_diffusion, length(season_blocks)))

      }

      time_step = 1 / seasons

      burn_steps <-
        burn_years * seasons # convert years to time steps

      season_foo <-
        function(x, seasons, time_step) {
          if (any(x > seasons)) {
            stop("habitat seasons must be less than or equal to the number of seasons")

          }

          y = x #/(x / seasons) - time_step

        }

      all_habitat_seasons <-
        purrr::map(season_blocks, ~ data.frame(season = .x)) |>
        purrr::list_rbind()

      if (!all(1:seasons %in% all_habitat_seasons$season)) {
        stop("all seasons must be represented in season_blocks")
      }

      season_blocks <-
        purrr::map(season_blocks,
                   season_foo,
                   time_step = time_step,
                   seasons = seasons)

      # if (!is.null(dim(habitat[[1]]))) {
      #   resolution <- c(nrow(habitat[[1]]), nrow(habitat[[2]]))
      #
      # }
      #
      # es <- prod(resolution)

      if (!is.na(scientific_name) & get_common_name == TRUE) {
        common_name <-
          taxize::sci2comm(scientific_name, db = "worms")[[1]][1]

      }

      if (is.na(scientific_name) &
          !is.na(common_name)) {
        scientific_name <-
          taxize::comm2sci(common_name, db = "worms")[[1]][1]

      }

      # check fishlife -------------
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
          stop("No match in taxize::classification or there was a problem connecting to the taxize database: check network connection/spelling or supply your own life history values")
        }
        fish_life <- fish_life %>%
          dplyr::mutate(fish_life_worked = purrr::map(life_traits, 'error') %>% purrr::map_lgl(is.null)) %>%
          dplyr::filter(fish_life_worked) %>%
          dplyr::mutate(life_traits = purrr::map(life_traits, c('result', "result"))) %>%
          tidyr::unnest(cols = life_traits) %>%
          dplyr::mutate(taxa = glue::glue('{genus} {species}')) %>%
          rlang::set_names(tolower)

        self$closest_taxa_match <- fish_life$closest_taxa_match

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


        if (is.na(length_50_mature)) {
          length_50_mature <- fish_life$lm

        }

        if (is.na(m)) {
          m <- fish_life$m

        }

      } #close fishlife query

      # if there still is no max age, guess it based on natural mortality
     if (is.na(max_age)){
        max_age <-  -log(0.05) / m
     }
      ages <- seq(min_age, max_age, by = time_step)

      self$ages <- ages

      length_at_age <-
        linf * (1 - exp(-vbk * (ages - t0)))

      # process weight

      weight_at_age <-
        weight_a * length_at_age ^ weight_b



      lmat_to_linf_ratio <- length_50_mature / linf

      if (lorenzen_m){

        m_at_age <- (length_at_age / max(length_at_age)) ^ -1

        m_at_age <- (m_at_age / mean(m_at_age)) * m

      } else {
        m_at_age <-
          rep(m, length(weight_at_age)) # place holder to allow for different m at age

      }

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


        # find first length at which 50%/95% of animals are fully mature

        length_bins <-
          as.numeric(colnames(length_at_age_key))

          length_50_mature <- length_bins[which.min((cumsum(as.matrix(length_at_age_key)[which.min((ages - (age_95_mature))^2),]) - 0.5)^2)[1]]

          length_95_mature <- length_bins[which.min((cumsum(as.matrix(length_at_age_key)[which.min((ages - (age_95_mature))^2),]) - 0.95)^2)[1]]

      } else if (is.na(age_mature) |
                 mat_mode == "length") {
        if (is.na(length_50_mature)) {
          length_50_mature <-  linf * lmat_to_linf_ratio
        }

        if (!is.na(length_95_mature)){
          delta_mature <- length_95_mature - length_50_mature
        }

        length_95_mature <-
          length_50_mature + delta_mature

        length_bins <-
          as.numeric(colnames(length_at_age_key))

        mat_at_bin <- ((1 / (1 + exp(-log(
          19
        ) * ((length_bins - length_50_mature) / (delta_mature)
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


      # if (is.na(length_50_mature)) {
      #
      #   length_95_mature <-
      #     length_50_mature + delta_mature
      #
      # }

      if (fec_form == "pups") {
        fec_at_age <- rep(pups, length(maturity_at_age))

        self$fec_at_age <- fec_at_age

      } else {
        fec_at_age <- weight_at_age ^ fec_expo

        self$fec_at_age <- fec_at_age
      }

      self$maturity_at_age <- maturity_at_age

      self$weight_at_age <- weight_at_age

      self$ssb_at_age <-
        maturity_at_age * weight_at_age


      self$max_hab_mult <- max_hab_mult
      # create habitat and movement matrices

      taxis_matrix <- habitat
      # reshape to vector, for some reason doesn't work inside function
      for (i in seq_along(taxis_matrix)) {
        taxis_matrix[[i]] <-
          tidyr::pivot_longer(as.data.frame(taxis_matrix[[i]]), tidyr::everything()) # need to use pivot_longer to match patch order from expand_grid

        taxis_matrix[[i]] <- as.numeric(taxis_matrix[[i]]$value)

        taxis_matrix[[i]] <- pmin(exp((time_step * outer(taxis_matrix[[i]], taxis_matrix[[i]], "-")) / sqrt(patch_area)),max_hab_mult) # convert habitat gradient into diffusion multiplier

      }

      self$taxis_matrix <- taxis_matrix


      diffusion_prep <- function(x,y, time_step, patch_area){

        x[!is.na(x)] <- 1 # this is here to allow for barriers; set diffusion to zero if there's a physical barrier

        z <- x * y * (time_step / patch_area)

      }

      self$diffusion_foundation <- purrr::map2(taxis_matrix, adult_diffusion,diffusion_prep, time_step = time_step, patch_area = patch_area) # prepare adult diffusion matrix account for potential land

      self$adult_diffusion <- adult_diffusion

      diffusion_and_taxis <- purrr::map2(self$diffusion_foundation, taxis_matrix, ~ .x * .y)

      inst_movement_matrix <-  purrr::pmap(list(multiplier = diffusion_and_taxis),
                               prep_movement,
                               resolution = self$resolution)

      # ideally, you would set things up with mean environmental conditions, but for now, set up a placeholder for movement ignoring taxis for unfished conditions...

      self$movement_matrix <-
        purrr::map(inst_movement_matrix,
                    ~ as.matrix(expm::expm((.x))))

      self$movement_seasons <- season_blocks

      rec_diff_foundation <- purrr::map2(taxis_matrix, recruit_diffusion,diffusion_prep, time_step = time_step,patch_area = patch_area) # prepare diffusion matrix account for potential land

      inst_recruit_move_matrix <-
        prep_movement(multiplier = rec_diff_foundation[[1]],
                      resolution = resolution)

      self$recruit_movement_matrix <-
        as.matrix(expm::expm(inst_recruit_move_matrix))


      # set up unfished recruitment by patch
      if (is.null(dim(recruit_habitat))) {
        recruit_habitat <- matrix(1, nrow = resolution[1], ncol = resolution[2])

      }

      recruit_habitat[is.na(recruit_habitat)] <- 0

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



      if (!is.na(ssb0)) {
        tune_ssb0 <-
          function(r0,
                   ssb0_target,
                   rec_habitat,
                   m_at_age,
                   max_age,
                   time_step,
                   patches,
                   length_at_age,
                   weight_at_age,
                   maturity_at_age) {

            tmp_r0s <- r0 * rec_habitat


            n_at_a <- matrix(NA, nrow = patches, ncol = length(m_at_age))

            n_at_a[,1] <- tmp_r0s

            for (i in 2:length(m_at_age)){

              n_at_a[,i] <- n_at_a[,i-1] * exp(-m_at_age[i-1] * time_step)

            } # fill in numbers at age

            n_at_a[,length(m_at_age)] <- n_at_a[,length(m_at_age)] / (1 - exp(-m_at_age[length(m_at_age)] * time_step))

            ssb0 <-
              sum(colSums(n_at_a) * fec_at_age * maturity_at_age)

            delta <- ((ssb0) - (ssb0_target)) ^ 2

          }
        tuned_r0 <- nlminb(
          r0,
          tune_ssb0,
          lower = 1e-3,
          ssb0_target = ssb0,
          rec_habitat = r0s$rec_habitat,
          m_at_age = m_at_age,
          max_age = max_age,
          time_step = time_step,
          patches = patches,
          length_at_age = length_at_age,
          weight_at_age = weight_at_age,
          maturity_at_age = maturity_at_age
        )

        local_r0s <- tuned_r0$par * r0s$rec_habitat

      }



      init_pop <- matrix(NA, nrow = patches, ncol = length(m_at_age))

      init_pop[,1] <- local_r0s
      for (i in 2:length(m_at_age)){

        init_pop[,i] <- init_pop[,i-1] * exp(-m_at_age[i-1] * time_step)

      } # fill in numbers at age

      init_pop[,length(m_at_age)] <- init_pop[,length(m_at_age)] / (1 - exp(-m_at_age[length(m_at_age)] * time_step))

      self$r0s <- local_r0s


      f_p_a <-
        matrix(0, nrow = patches, ncol = length(length_at_age))

      self$max_age <- max_age # in years

      self$min_age <- min_age # starting age for the model

      self$linf <- linf

      self$m <- m

      self$time_step <- time_step

      self$patch_area <- patch_area

      self$seasons <- seasons


      rec_form <- switch(EXPR = density_dependence,"global_habitat" = 0 , "local_habitat" = 1, "pre_dispersal" = 2,"post_dispersal" = 3)

      rec_form <- density_dependence

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

      self$spawning_seasons <- spawning_seasons

      self$sigma_r <- sigma_r

      self$rec_ac <- rec_ac

      unfished <- marlin::sim_fish(
        length_at_age = length_at_age,
        weight_at_age = weight_at_age,
        fec_at_age = fec_at_age,
        maturity_at_age = maturity_at_age,
        steepness = steepness,
        m_at_age = m_at_age,
        patches = prod(resolution),
        burn_steps = burn_steps,
        time_step = time_step,
        season = spawning_seasons[1],
        r0s = self$r0s,
        ssb0 = NA,
        ssb0_p = rep(-999, patches),
        f_p_a = f_p_a,
        movement_matrix =  self$movement_matrix,
        movement_seasons = self$movement_seasons,
        recruit_movement_matrix = self$recruit_movement_matrix,
        last_n_p_a = init_pop,
        tune_unfished = 1,
        rec_form = rec_form,
        spawning_seasons = self$spawning_seasons,
        rec_devs = rep(1, patches)
      )

      unfished$tmppop$ages <- ages

      unfished$tmppop$resolution <- self$resolution

      self$b0 <- sum(unfished$tmppop$b_p_a)

      self$ssb0 <- unfished$ssb0

      self$ssb0_p <- unfished$ssb0_p

      self$n_p_a_0 <- unfished$tmppop$n_p_a

      self$unfished <- unfished$tmppop

      self$ref_points <- NA

    },
  # close initialize
#' plot object of class fish
#'
#' @param type
#'
#' @return a plot of the life history values
  plot = function(type = 2) {
    tmp <- as.list(self)

    # ogives <- tmp[which(str_detect(names(tmp), "at_age$"))]

    ogives <- tmp[which(grepl("at_age$", names(tmp)))]


    tidy_ogives <-
      purrr::map_df(ogives, ~ data.frame(
        val_at_age = .x,
        age = seq(self$min_age, self$max_age, by = self$time_step)
      ), .id = "trait")

    tidy_ogives %>%
      ggplot2::ggplot(aes(age, val_at_age, color = trait)) +
      ggplot2::geom_line(show.legend = FALSE) +
      ggplot2::facet_wrap(~ trait, scales = "free_y") +
      marlin::theme_marlin() +
      ggplot2::scale_fill_manual(values = marlin::marlin_pal("diverging_fish")(length(unique(
        tidy_ogives$trait
      ))))

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
  #' @param adult_movement the adult movement matrix
  #' @param rec_devs externally supplied recruitment deviates
  #'
  #' @return the population in the next time step
  swim = function(burn_steps = 0,
                  season = 1,
                  f_p_a = NULL,
                  last_n_p_a = NULL,
                  adult_movement = NULL,
                  tune_unfished = 0,
                  rec_devs = NA) {

    if (is.null(f_p_a)) {
      f_p_a <-
        matrix(0,
               nrow = self$patches,
               ncol = length(self$length_at_age))

    }

    if (is.null(last_n_p_a)) {
      last_n_p_a <- self$n_p_a_0
    }

    if (all(is.na(rec_devs))){
      rec_devs <- exp(rnorm(self$patches, 0, self$sigma_r) - self$sigma_r^2/2)
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
      movement_matrix = adult_movement,
      movement_seasons = self$movement_seasons,
      recruit_movement_matrix = self$recruit_movement_matrix,
      f_p_a = f_p_a,
      last_n_p_a = last_n_p_a,
      tune_unfished = tune_unfished,
      rec_form = self$rec_form,
      spawning_seasons = self$spawning_seasons,
      rec_devs = rec_devs
    )
    pop$ages <- self$ages

    pop$resolution <- self$resolution
    return(pop)

  } # close swim

  ) # close public
) # close object

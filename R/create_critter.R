#' create_fish
#'
#' creates a fish list object with all the life history goodies
#'
#' @param common_name
#' @param scientific_name
#' @param linf
#' @param vbk
#' @param t0
#' @param max_age
#' @param weight_a
#' @param weight_b
#' @param length_50_mature
#' @param length_95_mature
#' @param age_50_mature
#' @param age_95_mature
#' @param age_mature
#' @param length_mature
#' @param m
#' @param steepness
#' @param density_dependence_form
#' @param adult_movement
#' @param larval_movement
#' @param query_fishlife
#' @param r0
#' @param cv_len
#' @param length_units
#' @param min_age
#' @param time_step
#' @param weight_units
#' @param delta_mature
#' @param price
#' @param sigma_r
#' @param rec_ac
#' @param cores
#' @param mat_mode
#' @param price_ac
#' @param price_cv
#'
#' @return a fish list object
#' @export
#'
#' @examples
#' \dontrun{
#' white_seabass = create_critter(scientific_name = "Atractoscion nobilis", query_fishlife = T)
#'}
create_critter <- function(common_name = 'white seabass',
                           scientific_name = NA,
                           linf = NA,
                           vbk = NA,
                           t0 = -0.1,
                           cv_len = 0.1,
                           length_units = 'cm',
                           min_age = 0,
                           max_age = NA,
                           time_step = 1,
                           weight_a = NA,
                           weight_b = NA,
                           weight_units = 'kg',
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
                           adult_movement = 2,
                           larval_movement = 2,
                           query_fishlife = T,
                           sigma_r = 0,
                           rec_ac = 0,
                           cores = 4,
                           mat_mode = "age",
                           default_wb = 2.8,
                           tune_weight = FALSE,
                           density_movement_modifier = 1,
                           linf_buffer = 1.2,
                           resolution = 25,
                           habitat = NA,
                           fished_depletion = 1) {
  fish <- list()


  if (!is.na(scientific_name)) {
    common_name <-
      taxize::sci2comm(scientific_name, db = "ncbi")[[1]][1]

  }

  if (is.na(scientific_name) & !is.na(common_name)) {
    scientific_name <-
      taxize::comm2sci(common_name, db = "worms")[[1]][1]

  }

  # check fishbase -------------
  if (is.na(scientific_name) == F & query_fishlife == T) {
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

      weight_fit <- broom::tidy(weight_fit) %>%
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

  # max_age <- ((-log(0.01)/m)) %>% floor()

  # if (is.na(vbk)){
  #
  #   vbk <- m / (lhi_groups$mean_m_v_k[lhi_groups$type == lhi_type])
  #
  # }
  # if (is.na(weight_a)){
  #
  #   weight_a <-lhi_groups$mean_wa[lhi_groups$type == lhi_type]
  #
  #   weight_b <-lhi_groups$mean_wb[lhi_groups$type == lhi_type]
  #
  # }

  length_at_age <-
    linf * (1 - exp(-vbk * (seq(
      min_age, max_age, by = time_step
    ) - t0)))

  # process weight

  weight_at_age <- weight_a * length_at_age ^ weight_b

  lmat_to_linf_ratio <- length_mature / linf

  #

  length_at_age_key <- generate_length_at_age_key(
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
       is.na(age_95_mature)) & is.na(age_mature) == F) {
    age_50_mature <- age_mature

    age_95_mature <-  age_50_mature + delta_mature

    maturity_at_age <-
      ((1 / (1 + exp(-log(
        19
      ) * ((seq(min_age, max_age, by = time_step) - age_50_mature) / (age_95_mature - age_50_mature)
      )))))

  } else if (is.na(age_mature) | mat_mode == "length") {
    if (is.na(length_mature)) {
      length_mature <-  linf * lmat_to_linf_ratio
    }

    length_bins <- as.numeric(colnames(length_at_age_key))

    mat_at_bin <- ((1 / (1 + exp(-log(
      19
    ) * ((length_bins - length_mature) / (delta_mature)
    )))))

    p_mat_at_age <- (as.matrix(length_at_age_key) %*% mat_at_bin)

    mat_at_age <-
      dplyr::tibble(age = seq(min_age, max_age, by = time_step),
                    mean_mat_at_age = p_mat_at_age)

    age_mature <-
      mat_at_age$age[mat_at_age$mean_mat_at_age >= 0.5][1]

    age_50_mature <- age_mature

    age_95_mature <-
      mat_at_age$age[mat_at_age$mean_mat_at_age >= 0.95][1]

    maturity_at_age <- mat_at_age$mean_mat_at_age
  }


  if (is.na(length_50_mature)) {
    length_50_mature <- length_mature

    length_95_mature <- length_50_mature + delta_mature

  }


  ssb_at_age <-  maturity_at_age * weight_at_age

  # create habitat and movement matrices

  if (!is.null(dim(habitat))) {
    resolution <- sqrt(nrow(habitat))

  }

  patches <- resolution ^ 2

  distance <-
    tidyr::expand_grid(x = 1:resolution, y = 1:resolution) %>%
    dist() %>%
    as.matrix()

  distance <- dnorm(distance, 0, adult_movement)

  distance <- distance / rowSums(distance)


  if ((!is.null(dim(habitat)))) {
    habitat <- habitat / rowSums(habitat_mat)

    move_mat <-  distance * habitat

    move_mat <- move_mat / rowSums(move_mat)

  } else {
    move_mat <- distance
  }

  move_mat <-
    t(move_mat) # needs to be transposed for use in population function

  # tune SSB0 and unfished equilibrium

  init_pop <-
    matrix(1, nrow = patches, ncol = length(length_at_age))

  init_pop[, 1] <- r0 / patches
  
  f_p_a <- matrix(0, nrow = patches, ncol = length(length_at_age))

  unfished <- marlin::sim_fish_pop(
    length_at_age = length_at_age,
    weight_at_age = weight_at_age,
    maturity_at_age = maturity_at_age,
    steepness = steepness,
    m = m,
    patches = resolution ^ 2,
    burn_steps = 100,
    r0 = r0,
    ssb0 = NA,
    f_p_a = f_p_a,
    movement = move_mat,
    last_n_p_a = init_pop,
    tune_unfished = 1
  )

  ssb0 <- unfished$ssb0

  n_p_a_0 <- unfished$tmppop$n_p_a

  unfished <- unfished$tmppop

  rm(list = c("sq", "f_p_a","weight_fit","distance"))
  fish <- list(mget(ls()))

  fish <- fish[[1]]


  return(fish)
}

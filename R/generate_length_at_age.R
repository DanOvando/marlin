#' generate_length_to_age_key
#'
#' produces an age by length bins matrix with probability of being in length bin at age
#'
#' @param cv the coefficient of variation of the length-at-age relationship (log-space)
#' @param k the
#' @param linf asymptotic length of the species in a von Bertalanffy growth function
#' @param t0 hypothetical age at which the fish would have length 0 (e.g. -0.5 years)
#' @param time_step the time step the model is running on (1 / seasons)
#' @param linf_buffer multiplier around linf to create length at age key, taking into account that some fish will be larger than Linf
#' @param min_age minimum age tracked in the model. Best to leave at 0, as the model does not explicitly track delays for recruitment
#' @param max_age maximum age tracked by the model (individuals this age or older are tracked in the plus group)
#'
#' @return a length-at-age key
#' @export
#'
generate_length_at_age_key <- function(min_age,
                                       max_age,
                                       length_bin_width = 1,
                                       growth_params = list(),
                                       growth_model = "von_bertalanffy",
                                       cv,
                                       time_step = 1,
                                       linf_buffer = 10) {

  if (growth_model == "von_bertalanffy") {
    mean_length_at_age <-
      growth_params$linf * (1 - exp(-growth_params$vbk * (
        seq(min_age, max_age, by = time_step) - growth_params$t0
      )))
  } else if (growth_model == "power") {
    mean_length_at_age <-
      growth_params$length_a *  (seq(min_age, max_age, by = time_step) - growth_params$t0)^growth_params$length_b
    
  }
  
  if (is.null(growth_params$linf)) {
    linf <- max(mean_length_at_age)
  } else {
    linf <- growth_params$linf
  }
  
  length_at_age_vars <- dplyr::tibble(
    age = seq(min_age, max_age, by = time_step),
    mean_length_at_age = mean_length_at_age,
    sigma_at_age = cv * mean_length_at_age
  ) #calculate standard deviation of length at age for each age bin

  # now calculate the probability of being in each length bin at each age
  p_length_at_age <-
    expand.grid(
      age = seq(min_age, max_age, by = time_step),
      length_bin = seq(0,(linf_buffer * linf), by = length_bin_width)
    ) %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(length_at_age_vars, by = 'age') %>%
    dplyr::arrange(age, length_bin)

  p_length_at_age <- p_length_at_age %>%
    dplyr::group_by(age) %>%
    dplyr::mutate(next_length_bin = dplyr::lead(length_bin, 1)) %>%
    dplyr::mutate(p_bin = ifelse(
      is.na(next_length_bin) == F,
      pnorm(next_length_bin, mean_length_at_age, sigma_at_age),
      1
    ) -
      pnorm(length_bin, mean_length_at_age, sigma_at_age))

}

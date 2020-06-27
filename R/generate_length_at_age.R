#' generate_length_to_age_key
#'
#' produces an age by length bins matrix with probability of being in length bin at age
#'
#' @param cv
#' @param k
#' @param linf
#' @param t0
#' @param time_step
#' @param linf_buffer
#' @param min_age
#' @param max_age
#'
#' @return a length-at-age key
#' @export
#'
generate_length_at_age_key <- function(min_age,
                                       max_age,
                                       cv,
                                       k,
                                       linf,
                                       t0,
                                       time_step = 1,
                                       linf_buffer = 10) {



  mean_length_at_age <-
    linf * (1 - exp(-k * (seq(
      min_age, max_age, by = time_step
    ) - t0)))

  length_at_age_vars <- dplyr::tibble(
    age = seq(min_age, max_age, by = time_step),
    mean_length_at_age = mean_length_at_age,
    sigma_at_age = cv * mean_length_at_age
  ) #calculate standard deviation of length at age for each age bin

  # now calculate the probability of being in each length bin at each age
  p_length_at_age <-
    expand.grid(
      age = seq(min_age, max_age, by = time_step),
      length_bin = 0:(linf_buffer * linf)
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

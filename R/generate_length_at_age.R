#' Generate Length-at-Age Key
#'
#' Produces an age-by-length-bin matrix giving the probability of being in each
#' length bin at each age, based on the specified growth model.
#'
#' @param min_age Numeric. Minimum age tracked in the model. Best left at 0,
#'   as the model does not explicitly track recruitment delays.
#' @param max_age Numeric. Maximum age tracked by the model (individuals this
#'   age or older are in the plus group).
#' @param length_bin_width Numeric. Width of each length bin in the key
#'   (default 1).
#' @param growth_params Named list of growth parameters. Contents depend on
#'   \code{growth_model}: for \code{"von_bertalanffy"}: \code{linf}, \code{vbk},
#'   \code{t0}; for \code{"power"}: \code{length_a}, \code{length_b}, \code{t0};
#'   for \code{"growth_cessation"}: \code{l0}, \code{rmax}, \code{k},
#'   \code{t50}.
#' @param growth_model Character. Growth model to use. One of
#'   \code{"von_bertalanffy"}, \code{"power"}, or \code{"growth_cessation"}.
#' @param cv Numeric. Coefficient of variation of length-at-age (log-space).
#' @param time_step Numeric. Time step the model is running on (1 / seasons).
#' @param linf_buffer Numeric. Multiplier around Linf used to set the upper
#'   bound of the length key, accounting for fish larger than Linf (default 10).
#'
#' @return A tibble with columns \code{age}, \code{length_bin},
#'   \code{mean_length_at_age}, \code{sigma_at_age}, \code{next_length_bin},
#'   and \code{p_bin} (probability of being in each length bin at each age).
#' @export
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
      growth_params$length_a * (seq(min_age, max_age, by = time_step) - growth_params$t0)^growth_params$length_b
  } else if (growth_model == "growth_cessation"){

    ages <- seq(min_age, max_age, by = time_step)
    mean_length_at_age <-
     (growth_params$l0 + growth_params$rmax * ((log(exp(-growth_params$k * growth_params$t50) + 1) - log(exp(growth_params$k * (ages - growth_params$t50)) + 1)) / growth_params$k + ages))


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
  ) # calculate standard deviation of length at age for each age bin

  # now calculate the probability of being in each length bin at each age
  p_length_at_age <-
    expand.grid(
      age = seq(min_age, max_age, by = time_step),
      length_bin = seq(0, (linf_buffer * linf), by = length_bin_width)
    ) %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(length_at_age_vars, by = "age") %>%
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

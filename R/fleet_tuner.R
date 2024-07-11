#' fleet tuner
#'
#' finds catchability (q) by fleet such that
#' target fished depletion is achieved
#'
#' @param qs vector of catchability coefficients
#' @param fauna fauna object
#' @param years number of years to tune
#' @param fleets fleet object
#'
#' @return objective function of fleet tuner
#' @export
#'
fleet_tuner <- function(qs,fauna, fleets, years = 50){


  cc <- 1

  # qs <- exp(log_qs)

  for (f in seq_along(fleets)){

    for (ff in seq_along(fauna)){

      fleets[[f]]$metiers[[ff]]$catchability <- qs[cc]

      if (all(fleets[[f]]$metiers[[ff]]$spatial_catchability == 0)) {
        # annoying step: if q = 0 from earlier, then this will be a matrix of zeros and can't get updated
        fleets[[f]]$metiers[[ff]]$spatial_catchability <-
          rep(1, length(fleets[[f]]$metiers[[ff]]$spatial_catchability))
      }

      mean_q <- mean(fleets[[f]]$metiers[[ff]]$spatial_catchability)

      mean_q <- ifelse(mean_q == 0, 1e-9, mean_q)

      fleets[[f]]$metiers[[ff]]$spatial_catchability <- fleets[[f]]$metiers[[ff]]$spatial_catchability  / mean_q * qs[cc]

      cc <- cc + 1
    }

  }

  storage <- simmar(fauna = fauna,
                    fleets = fleets,
                    years = years)

  # tmp <- purrr::map_dfr(storage[[length(storage)]], ~as.data.frame(.x$ssb_p_a), .id = "fauna")

  tmp <- purrr::map(storage[[length(storage)]], ~as.data.frame(.x$ssb_p_a)) |>
    purrr::list_rbind(names_to = "fauna")


  b_p <- rowSums(tmp[,2:ncol(tmp)], na.rm = TRUE)

  tmp <- data.frame(fauna = tmp$fauna, ssb = b_p) %>%
    dplyr::group_by(fauna) %>%
    dplyr::summarise(ssb = sum(ssb)) %>%
    dplyr::arrange(fauna)

  ssb0s <- purrr::map_dbl(fauna, "ssb0")

  ssb0s <- ssb0s[sort(names(ssb0s))]

  target_depletion <- purrr::map_dbl(fauna, "fished_depletion")

  target_depletion <- target_depletion[sort(names(target_depletion))]

  tmp$depletion <- tmp$ssb / ssb0s

  ss <- sum((log(tmp$depletion) - log(target_depletion))^2)

  return(ss)


}

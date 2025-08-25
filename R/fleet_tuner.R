#' fleet tuner
#'
#' finds catchability (q) by fleet such that
#' target fished depletion is achieved
#'
#' @param fauna fauna object
#' @param years number of years to tune
#' @param fleets fleet object
#' @param log_fs log instantaneous fishing mortality per critter
#' @param e_fl baseline effort per fleet
#'
#' @return objective function of fleet tuner
#' @export
#'
fleet_tuner <- function(log_fs, fauna, fleets,e_fl, years = 50) {

  fs <- exp(log_fs)

  tfleets <- fleets

  for (i in length(tfleets)){


    for (j in 1:length(tfleets[[i]]$metiers)){

      tfleets[[i]]$metiers[[j]] <- fleets[[i]]$metiers[[j]]$clone(deep = TRUE)

    }

  }


  for (f in seq_along(tfleets)) {
    for (ff in seq_along(fauna)) {

      f_critter <- fs[ff]

      f_metier <-  tfleets[[f]]$metiers[[ff]]$p_explt * f_critter

      metier_q <- f_metier / e_fl[f]

      # print(e_fl)
      #
      # print(metier_q)

      tfleets[[f]]$metiers[[ff]]$catchability <- metier_q

      if (all(tfleets[[f]]$metiers[[ff]]$spatial_catchability == 0)) {
        # annoying step: if q = 0 from earlier, then this will be a matrix of zeros and can't get updated
        tfleets[[f]]$metiers[[ff]]$spatial_catchability <-
          rep(1, length(tfleets[[f]]$metiers[[ff]]$spatial_catchability))
      }

      mean_q <- mean(tfleets[[f]]$metiers[[ff]]$spatial_catchability)

      mean_q <- ifelse(mean_q == 0, 1e-9, mean_q)

      tfleets[[f]]$metiers[[ff]]$spatial_catchability <- tfleets[[f]]$metiers[[ff]]$spatial_catchability / mean_q * metier_q

    }
  }

  storage <- simmar(
    fauna = fauna,
    fleets = tfleets,
    years = years
  )

  # tmp <- purrr::map_dfr(storage[[length(storage)]], ~as.data.frame(.x$ssb_p_a), .id = "fauna")

  tmp <- purrr::map(storage[[length(storage)]], ~ as.data.frame(.x$ssb_p_a)) |>
    purrr::list_rbind(names_to = "fauna")


  b_p <- rowSums(tmp[, 2:ncol(tmp)], na.rm = TRUE)

  tmp <- data.frame(fauna = tmp$fauna, ssb = b_p) %>%
    dplyr::group_by(fauna) %>%
    dplyr::summarise(ssb = sum(ssb)) %>%
    dplyr::arrange(fauna)

  ssb0s <- purrr::map_dbl(fauna, "ssb0")

  ssb0s <- ssb0s[sort(names(ssb0s))]

  target_depletion <- purrr::map_dbl(fauna, "depletion")

  target_depletion <- target_depletion[sort(names(target_depletion))]

  tmp$depletion <- tmp$ssb / ssb0s

  ss <- sum((log(tmp$depletion) - log(target_depletion))^2)

  rm(storage)

  return(ss)
}

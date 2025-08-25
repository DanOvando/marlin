#' tune_fleets tunes parameters of the fleet model to achieve desired initial conditions.
#' note that this is not exact: post-tuning values will not perfectly match inputs since
#' for example some tuning steps depend on prior tuning step, making it difficult to tune everything
#' at once.
#'
#' @param fauna a fauna object
#' @param fleets a fleet object
#' @param years the number of years to tune over
#' @param tune_type one of 'f' or 'depletion' to tune catchability to achieve a desired fishing mortality rate (f) or a target depletion (B/B0)
#' @param tune_costs TRUE or FALSE to tune costs to a target cost to revenue ratio
#' @param fine_tune_costs TRUE or FALSE
#'
#' @return tuned fleet object
#' @export
tune_fleets <- function(fauna,
                        fleets,
                        years = 50,
                        tune_type = "f",
                        tune_costs = TRUE,
                        fine_tune_costs = TRUE) {

  tfleets <- fleets

  if (tune_type == "explt"){
    tune_type = "f"
  }

  for (i in length(tfleets)){


    for (j in 1:length(tfleets[[i]]$metiers)){

      tfleets[[i]]$metiers[[j]] <- fleets[[i]]$metiers[[j]]$clone(deep = TRUE)

    }

  }

  fleet_names <- names(tfleets)

  fauni <- names(fauna)

  og_fleet_model <- character(length = length(tfleets))
  # names(og_fleet_model) <- names(tfleets)

  for (f in names(tfleets)) {
    og_fleet_model[f] <- tfleets[[f]]$fleet_model
    tfleets[[f]]$fleet_model <- "constant_effort"
  }

  # normalize p_explt to make sure it sums to 1

  for (s in fauni) {
    p_explts <- purrr::map_dbl(tfleets, c("metiers", s, "p_explt"))

    p_explts <- p_explts / ifelse(sum(p_explts) > 0, sum(p_explts), 1e-6)

    for (f in fleet_names) {
      tfleets[[f]]$metiers[[s]]$p_explt <- as.numeric(p_explts[f])
    } # close fleet loop
  } # close fauna loop


  # find patches within fishing grounds where b0_p is > 0

  fauni <- names(fauna)

  fleeti <- names(tfleets)

  for (s in fauni) {

    e_fl <- rep(NA, length(tfleets))
    # names(e_fl) <- names(tfleets)
    # evenly divide up fishing effort into fishing grounds
    # find intersection with viable habitat for the species in question (b0_p > 0)
    # really, calculate proportion of fishing grounds that have viable habitat for species in question

    b0_p <- fauna[[s]]$b0_p

    for (fl in seq_along(fleeti)){

      e_fl[fl] <- (tfleets[[fl]]$base_effort / sum(tfleets[[fl]]$fishing_grounds$fishing_ground > 0)) * mean(b0_p[tfleets[[fl]]$fishing_grounds$fishing_ground > 0] > 0)


    }

    p_explt <-
      purrr::map_dbl(tfleets, c("metiers", s, "p_explt"))

    explt_by_fleet <- (fauna[[s]]$init_explt) * p_explt

    catchability <- explt_by_fleet / e_fl

    if (any(catchability > 1) & tune_type != "depletion"){

      stop("Desired exploitation rate if not possible given supplied effort levels (q >=1 ); try increasing base_effort")

    }

    for (f in fleeti) {

      tfleets[[f]]$metiers[[s]]$catchability <- catchability[f]

      if (all(tfleets[[f]]$metiers[[s]]$spatial_catchability == 0)) {
        # annoying step: if q = 0 from earlier, then this will be a matrix of zeros and can't get updated
        tfleets[[f]]$metiers[[s]]$spatial_catchability <-
          rep(1, length(tfleets[[f]]$metiers[[s]]$spatial_catchability))
      }

      mean_q <- mean(tfleets[[f]]$metiers[[s]]$spatial_catchability)

      mean_q <- ifelse(mean_q == 0, 1e-9, mean_q)

      tfleets[[f]]$metiers[[s]]$spatial_catchability <- (tfleets[[f]]$metiers[[s]]$spatial_catchability / mean_q) * catchability[f]
    } # close internal fleet loop
  } # close fauna loop

  if (tune_type == "depletion") {


    log_fs <-
      optim(
        par = rep(3, length(fauna)),
        fleet_tuner,
        fleets = tfleets,
        e_fl = e_fl,
        fauna = fauna,
        years = years,
        upper = rep(log(7), length(fauna)),
        method = "L-BFGS-B"
      )

    if (log_fs$value > .1){
      warning("tune_fleets failed to match desired depletion levels, check whether target depletion is plausible given supplied selectivities, fishing grounds, and p_explt.")
    }

    for (f in seq_along(tfleets)) {
      for (ff in seq_along(fauna)) {
        f_critter <- exp(log_fs$par[ff])

        f_metier <-  tfleets[[f]]$metiers[[ff]]$p_explt * f_critter

        metier_q <- f_metier / e_fl[f]

        tfleets[[f]]$metiers[[ff]]$catchability <- metier_q

        if (all(tfleets[[f]]$metiers[[ff]]$spatial_catchability == 0)) {
          # annoying step: if q = 0 from earlier, then this will be a matrix of zeros and can't get updated
          tfleets[[f]]$metiers[[ff]]$spatial_catchability <-
            rep(1, length(tfleets[[f]]$metiers[[ff]]$spatial_catchability))
        }

        mean_q <- mean(tfleets[[f]]$metiers[[ff]]$spatial_catchability)

        mean_q <- ifelse(mean_q == 0, 1e-9, mean_q)

        tfleets[[f]]$metiers[[ff]]$spatial_catchability <- tfleets[[f]]$metiers[[ff]]$spatial_catchability / mean_q * metier_q

      } # close internal fauna loop
    } # close fleet loop
  } # close depletion

  # storage <- simmar(
  #   fauna = fauna,
  #   tfleets = tfleets,
  #   years = years
  # )
  #
  # revenue <-
  #   purrr::map(
  #     storage[[1]],
  #     ~ data.frame(expand.grid(dimnames(.x$r_p_a_fl)), value = as.vector(.x$r_p_a_fl)) |>
  #       purrr::set_names("patch", "age", "fleet", "revenue") |>
  #       dplyr::mutate(across(patch:age, ~ as.numeric(as.character(
  #         .x
  #       ))))
  #   ) %>%
  #   dplyr::bind_rows(.id = "critter") %>%
  #   dplyr::group_by(fleet) %>%
  #   dplyr::summarise(revenue = sum(revenue, na.rm = TRUE))


  if (tune_costs) {
    init_sim <- simmar(
      fauna = fauna,
      fleets = tfleets,
      years = years
    )

    eq <- init_sim[[length(init_sim)]]


    revenue <-
      purrr::map(
        eq,
        ~ tibble::rownames_to_column(data.frame(revenue = (colSums(.x$r_p_fl, na.rm = TRUE))), "fleet")
      ) %>%
      purrr::list_rbind(names_to = "critter") |>
      dplyr::group_by(fleet) %>%
      dplyr::summarise(revenue = sum(revenue))

    effort <- purrr::map(
      eq[1],
      ~ data.frame(.x$e_p_fl) %>% dplyr::mutate(patch = 1:nrow(.))
    ) %>%
      purrr::list_rbind(names_to = "critter") |>
      tidyr::pivot_longer(-c(critter, patch), names_to = "fleet", values_to = "effort") # effort is the same for all critters per fleet so only selecting first entry


    effort <- purrr::map(
      eq[1],
      ~ data.frame(.x$e_p_fl) %>% dplyr::mutate(patch = 1:nrow(.))
    ) %>%
      purrr::list_rbind(names_to = "critter") |>
      tidyr::pivot_longer(-c(critter, patch), names_to = "fleet", values_to = "effort") # effort is the same for all critters per fleet so only selecting first entry



    cost_per_patch <-
      purrr::map(
        tfleets,
        ~ data.frame(
          patch = 1:length(.x$cost_per_patch),
          cost = .x$cost_per_patch
        )
      ) |>
      purrr::list_rbind(names_to = "fleet")

    base_cr_ratio <- purrr::map(tfleets, ~ data.frame(base_cr_ratio = .x$cr_ratio), .id = "fleet") |>
      purrr::list_rbind(names_to = "fleet")

    fleet_cost_expos <- purrr::map(tfleets, ~ data.frame(beta = .x$effort_cost_exponent)) |>
      purrr::list_rbind(names_to = "fleet")

    effort_and_costs <- effort %>%
      dplyr::left_join(cost_per_patch, by = c("fleet", "patch")) %>%
      dplyr::group_by(fleet) %>%
      dplyr::summarise(
        travel_cost = sum(effort * cost),
        effort = sum(effort)
      ) %>%
      dplyr::left_join(revenue, by = "fleet") %>%
      dplyr::left_join(base_cr_ratio, by = "fleet") %>%
      dplyr::left_join(fleet_cost_expos, by = "fleet") %>%
      dplyr::mutate(cost_per_unit_effort = (base_cr_ratio * revenue) / (effort^
                                                                          beta + travel_cost)) %>%
      dplyr::mutate(profits = revenue - (cost_per_unit_effort * (effort^
                                                                   beta + travel_cost)))

    for (f in names(tfleets)) {
      tfleets[[f]]$cost_per_unit_effort <- effort_and_costs$cost_per_unit_effort[effort_and_costs$fleet == f]

    }
  }

  for (f in names(tfleets)) {
    tfleets[[f]]$fleet_model <- og_fleet_model[f]
  }



  return(tfleets)
}

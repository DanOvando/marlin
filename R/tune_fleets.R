#' tune_fleets tunes parameters of the fleet model to achieve desired initial conditions.
#' note that this is not exact: post-tuning values will not perfectly match inputs since
#' for example some tuning steps depend on prior tuning step, making it difficult to tune everything
#' at once.
#'
#' @param fauna a fauna object
#' @param fleets a fleet object
#' @param years the number of years to tune over
#' @param tune_type one of 'explt' or 'depletion' to tune to an exploitation rate or a target depletion (B/B0)
#' @param tune_costs TRUE or FALSE to tune costs to a target cost to revenue ratio
#' @param fine_tune_costs TRUE or FALSE
#'
#' @return tuned fleet object
#' @export
tune_fleets <- function(fauna,
                        fleets,
                        years = 50,
                        tune_type = "explt",
                        tune_costs = TRUE,
                        fine_tune_costs = TRUE) {

  # might be best to define costs based on unfished rather than EQ to get around weird tuning behavior with tuning assuming revenue distribution but fleet model assuming IFD

  fleet_names <- names(fleets)

  fauni <- names(fauna)

  og_fleet_model <- character(length = length(fleets))
  names(og_fleet_model) <- names(fleets)

  for (f in names(fleets)){
    og_fleet_model[f] <- fleets[[f]]$fleet_model
    fleets[[f]]$fleet_model <- "constant_effort"
  }

  # normalize p_explt to make sure it sums to 1

  for (s in fauni) {
    p_explts <- purrr::map_dbl(fleets,c("metiers",s,"p_explt"))

    p_explts <- p_explts / ifelse(sum(p_explts) > 0, sum(p_explts), 1e-6)

    for ( f in fleet_names){

      fleets[[f]]$metiers[[s]]$p_explt <- as.numeric(p_explts[f])

    } # close fleet loop

  } # close fauna loop

  storage <- simmar(fauna = fauna,
                      fleets = fleets,
                      years = years)





  revenue <-
    purrr::map(
      storage[[1]],
      ~  data.frame(expand.grid(dimnames(.x$r_p_a_fl)), value = as.vector(.x$r_p_a_fl)) |>
        purrr::set_names("patch", "age", "fleet", "revenue") |>
        dplyr::mutate(across(patch:age, ~ as.numeric(as.character(
          .x
        ))))
    ) %>%
    dplyr::bind_rows(.id = "critter") %>%
    dplyr::group_by(fleet) %>%
    dplyr::summarise(revenue = sum(revenue, na.rm = TRUE))

    fauni <- names(fauna)

    fleeti <- names(fleets)

    for (s in fauni) {
      e_p_fl <- storage[[length(storage)]][[1]]$e_p_fl

      b_p = rowSums(storage[[length(storage)]][[s]]$b_p_a)

      weights <- b_p / max(b_p)

      e_fl <- colSums((e_p_fl * weights)) / sum(weights)  # calculate the total effort weighted by biomass of that species in patches.

      p_explt <-
        purrr::map_dbl(fleets, c("metiers",s, "p_explt"))[names(e_p_fl)]

      explt_by_fleet <- (fauna[[s]]$init_explt)  * p_explt

      catchability <-  explt_by_fleet / e_fl
      for (f in fleeti) {
        fleets[[f]]$metiers[[s]]$catchability <- catchability[f]

        if (all(fleets[[f]]$metiers[[s]]$spatial_catchability == 0)) {
          # annoying step: if q = 0 from earlier, then this will be a matrix of zeros and can't get updated
          fleets[[f]]$metiers[[s]]$spatial_catchability <-
            rep(1, length(fleets[[f]]$metiers[[s]]$spatial_catchability))
        }

        mean_q <- mean(fleets[[f]]$metiers[[s]]$spatial_catchability)

        mean_q <- ifelse(mean_q == 0, 1e-9, mean_q)

        fleets[[f]]$metiers[[s]]$spatial_catchability <-    (fleets[[f]]$metiers[[s]]$spatial_catchability  / mean_q) * catchability[f]

      } # close internal fleet loop

    } # close fauna loop



  if (tune_type == "depletion") {

    fleet_fauna <- length(fauni) * length(fleeti)

    qs <- vector(mode = "double", length = (fleet_fauna))

    cc <- 0
    for (f in fauni){

      for (ff in fleeti){

        cc <- cc + 1

        qs[cc] <- fleets[[ff]]$metiers[[f]]$catchability

      }

    }

    qs <-
      optim(
        par = qs,
        fleet_tuner,
        fleets = fleets,
        fauna = fauna,
        years = years,
        lower = rep(0, length(fauna) * length(fleets)),
        upper = rep(20, length(fauna) * length(fleets)),
        method = "L-BFGS-B"
      )
    cc <- 1

    for (f in seq_along(fleets)) {
      for (ff in seq_along(fauna)) {
        fleets[[f]]$metiers[[ff]]$catchability <- qs$par[cc]


        if (all(fleets[[f]]$metiers[[ff]]$spatial_catchability == 0)) {
          # annoying step: if q = 0 from earlier, then this will be a matrix of zeros and can't get updated
          fleets[[f]]$metiers[[ff]]$spatial_catchability <-
            rep(1, length(fleets[[f]]$metiers[[ff]]$spatial_catchability))
        }

        mean_q <- mean(fleets[[f]]$metiers[[ff]]$spatial_catchability)

        mean_q <- ifelse(mean_q == 0, 1e-9, mean_q)


        fleets[[f]]$metiers[[ff]]$spatial_catchability <-fleets[[f]]$metiers[[ff]]$spatial_catchability  / mean_q * qs$par[cc]
        cc <- cc + 1
      } # close internal fauna loop

    } # close fleet loop
  } # close depletion

    if (tune_costs){

      init_sim <- simmar(fauna = fauna,
                      fleets = fleets,
                      years = years)

    eq <- init_sim[[length(init_sim)]]


    revenue <-
      purrr::map(eq,
                     ~ tibble::rownames_to_column(data.frame( revenue = (colSums(.x$r_p_fl, na.rm = TRUE))), "fleet")) %>%
      purrr::list_rbind(names_to = "critter") |>
      dplyr::group_by(fleet) %>%
      dplyr::summarise(revenue = sum(revenue))

    effort <-    purrr::map(eq[1],
                                ~ data.frame(.x$e_p_fl) %>% dplyr::mutate(patch = 1:nrow(.))) %>%
      purrr::list_rbind(names_to = "critter") |>
      tidyr::pivot_longer(-c(critter, patch), names_to = "fleet", values_to = "effort") # effort is the same for all critters per fleet so only selecting first entry


    effort <-    purrr::map(eq[1],
                                ~ data.frame(.x$e_p_fl) %>% dplyr::mutate(patch = 1:nrow(.))) %>%
      purrr::list_rbind(names_to = "critter") |>
      tidyr::pivot_longer(-c(critter, patch), names_to = "fleet", values_to = "effort") # effort is the same for all critters per fleet so only selecting first entry



    cost_per_patch <-
      purrr::map(fleets,
                     ~ data.frame(
                       patch = 1:length(.x$cost_per_patch),
                       cost = .x$cost_per_patch
                     )) |>
      purrr::list_rbind(names_to = "fleet")

    base_cr_ratio <- purrr::map(fleets, ~data.frame(base_cr_ratio = .x$cr_ratio), .id = "fleet") |>
      purrr::list_rbind(names_to = "fleet")

    fleet_cost_expos <- purrr::map(fleets, ~ data.frame(beta = .x$effort_cost_exponent)) |>
      purrr::list_rbind(names_to = "fleet")

    effort_and_costs <- effort %>%
      dplyr::left_join(cost_per_patch, by = c("fleet", "patch")) %>%
      dplyr::group_by(fleet) %>%
      dplyr::summarise(travel_cost = sum(effort * cost),
                       effort = sum(effort)) %>%
      dplyr::left_join(revenue, by = "fleet") %>%
      dplyr::left_join(base_cr_ratio, by = "fleet") %>%
      dplyr::left_join(fleet_cost_expos, by = "fleet") %>%
      dplyr::mutate(cost_per_unit_effort = (base_cr_ratio * revenue) / (effort ^
                                                                          beta + travel_cost)) %>%
      dplyr::mutate(profits =  revenue - (cost_per_unit_effort * (effort ^
                                                                    beta + travel_cost)))

    for (f in names(fleets)){
      fleets[[f]]$cost_per_unit_effort <- effort_and_costs$cost_per_unit_effort[effort_and_costs$fleet == f]

      fleets[[f]]$fleet_model <-    og_fleet_model[f]

    }

    }



  return(fleets)
}

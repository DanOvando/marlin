#' Internal function to fine-tune cost parameters
#'
#' @param log_cost_mult log cost per unit effort multiplier
#' @param target the target to be achieved
#' @param fauna a list of fauna
#' @param fleets a list of fleets
#' @param years the number of years to run the optimization
#' @param tune_type what is being tuned for
#'
#' @return fine-tuned costs
#' @export
#'
fine_tune_costs <-
  function(log_cost_mult,
           target,
           fauna,
           fleets,
           years = 25,
           tune_type) {
    fleet_models <- purrr::map_chr(fleets, "fleet_model")

    oa_fleets <- names(fleet_models)[fleet_models == "open_access"]

    tmp_fleets <- fleets

    for (f in seq_along(oa_fleets)) {
      tmp_fleets[[oa_fleets[f]]]$cost_per_unit_effort <-
        exp(log_cost_mult[f])
    }

    sim <- simmar(
      fauna = fauna,
      fleets = tmp_fleets,
      years = years
    )

    if (tune_type == "depletion") {
      depletion <-
        purrr::map_df(sim[[length(sim)]], ~ sum(.x$ssb_p_a) / .x$ssb0) %>% # depletion of each species
        tidyr::pivot_longer(tidyselect::everything(),
          names_to = "critter",
          values_to = "depletion"
        )

      ss_frame <- depletion %>%
        dplyr::left_join(target, by = "critter")

      out <- sum((ss_frame$depletion - ss_frame$target)^2)
    } else if (tune_type == "explt") {
      for (s in names(fauna)) {
        e_p_fl <- sim[[length(sim)]][[1]]$e_p_fl

        b_p <- rowSums(sim[[length(sim)]][[s]]$b_p_a)

        weights <- b_p / max(b_p)

        e_fl <-
          colSums((e_p_fl * weights)) / sum(weights) # calculate the total effort weighted by biomass of that species in patches.

        p_explt <-
          purrr::map_dbl(fleets, c("metiers", s, "p_explt"))[names(e_p_fl)]

        explt_by_fleet <- (fauna[[s]]$init_explt) * p_explt

        # catchability <-  log(1 - explt_by_fleet) / -e_fl

        catchability <- explt_by_fleet / e_fl
      }
    }


    return(out)
  }

# target <- data.frame(critter = "bigeye", target = 0.2)
#
# a <- Sys.time()
# test <-
#   nlminb(
#     log(0.001),
#     fine_tune_costs,
#     fauna = fauna,
#     fleets = fleets,
#     target = target,
#     tune_type = "depletion",
#     upper = log(10)
#   )
# Sys.time() - a
#
# test_fleets <- fleets
#
# fleet_models <- purrr::map_chr(test_fleets, "fleet_model")
#
# oa_fleets <- names(fleet_models)[fleet_models == "open_access"]
#
# for (f in seq_along(oa_fleets)) {
#   test_fleets[[oa_fleets[f]]]$cost_per_unit_effort <- exp(test$par[f])
#
# }
#
# test_sim <- simmar(fauna = fauna,
#                    fleets = test_fleets)
#
# pt <- process_marlin(test_sim, keep_age = FALSE)
#
#
# plot_marlin(pt)
#
#
# effort <-
#   map_df(test_sim, ~ data.frame(effort = sum(.x$bigeye$e_p_fl)), .id = "step") %>%
#   mutate(step = as.numeric(step))
#
# effort %>%
#   ungroup() %>%
#   ggplot(aes(step, effort)) +
#   geom_line()

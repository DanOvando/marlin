#' Optimize MPA network
#'
#' @param fauna fauna object
#' @param fleets fleet object
#' @param alpha the weight given to conservation relative to economics (range 0 to 1)
#' @param max_prop_mpa max proportion of cells to place in MPA
#' @param resolution the resolution of the simulated system
#' @param prop_sampled the proportion of cells to sample in each SIR iteration
#' @param workers number of workers for parallel process
#' @param starting_conditions starting conditions of the simmar object
#' @param objective one of max_ssb to maximize spawning stock biomass or min_loss to prioritize not having any losses
#'
#' @return a list with results of MPA optimization
#' @export
#'
optimize_mpa <-
  function(fauna,
           fleets,
           starting_conditions = NA,
           alpha = 0.33,
           max_prop_mpa = 1,
           resolution,
           prop_sampled = .2,
           max_delta = 2,
           workers = 6,
           bio_objective = "max_ssb",
           econ_objective = "yield",
           work_backwards = TRUE,
           patches_at_a_time = 1) {
    if (length(resolution) == 1) {
      resolution <- rep(resolution, 2)
    }
    future::plan(future::multisession, workers = workers)

    on.exit(future::plan(future::sequential))

    patches <- prod(resolution)

    samps <- round(prop_sampled * patches)

    mpas <-
      tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) %>%
      dplyr::mutate(patch = 1:nrow(.))

    if (max_prop_mpa < 1) {
      work_backwards <- FALSE
      message("Setting work_backwards to FALSE since max_prop_mpa < 1")
    }

    # if work_backwards, start at 100% coverage and remove MPAs. Otherwise, start at 0% and grow. Working backwards produces more stable results
    if (work_backwards) {
      mpas$mpa <- TRUE
    } else {
      mpas$mpa <- FALSE
    }

    max_patches_protected <- round(patches * max_prop_mpa)

    candidate_patches <- 1:patches

    is <- seq(1, max_patches_protected, by = patches_at_a_time)

    patch_value <- dplyr::tibble(patch = 1:patches, obj_value = 0.5)

    results <- vector(mode = "list", length = length(is))

    mpa_network <-
      vector(mode = "list", length = length(is))

    # set up open access

    calc_objective_function <-
      function(candidate_patch,
               fauna,
               fleets,
               mpas,
               starting_conditions,
               econ_objective = "yield",
               bio_objective = "max_ssb",
               work_backwards = TRUE) {
        tmp_mpas <- mpas

        if (work_backwards) {
          tmp_mpas$mpa[tmp_mpas$patch %in% candidate_patch] <- FALSE
        } else {
          tmp_mpas$mpa[tmp_mpas$patch %in% candidate_patch] <- TRUE
        }

        sim_mpa <- simmar(
          fauna = fauna,
          fleets = fleets,
          years = years,
          manager = list(mpas = list(
            locations = tmp_mpas,
            mpa_year = 1
          )),
          initial_conditions = starting_conditions
        )


        res <-
          sim_mpa[[length(sim_mpa)]] # for now, just calculate in the final timestep ## make this in final year

        biodiv_mpa <-
          (purrr::map_dbl(res, ~ sum(.x$ssb_p_a) / .x$ssb0)) # calculate biodiversity component of objective function

        biodiv_sq <-
          purrr::map_dbl(starting_conditions, ~ sum(.x$ssb_p_a) / .x$ssb0) # calculate biodiversity component of objective function


        delta_biodiv <- biodiv_mpa - biodiv_sq

        if (bio_objective == "max_ssb") {
          biodiv <-
            sum(purrr::map_dbl(res, ~ sum(.x$ssb_p_a) / .x$ssb0)) # calculate biodiversity component of objective function
        } else if (bio_objective == "min_loss") {
          biodiv <- sum((biodiv_mpa - biodiv_sq) >= 0)
        }

        revenue_mpa <-
          (purrr::map_dbl(res, ~ sum(.x$r_p_a_fl, na.rm = TRUE)))

        profit_mpa <-
          (purrr::map_dbl(res, ~ sum(.x$prof_p_fl, na.rm = TRUE)))

        yield_mpa <-
          (purrr::map_dbl(res, ~ sum(.x$c_p_a, na.rm = TRUE)))


        if (econ_objective == "profits") {
          econ <- sum(profit_mpa, na.rm = TRUE)
        } else if (econ_objective == "revenue") {
          econ <- sum(revenue_mpa, na.rm = TRUE)
        } else if (econ_objective == "yield") {
          econ <- sum(yield_mpa, na.rm = TRUE)
        }

        out <- dplyr::tibble(biodiv = biodiv, econ = econ)
      } # close objective function


    for (i in seq_along(is)) {
      # determine marginal objective function value of each sampled cell
      #
      marginal_values <-
        furrr::future_map_dfr(
          candidate_patches,
          calc_objective_function,
          bio_objective = bio_objective,
          econ_objective = econ_objective,
          fauna = fauna,
          fleets = fleets,
          mpas = mpas,
          starting_conditions = starting_conditions,
          .options = furrr::furrr_options(seed = 42),
          .progress = FALSE
        )

      marginal_values$patch <- candidate_patches # assign patches

      marginal_values$obj_value <-
        alpha * scales::rescale(marginal_values$biodiv) + (1 - alpha) * scales::rescale(marginal_values$econ) # calculate objective function. Rescaling means that alpha dictates the weithing of a given percent rank of biodiversity relative to a relative percent rank of economics

      marginal_values <- marginal_values %>%
        dplyr::arrange(patch) # make sure marginal values are ordered by patches


      top_patches <- marginal_values %>%
        dplyr::arrange(dplyr::desc(obj_value))

      top_patches <- top_patches$patch[1:patches_at_a_time]

      patch_value$obj_value[patch_value$patch %in% marginal_values$patch] <-
        marginal_values$obj_value

      # update MPA locations

      if (work_backwards) {
        mpas$mpa[mpas$patch %in% top_patches] <- FALSE
      } else {
        mpas$mpa[mpas$patch %in% top_patches] <- TRUE
      }


      mpa_network[[i]] <- mpas

      # store mpa network results

      tmp_result <- simmar(
        fauna = fauna,
        fleets = fleets,
        years = years,
        manager = list(mpas = list(
          locations = mpas,
          mpa_year = 1
        )),
        initial_conditions = starting_conditions
      )

      res <-
        tmp_result[[length(tmp_result)]] # for now, just calculate in the final timestep

      biodiv <-
        (purrr::map_dbl(res, ~ sum(.x$ssb_p_a) / .x$ssb0)) # calculate biodiversity component of objective function

      biodiv_sq <-
        purrr::map_dbl(starting_conditions, ~ sum(.x$ssb_p_a) / .x$ssb0) # calculate biodiversity component of objective function


      delta_biodiv <- biodiv - biodiv_sq


      revenue_mpa <-
        (purrr::map_dbl(res, ~ sum(.x$r_p_a_fl, na.rm = TRUE)))

      profit_mpa <-
        (purrr::map_dbl(res, ~ sum(.x$prof_p_fl, na.rm = TRUE)))

      yield_mpa <-
        (purrr::map_dbl(res, ~ sum(.x$c_p_a, na.rm = TRUE)))


      if (econ_objective == "profits") {
        econ <- profit_mpa
      } else if (econ_objective == "revenue") {
        econ <- revenue_mpa
      } else if (econ_objective == "yield") {
        econ <- yield_mpa
      }



      out <-
        list(
          biodiv = dplyr::tibble(critter = names(biodiv), ssb_v_ssb0 = biodiv),
          econ = dplyr::tibble(critter = names(econ), econ = econ),
          yield = dplyr::tibble(critter = names(yield_mpa), yield = yield_mpa)
        )
      # keep just biomass and catch for now: you can always use the MPA layer at that step to recreate the whole sim if you need

      results[[i]] <- out

      # update candidate cells
      #
      if (work_backwards) {
        if (sum(mpas$mpa) > 1) {
          candidate_patches <-
            sample(patch_value$patch[mpas$mpa],
              ceiling(sum(mpas$mpa) * prop_sampled),
              prob = patch_value$obj_value[mpas$mpa] + 1e-3
            )
        } else {
          candidate_patches <- patch_value$patch[mpas$mpa]
        }
      } else {
        candidate_patches <-
          sample(patch_value$patch[!mpas$mpa],
            ceiling(sum(!mpas$mpa) * prop_sampled),
            prob = patch_value$obj_value[!mpas$mpa] + 1e-3
          )
      }
      message(
        glue::glue(
          "{scales::percent(is[i] / max_patches_protected)} done designing MPA network"
        )
      )
    } # close MPA size loop

    # run counterfactual experiment (world with no MPA)


    tmp_result <- simmar(
      fauna = fauna,
      fleets = fleets,
      years = years,
      initial_conditions = starting_conditions
    )

    res <-
      tmp_result[[length(tmp_result)]] # for now, just calculate in the final timestep

    biodiv <-
      (purrr::map_dbl(res, ~ sum(.x$ssb_p_a) / .x$ssb0)) # calculate biodiversity component of objective function

    biodiv_sq <-
      purrr::map_dbl(starting_conditions, ~ sum(.x$ssb_p_a) / .x$ssb0) # calculate biodiversity component of objective function


    delta_biodiv <- biodiv - biodiv_sq


    revenue_sq <-
      (purrr::map_dbl(res, ~ sum(.x$r_p_a_fl, na.rm = TRUE)))

    profit_sq <-
      (purrr::map_dbl(res, ~ sum(.x$prof_p_fl, na.rm = TRUE)))

    yield_sq <-
      (purrr::map_dbl(res, ~ sum(.x$c_p_a, na.rm = TRUE)))


    if (econ_objective == "profits") {
      econ_sq <- profit_sq
    } else if (econ_objective == "revenue") {
      econ_sq <- revenue_sq
    } else if (econ_objective == "yield") {
      econ_sq <- yield_sq
    }



    out <-
      list(
        biodiv = dplyr::tibble(critter = names(biodiv), ssb_v_ssb0 = biodiv),
        econ = dplyr::tibble(critter = names(econ_sq), econ = econ_sq),
        yield = dplyr::tibble(critter = names(yield_sq), yield = yield_sq)
      )

    baseline_outcome <-
      dplyr::tibble(p_protected = 0, out = list(out))


    if (work_backwards) {
      patch_order <- rev(is)
    } else {
      patch_order <- is
    }

    # process results
    outcomes <- baseline_outcome %>%
      dplyr::bind_rows(dplyr::tibble(
        p_protected = patch_order,
        out = results
      )) %>%
      dplyr::mutate(
        bio = purrr::map(out, "biodiv"),
        econ = purrr::map(out, ~ .x$econ %>% dplyr::select(-critter)),
        yield = purrr::map(out, ~ .x$yield %>% dplyr::select(-critter))
      ) %>%
      tidyr::unnest(cols = c(bio, econ, yield), names_repair = "universal")

    objective <- outcomes %>%
      dplyr::group_by(p_protected) %>%
      dplyr::summarise(
        ssb = sum(ssb_v_ssb0),
        econ = sum(econ),
        yield = sum(yield)
      ) %>%
      dplyr::mutate(obj = alpha * scales::rescale(ssb) + (1 - alpha) * (scales::rescale(econ)))

    out_mpa_network <-
      dplyr::tibble(
        p_protected = patch_order,
        mpa = mpa_network
      ) %>%
      tidyr::unnest(cols = mpa)

    return(list(
      outcomes = outcomes,
      objective = objective,
      mpa_network = out_mpa_network
    ))
  } # close optimize_mpa

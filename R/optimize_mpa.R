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
           objective = "max_ssb") {
    future::plan(future::multisession, workers = workers)
    
    on.exit(future::plan(future::sequential))
    
    # fauna <- casestudy$fauna[[1]]
    #
    # fleets <- casestudy$fleet[[1]]
    #
    # alpha <- 0
    #
    # percent_mpa <- 0.3
    #
    # prop_sampled <- 0.2
    #
    # max_prop_mpa <-  1
    
    patches <- resolution ^ 2
    
    samps <- round(prop_sampled * patches)
    
    mpas <- tidyr::expand_grid(x = 1:resolution, y = 1:resolution) %>%
      dplyr::mutate(patch = 1:nrow(.)) %>%
      dplyr::mutate(mpa = TRUE)
    
    max_patches_protected <- round(patches * max_prop_mpa)
    
    candidate_patches <- 1:patches
    
    patch_value <- dplyr::tibble(patch = 1:patches, obj_value = 0.5)
    
    results <- vector(mode = "list", length = max_patches_protected)
    
    mpa_network <-
      vector(mode = "list", length = max_patches_protected)
    
    # set up open access
    
    calc_objective_function <-
      function(candidate_patch, fauna, fleets, mpas,starting_conditions) {
        tmp_mpas <- mpas
        
        tmp_mpas$mpa[tmp_mpas$patch %in% candidate_patch] <- FALSE
        
        sim_mpa <- simmar(
          fauna = fauna,
          fleets = fleets,
          years = years,
          manager = list(mpas = list(locations = tmp_mpas,
                      mpa_year = 1)),
          initial_conditions = starting_conditions
        )
        
        
        res <-
          sim_mpa[[length(sim_mpa)]] # for now, just calculate in the final timestep
        
        biodiv_mpa <-
          (purrr::map_dbl(res, ~ sum(.x$ssb_p_a) / .x$ssb0)) # calculate biodiversity component of objective function
        
        biodiv_sq <-
          purrr::map_dbl(starting_conditions, ~ sum(.x$ssb_p_a) / .x$ssb0) # calculate biodiversity component of objective function
        
        
        delta_biodiv <- biodiv_mpa - biodiv_sq
        
        if (objective == "max_ssb"){
          biodiv <-
            sum(purrr::map_dbl(res, ~ sum(.x$ssb_p_a) / .x$ssb0)) # calculate biodiversity component of objective function
        } else if (objective == "min_loss") {
          
          
          biodiv <- sum((biodiv_mpa - biodiv_sq) >= 0)
          
        }
        
        econ_sq <-
          (purrr::map_dbl(starting_conditions, ~ sum(.x$r_p_a_fl, na.rm = TRUE)))
        
        profit_mpa <-  (purrr::map_dbl(res, ~ sum(.x$prof_p_fl, na.rm = TRUE)))
        
        econ <- sum(profit_mpa, na.rm = TRUE)
        
        out <- dplyr::tibble(biodiv = biodiv, econ = econ)
        
      }
    
    for (i in 1:max_patches_protected) {
      # determine marginal objective function value of each sampled cell
      #

      marginal_values <-
        furrr::future_map_dfr(
          candidate_patches,
          calc_objective_function,
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
      
      top_patch <-
        marginal_values$patch[marginal_values$obj_value == max(marginal_values$obj_value)][1] # find the next best patch to add
      # update marginal benefit surface
      
      patch_value$obj_value[patch_value$patch %in% marginal_values$patch] <-
        marginal_values$obj_value
      
      # update MPA locations
      
      mpas$mpa[mpas$patch == top_patch] <- FALSE
      
      mpa_network[[i]] <- mpas
      
      # store mpa network results
      
      tmp_result <- simmar(
        fauna = fauna,
        fleets = fleets,
        years = years,
        manager = list(mpas = list(locations = mpas,
                    mpa_year = 1)),
        initial_conditions = starting_conditions
      )
      
      res <-
        tmp_result[[length(tmp_result)]] # for now, just calculate in the final timestep
      
      biodiv <-
        (purrr::map_dbl(res, ~ sum(.x$ssb_p_a) / .x$ssb0)) # calculate biodiversity component of objective function
      
      biodiv_sq <-
        purrr::map_dbl(starting_conditions, ~ sum(.x$ssb_p_a) / .x$ssb0) # calculate biodiversity component of objective function
      
      
      delta_biodiv <- biodiv - biodiv_sq
      
      econ_sq <-
        (purrr::map_dbl(starting_conditions, ~ sum(.x$r_p_a_fl, na.rm = TRUE)))
      
      profit_mpa <-  (purrr::map_dbl(res, ~ sum(.x$prof_p_fl, na.rm = TRUE)))
      
      econ <- profit_mpa
      
      yield_mpa <-  (purrr::map_dbl(res, ~ sum(.x$c_p_fl, na.rm = TRUE)))
      
      
      out <-
        list(
          biodiv = dplyr::tibble(critter = names(biodiv), ssb_v_ssb0 = biodiv),
          econ = dplyr::tibble(critter = names(econ), econ = econ),
          yield = dplyr::tibble(critter = names(yield_mpa), yield = yield_mpa)
        )
      # keep just biomass and catch for now: you can always use the MPA layer at that step to recreate the whole sim if you need
      
      results[[i]] <- out
      
      # update candidate cells
      
      if (sum(mpas$mpa) > 1) {
        candidate_patches <-
          sample(patch_value$patch[mpas$mpa],
                 ceiling(sum(mpas$mpa) * prop_sampled),
                 prob = patch_value$obj_value[mpas$mpa] + 1e-3)
        
      } else {
        candidate_patches <- patch_value$patch[mpas$mpa]
      }
      
      message(glue::glue("{scales::percent(i / max_patches_protected)} done"))
    } # close MPA size loop
    
    # run counterfactual experiment (world with no MPA)
    
    
    tmp_result <- simmar(fauna = fauna,
                         fleets = fleets,
                         years = years,
                         initial_conditions = starting_conditions)
    
    res <-
      tmp_result[[length(tmp_result)]] # for now, just calculate in the final timestep
    
    biodiv <-
      (purrr::map_dbl(res, ~ sum(.x$ssb_p_a) / .x$ssb0)) # calculate biodiversity component of objective function
    
    biodiv_sq <-
      purrr::map_dbl(starting_conditions, ~ sum(.x$ssb_p_a) / .x$ssb0) # calculate biodiversity component of objective function
    
    
    delta_biodiv <- biodiv - biodiv_sq
    
    econ_sq <-
      (purrr::map_dbl(starting_conditions, ~ sum(.x$r_p_a_fl, na.rm = TRUE)))
    
    profit_sq <-  (purrr::map_dbl(res, ~ sum(.x$prof_p_fl, na.rm = TRUE)))
    
    econ <- profit_sq
    
    yield <-  (purrr::map_dbl(res, ~ sum(.x$c_p_fl, na.rm = TRUE)))
    
    out <-
      list(
        biodiv = dplyr::tibble(critter = names(biodiv), ssb_v_ssb0 = biodiv),
        econ = dplyr::tibble(critter = names(econ), econ = econ),
        yield = dplyr::tibble(critter = names(yield), yield = yield)
        
      )
    
    baseline_outcome <- dplyr::tibble(p_protected = 0, out = list(out))
    
    # process results
    
    outcomes <- baseline_outcome %>%
      dplyr::bind_rows(dplyr::tibble(
        p_protected = max_patches_protected:1,
        out = results
      )) %>%
      dplyr::mutate(bio = purrr::map(out, "biodiv"),
             econ = purrr::map(out, "econ"),
             yield = purrr::map(out, "yield")) %>%
      tidyr::unnest(cols = c(bio, econ, yield), names_repair = "universal") %>%
      dplyr::select(-`critter...5`,-`critter...7`) %>%
      dplyr::rename(critter = `critter...3`)
    
    objective <-  outcomes %>%
      dplyr::group_by(p_protected) %>%
      dplyr::summarise(ssb = sum(ssb_v_ssb0),
                econ = sum(econ),
                yield = sum(yield)) %>%
      dplyr:: mutate(obj = alpha * scales::rescale(ssb) + (1 - alpha) * (scales::rescale(econ)))
    
    out_mpa_network <-
      dplyr::tibble(p_protected = max_patches_protected:1,
             mpa = mpa_network) %>%
      tidyr::unnest(cols = mpa)
    
    return(list(
      outcomes = outcomes,
      objective = objective,
      mpa_network = out_mpa_network
    ))
  } # close optimize_mpa
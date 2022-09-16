#' Create Fleet
#' 
#' Creates a fleet object, mostly by adding in 
#' selectivity at age for each fleet and species
#'
#' @param base_effort base effort for the fleet
#' @param mpa_response one of "stay" or "leave" indicating response of vessels that used to fish in MPA to MPA
#' @param cr_ratio cost to revenue ratio at initial conditions (1 implies OA equilibrium, total profits = 0)
#' @param spatial_allocation spatial effort allocation strategy (ideal_free or revenue)
#' @param metiers a list of metiers
#' @param fleet_model which fleet model to use, one of "constant effort" or "open access"
#' @param profit_sensitivity the profit sensitivity of the open access model
#' @param cost_per_unit_effortt the cost per unit effort in the open access model
#'
#' @return a fleet object
#' @export
#'
create_fleet <-
  function(metiers,
           base_effort = 0,
           mpa_response = "stay",
           fleet_model = "constant effort",
           responsiveness = 0.5,
           cost_per_unit_effort = NA,
           spatial_allocation = "rpue",
           effort_cost_exponent = 1.3,
           ports = NULL,
           cost_per_distance = 1,
           cr_ratio = 1,
           resolution) {
    # idea: each fleet has a list of fauna inside of it specifying the price, selectivity, q for that species
    
    
    if (is.null(ports)){
      
      cost_per_patch <- rep(0, resolution^2)
      
    } else {
      
      patches <- tidyr::expand_grid(x = 1:resolution, y = 1:resolution) %>% 
        dplyr::mutate(patch = 1:nrow(.)) # extra step to make sure patch ordering is consistent
      
      ports <- ports %>% 
        dplyr::left_join(patches, by = c("x","y"))
      
      # calculate the distance between each of the patches
      port_distance <- distance <-
        tidyr::expand_grid(x = 1:resolution, y = 1:resolution) %>%
        dist(diag = TRUE) %>%
        as.matrix()
      
      
      port_distances <- apply(port_distance[ports$patch,],2,min) # calculate the distance from each patch to the port patches, then find the minimum distance
      
      cost_per_patch <- port_distances * cost_per_distance # calculate total travel cost per patch
      
      
    }
    
   fleet <- list(metiers = metiers, 
                 base_effort = base_effort, 
                 mpa_response = mpa_response,
                 cr_ratio = cr_ratio,
                 cost_per_unit_effort = cost_per_unit_effort,
                 responsiveness = responsiveness,
                 spatial_allocation = spatial_allocation,
                 fleet_model = fleet_model,
                 effort_cost_exponent = effort_cost_exponent,
                 cost_per_patch = cost_per_patch)

    
    return(fleet)
    
  } # close function
#' Create Fleet
#' 
#' Creates a fleet object, mostly by adding in 
#' selectivity at age for each fleet and species
#'
#' @param fleets 
#' @param fauna 
#' @param base_effort base effort for the fleet
#' @param mpa_response one of "stay" or "leave" indicating response of vessels that used to fish in MPA to MPA
#' @param cr_ratio cost to revenue ratio at initial conditions (1 implies OA equilibrium, total profits = 0)
#' @param spatial_allocation spatial effort allocation strategy (ideal_free or revenue)
#'
#' @return a fleet object
#' @export
#'
create_fleet <-
  function(metiers,
           base_effort = 0,
           mpa_response = "stay",
           spatial_allocation = "revenue",
           cr_ratio = 0.9) {
    # idea: each fleet has a list of fauna inside of it specifying the price, selectivity, q for that species
    
   fleet <- list(metiers = metiers, 
                 base_effort = base_effort, 
                 mpa_response = mpa_response,
                 cr_ratio = cr_ratio,
                 cost_per_unit_effort = NA,
                 spatial_allocation = spatial_allocation)

    
    return(fleet)
    
  } # close function
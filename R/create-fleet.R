#' Create Fleet
#' 
#' Creates a fleet object, mostly by adding in 
#' selectivity at age for each fleet and species
#'
#' @param fleets 
#' @param fauna 
#' @param base_effort 
#'
#' @return
#' @export
#'
create_fleet <-
  function(metiers,
           base_effort = 0,
           mpa_response = "stay") {
    # idea: each fleet has a list of fauna inside of it specifying the price, selectivity, q for that species
    
   fleet <- list(metiers = metiers, 
                 base_effort = base_effort, 
                 mpa_response = mpa_response)

    
    return(fleet)
    
  } # close function
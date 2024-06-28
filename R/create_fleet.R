#' Create Fleet
#'
#' Creates a fleet object, mostly by adding in
#' selectivity at age for each fleet and species
#'
#' @param base_effort base effort for the fleet
#' @param mpa_response one of "stay" or "leave" indicating response of vessels that used to fish in MPA to MPA
#' @param cr_ratio cost to revenue ratio at initial conditions (1 implies OA equilibrium, total profits = 0)
#' @param spatial_allocation spatial effort allocation strategy ('revenue','rpue','profit','ppue')
#' @param metiers a list of metiers
#' @param fleet_model which fleet model to use, one of "constant_effort" or "open_access"
#' @param profit_sensitivity the profit sensitivity of the open access model
#' @param cost_per_unit_effortt the cost per unit effort in the open access model
#'
#' @return a fleet object
#' @export
#'
create_fleet <-
  function(metiers,
           mpa_response = "stay",
           fleet_model = "constant_effort",
           responsiveness = 0.5,
           cost_per_unit_effort = 1,
           spatial_allocation = "rpue",
           effort_cost_exponent = 1,
           ports = NULL,
           cost_per_distance = 1,
           cr_ratio = 1,
           resolution,
           base_effort = NULL,
           fishing_grounds = NULL) {

    fleet_model <- gsub(" ","_", fleet_model) # in case someone used spaces accidentally (like dumbass old dan)

    if (length(resolution) == 1){
      resolution <- rep(resolution,2)
    }
    if (is.null(base_effort)){
      base_effort <- prod(resolution)
    }

    if (is.null(ports)){

      cost_per_patch <- rep(0, prod(resolution))

    } else {
      patches <- tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) %>%
        dplyr::mutate(patch = 1:nrow(.)) # extra step to make sure patch ordering is consistent

      if (any(ports$x > resolution[1]) | any(ports$y > resolution[2])){
        stop("one or more port locations is outside of spatial grid")
      }

      ports <- ports %>%
        dplyr::left_join(patches, by = c("x","y"))

      # calculate the distance between each of the patches
      port_distance <- distance <-
        tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) %>%
        dist(diag = TRUE) %>%
        as.matrix()

      port_distances <- apply(matrix(port_distance[ports$patch,], nrow = length(ports$patch)),2,min) # calculate the distance from each patch to the port patches, then find the minimum distance

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
                 cost_per_patch = cost_per_patch,
                 fishing_grounds = fishing_grounds)


    return(fleet)

  } # close function

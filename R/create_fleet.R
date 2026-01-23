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
#' @param responsiveness how responsive the fleet is to profits when fleet_model = "open_access"
#' @param cost_per_unit_effort the cost per unit effort
#' @param effort_cost_exponent exponent of costs
#' @param ports location of fishing ports
#' @param cost_per_distance cost per unit distance
#' @param resolution spatial resolution of the simulated seascape
#' @param patch_area the area of each patch (KM^2^)
#' @param fishing_grounds the location of fishing grounds (TRUE or FALSE)
#' @param fleet_model which fleet model to use, one of "constant_effort" or "open_access" or constant catch
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
           patch_area = 1,
           base_effort = NULL,
           fishing_grounds = NULL) {
    fleet_model <- stringr::str_replace_all(fleet_model, " ", "_") # in case someone used spaces accidentally (like dumbass old dan)

    if (length(resolution) == 1) {
      resolution <- rep(resolution, 2)
    }
    if (is.null(base_effort)) {
      base_effort <- prod(resolution)
    }

    if (is.null(fishing_grounds)){

      fishing_grounds <- tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) |>
        dplyr::mutate(fishing_ground = TRUE) |>
        dplyr::arrange(x,y)

    }

    fishing_grounds <- fishing_grounds |>
      dplyr::arrange(x,y)

    if (nrow(fishing_grounds) != prod(resolution)){
      stop("supplied fishing_grounds do not match the spatial dimensions of the simualted domain. Make sure that number of rows and columns match supplied resolution.")
    }

    if (is.null(ports)) {
      cost_per_patch <- rep(0, prod(resolution))
    } else {
      patches <- tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) %>%
        dplyr::mutate(patch = 1:nrow(.)) # extra step to make sure patch ordering is consistent

      if (any(ports$x > resolution[1]) | any(ports$y > resolution[2])) {
        stop("one or more port locations is outside of spatial grid")
      }

      ports <- ports %>%
        dplyr::left_join(patches, by = c("x", "y"))

      # calculate the distance between each of the patches
      port_distance <- distance <-
        tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) %>%
        dist(diag = TRUE) %>%
        as.matrix() * sqrt(patch_area)

      port_distances <- apply(matrix(port_distance[ports$patch, ], nrow = length(ports$patch)), 2, min) # calculate the distance from each patch to the port patches, then find the minimum distance

      cost_per_patch <- port_distances * cost_per_distance # calculate total travel cost per patch
    }

    fleet <- list(
      metiers = metiers,
      base_effort = base_effort,
      mpa_response = mpa_response,
      cr_ratio = cr_ratio,
      cost_per_unit_effort = cost_per_unit_effort,
      responsiveness = responsiveness,
      spatial_allocation = spatial_allocation,
      fleet_model = fleet_model,
      effort_cost_exponent = effort_cost_exponent,
      cost_per_patch = cost_per_patch,
      fishing_grounds = fishing_grounds
    )


    return(fleet)
  } # close function

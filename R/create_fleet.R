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
#' @param cost_per_unit_effort the cost per unit effort (deprecated - will be calibrated by tune_fleets)
#' @param effort_cost_exponent exponent of effort costs (gamma), controls congestion/convexity (default 1.2)
#' @param travel_fraction fraction of total costs from travel at equilibrium (default 0). Controls relative importance of spatial cost heterogeneity.
#' @param ports location of fishing ports
#' @param cost_per_distance cost per unit distance (deprecated - use travel_fraction instead)
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
           effort_cost_exponent = 1.2,
           travel_fraction = 0,
           ports = NULL,
           cost_per_distance = 1,
           cr_ratio = 1,
           resolution,
           patch_area = 1,
           base_effort = NULL,
           fishing_grounds = NULL,
           eta = 0.1) {
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

    # Count number of open patches (fishing_ground > 0)
    n_open_patches <- sum(fishing_grounds$fishing_ground > 0)

    # Set default effort_reference as mean per-patch effort across open patches
    effort_reference <- base_effort / n_open_patches

    # Initialize cost_per_patch based on ports or as uniform
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
        as.matrix() * sqrt(mean(patch_area))

      port_distances <- apply(matrix(port_distance[ports$patch, ], nrow = length(ports$patch)), 2, min) # calculate the distance from each patch to the port patches, then find the minimum distance

      cost_per_patch <- port_distances * cost_per_distance # calculate total travel cost per patch
    }

    # Normalize cost_per_patch to mean 1 across open patches
    # This makes it a pure spatial shape, interpretable with travel_fraction
    open_patches <- fishing_grounds$fishing_ground > 0

    # Handle edge case: if cost_per_patch is all zeros, convert to travel_fraction = 0
    if (all(cost_per_patch == 0)) {
      travel_fraction <- 0
      cost_per_patch_normalized <- rep(1, length(cost_per_patch))
    } else {
      # Normalize to mean 1 across open patches
      mean_cost_open <- mean(cost_per_patch[open_patches])
      if (mean_cost_open > 0) {
        cost_per_patch_normalized <- cost_per_patch / mean_cost_open
      } else {
        cost_per_patch_normalized <- rep(1, length(cost_per_patch))
      }
    }

    # Compute travel_weight (theta) from travel_fraction
    # theta = travel_fraction / (1 - travel_fraction)
    # This parameterization ensures travel costs are travel_fraction of total costs
    if (travel_fraction < 0 || travel_fraction >= 1) {
      stop("travel_fraction must be in [0, 1)")
    }

    if (travel_fraction == 0) {
      travel_weight <- 0
    } else {
      travel_weight <- travel_fraction / (1 - travel_fraction)
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
      travel_fraction = travel_fraction,
      travel_weight = travel_weight,
      cost_per_patch = cost_per_patch,  # Keep original for backwards reference
      cost_per_patch_normalized = cost_per_patch_normalized,
      effort_reference = effort_reference,
      fishing_grounds = fishing_grounds,
      eta = eta
    )


    return(fleet)
  } # close function

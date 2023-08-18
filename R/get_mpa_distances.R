#' Measures the distance of each cell to the nearest MPA edge and to all MPA cells
#'
#' @param mpa_locations a dataframe with at least column coordinate columns x,y,and mpa (TRUE or FALSE)
#' @param patch_area the area of each patch
#' @param resolution the resolution of the simulation system. Can supply a vector `c(x,y)` to denote an X by Y system, or one number to denote an X by X system
#'
#' @return a dataframe with two columns added.
#' `distance_to_mpa_edge` measures the distance from each patch to the nearest MPA edge in units of sqrt(patch_area), with negative values indicating areas inside an MPA
#' `total_mpa_distance` measure to total distance to all MPA cells from each patch, in units of sqrt(patch_area)
#' @export
#'
get_distance_to_mpas <-
  function(mpa_locations, resolution, patch_area = 10) {
    if (length(resolution) == 1) {
      resolution <- rep(resolution, 2)
    }
    
    # prepare MPA locations
    mpa_locations <- mpa_locations |>
      filter(mpa) |>
      dplyr::mutate(patch_name = paste(x, y, sep = "_"))
    
    # set up patch grid in case MPA locations just has MPAs in it
    patch_grid <-
      tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) |>
      dplyr::arrange(x) |>
      dplyr::mutate(
        patch_name = paste(x, y, sep = "_"),
        patch = 1:length(x),
        mpa = patch_name %in% mpa_locations$patch_name
      )
    
    out <- patch_grid
    
    # measure distance to nearest MPA border, or return Inf if there are no MPAs
    if (any(patch_grid$mpa) & any(!patch_grid$mpa)) {
      # calculate the euclidean distance between each patch
      patch_distances <- patch_grid |>
        select(x, y) |>
        dist(diag = TRUE) |>
        as.matrix()
      
      patch_distances <-
        patch_distances * sqrt(patch_area) # convert distances into the units of the system
      
      total_mpa_distance <-
        rowSums(patch_distances[, patch_grid$patch[patch_grid$mpa]]) # find the total distance to every MPA patch
      
      nearest_mpa <-
        apply(patch_distances[, patch_grid$patch[patch_grid$mpa]], 1, min) # find the distance to the nearest MPA edge
      
      nearest_fished <-
        apply(patch_distances[, patch_grid$patch[!patch_grid$mpa]], 1, min) # find the distance to the nearest fished edge
      
      distance_to_edge <-
        nearest_mpa - nearest_fished # calculate distance to nearest MPA edge, negative being inside MPA
      
      out$distance_to_mpa_edge <- distance_to_edge
      
      out$total_mpa_distance <- total_mpa_distance
    } else if (all(!patch_grid$mpa)){
      out$distance_to_mpa_edge <- Inf
      
      out$total_mpa_distance <- Inf
    } else {
      out$distance_to_mpa_edge <- 0
      
      out$total_mpa_distance <- 0
    }
    
    return(out)
    
    
  }
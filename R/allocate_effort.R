#' Allocate Effort
#'
#' This function takes an objective function in space for a given fleet
#' and returns an allocator of total effort based on that objective function
#'
#' @param obj_p a vector of length p with the objective function value in each patch
#' @param beta the responsiveness of allocation to differences in objective function, lower values make the surface smoother. >0
#' @param eps a small amount of exploratory fishing to apply everywhere available to the fleet to prevent pooling
#' @param fleet_fishable a boolean vector of length p indicating whether a given patch is open to fishing for that fleet
#'
#' @returns a numeric vector of length p that sums to 1 that can be used to distribute the total pool of effort
#' @export
#'
#' @examples
allocate_effort <- function(obj_p, beta = 0.1, eps = 0.01, fleet_fishable) {
  # robust scale, convert to median absolute deviations from the median
  z <- obj_p - median(obj_p)
  z <- z / (mad(obj_p) + 1e-6) # robust scale, constant prevents 0
  z <- pmax(pmin(z, 6), -6) # prevent extreme spikes
  w <- exp(beta * (z - max(z))) # softmax transformation


  if (sum(w) == 0) {
    w <- fleet_fishable # if numerical issues cause all w to be 0
  } else {
    w <- w * fleet_fishable
  }

  alloc <- w / sum(w) # normalize

  alloc <- (1 - eps) * alloc + eps * (fleet_fishable / sum(fleet_fishable)) # add a small bit of exploratory fishing and normalize again

  return(alloc)
}

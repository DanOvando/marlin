#' Tune Diffusion
#'
#' Tunes the diffusion parameters to achieve a desired "home range". The home range is set by
#' finding the diffusion rates that results in 95% or less of animals being within "home range" distance
#' of their starting patch after one year.
#'
#' @param home_range the desired home range, defined as the linear distance such that 95% of animals travel that distance or less from their starting location
#' @param mode one of "opt" to optimize diffusion rate or "plot" to return diffusion matrix
#'
#' @return The tuned diffusion rate
#' @export
#'
#' @examples
#'
#' home_range = 42
#' diffusion_rate <- tune_diffusion(home_range)
#' diffusion_rate
tune_diffusion <- function(home_range, mode = "opt") {

  resolution <- c(2,1)

  patch_area <- home_range^2

  patch_width = sqrt(patch_area)


  patches <- tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) |>
    dplyr::mutate(patch = 1:dplyr::n())

  centroid <- 1 # which(patches$x == round(mean(patches$x)) & patches$y == round(mean(patches$y)))

  distance <- patches |>
    dplyr::select(-patch) |>
    dist() |>
    as.matrix()

  distance_from_centroid <- distance[centroid, ] * sqrt(patch_area)

  max_distance <- max(distance_from_centroid)

  patches$distance_from_centroid <- distance_from_centroid

  foo <- function(log_diffusion_rate, home_range, patch_area, delta_t = 1, mode = "opt") {

    diffusion_rate <- exp(log_diffusion_rate)

    delta_d <- sqrt(patch_area) # length of a cell side in km

    # Mark adjacent cells
    adjacent <- distance

    adjacent[adjacent != 1] <- 0

    diffusion_matrix <- adjacent * diffusion_rate * (delta_t / delta_d^2) # adjacent times diffusion rate times simulation units

    # fill in diagonal
    diag(diffusion_matrix) <- -1 * colSums(diffusion_matrix)

    # ensure matrix class
    diffusion_matrix <- as.matrix(diffusion_matrix)

    # assign 100 individuals to the centroid patch
    next_patches <- patches

    next_patches$n_start <- 0

    next_patches$n_start[centroid] <- 1

    next_patches$n_next <- as.numeric(next_patches$n_start %*% as.matrix(expm::expm(diffusion_matrix, method = "R_Eigen"))) # calculate the realized diffusion over one time step of the model

    delta <- (next_patches$n_next[length(next_patches$n_next)] - 0.05)^2

    next_patches <- next_patches |>
      dplyr::arrange(distance_from_centroid) |>
      dplyr::mutate(
        cdist = cumsum(n_next) / sum(n_next)
      )
    #

    if (mode == "opt"){
      out <- delta

    } else {
      out <- next_patches |>
        dplyr::mutate(home_range = home_range,
                      diffusion_rate = diffusion_rate)
    }

    return(out)
  }

  best_guess <- optimise(foo, c(log(1e-12), log(1e9)), home_range = home_range, patch_area = patch_area)


  if (mode == "opt"){

    out <- exp(best_guess$minimum)

  } else {
    out <- foo(best_guess$minimum, home_range = home_range, patch_area = patch_area, mode = "plot")

  }


  return(out)
}

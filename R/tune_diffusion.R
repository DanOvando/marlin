#' Tune Diffusion
#'
#' Tunes the diffusion parameters to achieve a desired "home range". The home range is set by
#' finding the diffusion rates that results in 95% or less of animals being within "home range" distance
#' of their starting patch after one year.
#'
#' @param home_range the desired home range, defined as the linear distance such that 95% of animals travel that distance or less from their starting location
#' @param resolution the resolution of the system
#' @param patch_area the area of each patch in the system
#'
#' @return The tuned diffusion rate
#' @export
#'
#' @examples
#'
#' home_range = 42
#' diffusion_rate <- tune_diffusion(home_range)
#' diffusion_rate
tune_diffusion <- function(home_range, resolution = c(500, 1), patch_area = 4) {
  patches <- tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) |>
    dplyr::mutate(patch = 1:dplyr::n())


  centroid <- 1 # which(patches$x == round(mean(patches$x)) & patches$y == round(mean(patches$y)))

  distance <- patches |>
    dplyr::select(-patch) |>
    dist() |>
    as.matrix()

  distance_from_centroid <- distance[centroid, ] * sqrt(patch_area)

  max_distance <- max(distance_from_centroid)

  if (home_range > max_distance) {
    warning(glue::glue("supplied home range of {round(home_range,2)} is greater than maximum linear distance from centroid of patch system ({round(max_distance,2)});
          consider increasing the resolution or patch area, or reducing the home range. `home_range` has been reduced to 95% of max distance"))
    home_range <- max_distance * .95
  }

  patches$distance_from_centroid <- distance_from_centroid
  #
  # patches |>
  #   ggplot(aes(x,y,fill = distance_from_centroid)) +
  #   geom_tile()

  foo <- function(log_diffusion_rate, home_range, delta_t = 1) {
    # diffusion_rate <- 10

    diffusion_rate <- exp(log_diffusion_rate)
    # diffusion_rate <- exp(best_guess$minimum)
    # diffusion_rate <- 10000
    # diffusion_rate <- best_guess


    delta_d <- sqrt(patch_area) # length of a cell side in km

    # patch_area <- $patch_area # in units of km^2/patch

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

    next_patches$n_start[centroid] <- 100

    next_patches$n_next <- as.numeric(next_patches$n_start %*% as.matrix(expm::expm(diffusion_matrix, method = "R_Eigen"))) # calculate the realized diffusion over one time step of the model

    # next_patches$n_nex2 <- as.numeric(next_patches$n_start %*% as.matrix(expm::expm(diffusion_matrix))) # calculate the realized diffusion over one time step of the model

    #
    # next_patches |>
    #   ggplot(aes(n_next,n_nex2)) +
    #   geom_point()

    next_patches <- next_patches |>
      dplyr::arrange(distance_from_centroid) |>
      dplyr::mutate(
        cdist = cumsum(n_next) / sum(n_next),
        delta = (cdist - 0.95)^2
      )
    #

    # next_patches |>
    #   ggplot(aes(distance_from_centroid, cdist)) +
    #   geom_line() +
    #   geom_hline(yintercept = 0.95) +
    #   geom_vline(xintercept = home_range) +
    #   scale_y_continuous(limits = c(0, 1.1))


    diffusion_distance <- next_patches$distance_from_centroid[which.min(next_patches$delta)]

    delta <- as.numeric((diffusion_distance - home_range)^2)

    return(delta)
  }

  # home_range <- max_distance *.01

  # delta <- 999
  #
  # tol <- 0.1
  #
  # range_vals <-log(c(.01,50000))
  #
  # last_min <- 1000
  #
  # while (delta >= tol){
  #
  #   possible_values = seq(range_vals[1], range_vals[2], length.out = 25)
  #   ss = purrr::map_dbl(possible_values,foo, home_range = home_range, .progress = TRUE)
  #   plot(possible_values, ss)
  #
  #   best_guess <- possible_values[which.min(ss)]
  #
  #   range_vals <- c(best_guess*.5, best_guess * 10)
  #
  #   delta <- (min(ss) - last_min)^2
  #
  #   last_min <- min(ss)
  #
  #
  # }

  best_guess <- optimise(foo, c(log(1e-6), log(100000)), home_range = home_range)

  diffusion_rate <- exp(best_guess$minimum)
}

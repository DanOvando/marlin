#' Simulate species habitats with correlations across space and among species
#'
#' @param critters a vector of critters names of length n_species representing the number of species to be simulated
#' @param kp a rate parameter that governs the rate at which spatial cells become decorrelated with distance. Lower values mean a smoother (more correlated) habitat map
#' @param critter_correlations leave as NA to randomly simulate correlations across species. Otherwise, a n_species x n_species correlation matrix
#' @param resolution the resolution of the system
#' @param patch_area
#'
#' @return a list with simulated critter habitats
#' @export
#'
#' @examples
#'
#' n_species <- 6
#'
#' critters <- paste0(fruit[1:n_species], "_fish")
#'
#' resolution <- c(4, 20)
#'
#' patch_area <- 4
#'
#' kp <- .1
#'
#'
#' species_distributions <- sim_habitat(critters = critters, resolution = resolution, patch_area = patch_area, kp = kp)
sim_habitat <-
  function(critters,
           kp,
           critter_correlations = NA,
           resolution,
           patch_area,
           rescale_habitat = TRUE,
           max_delta = 3,
           max_abs_cor = 1,
           output = "df") {
    if (length(resolution) == 1) {
      resolution <- rep(resolution, 2)
    }

    n_species <- length(critters)

    patches <- prod(resolution) # the total number of patches

    patch_width <- sqrt(patch_area) # the width of a patch

    grid <-
      tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) # set up spatial grid

    # calculate the distance
    distances <-
      grid |>
      dist() |>
      as.matrix() * patch_width

    # following Thorson & Barnett 2017 doi:10.1093/icesjms/fsw193
    # kp = 1e-1 # the rate parameter; closer to 0, the farther correlation among patches persists

    n <- 1 # set N to 1 per paper

    H <- 1 # matrix(rep(1,4), ncol = 2) # placeholder for now since I can't figure out how to implement anisotropy uet

    sigma <- 1 # placeholder which I believe should be left at 1

    # set up the spatial correlation matrix
    spatial_correlations <- matrix(0, nrow = (patches), ncol = (patches))

    for (x in 1:ncol(distances)) {
      for (y in 1:nrow(distances)) {
        if (distances[y, x] != 0) {
          spatial_correlations[y, x] <-
            sigma^2 / (2^(n - 1) * gamma(n)) * (kp * abs(distances[y, x] * H))^
              n * besselK(kp * abs(distances[y, x] * H), 1)
        } else {
          spatial_correlations[y, x] <- 1
        }
      }
    }

    # create species correlation matrix

    if (all(is.na(critter_correlations))) {
      n_species_cores <- n_species * (n_species + 1) / 2 - n_species
      # n_species * (n_species + 1) / 2 # formula for the number of elements in the upper triangle of an n x n matric

      core_matrix <- matrix(0, nrow = n_species, ncol = n_species)

      species_cores <-
        runif(n_species_cores,
          min = -max_abs_cor,
          max = max_abs_cor
        ) # randomly generate correlations among species
      # Fill in the upper triangle of the matrix
      core_matrix[upper.tri(core_matrix)] <- species_cores

      lower_triangle <- t(core_matrix)

      critter_correlations <- core_matrix + lower_triangle

      diag(critter_correlations) <- 1
    }

    species_x_space <-
      as.matrix(Matrix::nearPD(
        Matrix::kronecker(spatial_correlations, critter_correlations)
      )$mat)

    # generate a random species distribution for each species in space
    habitats <-
      MASS::mvrnorm(n = 1, rep(0, ncol(species_x_space)), Sigma = species_x_space)

    species_distributions <-
      tidyr::expand_grid(critter = critters, patch = 1:patches) |>
      cbind(grid) |>
      dplyr::arrange(patch, rev(critter)) |>
      dplyr::mutate(habitat = habitats)

    if (rescale_habitat) {
      species_distributions$habitat <- scales::rescale(species_distributions$habitat, to = c(0, log(max_delta)))
    }

    check_species_cores <- species_distributions |>
      dplyr::select(patch, critter, habitat) |>
      tidyr::pivot_wider(names_from = "critter", values_from = "habitat")

    final_species_cores <- cor(check_species_cores[, -1])


    if (output == "list") {
      critter_habitats <- species_distributions |>
        dplyr::group_by(critter) |>
        tidyr::nest() |>
        dplyr::mutate(
          habitat = purrr::map(
            data,
            \(x) x |>
              dplyr::select(-patch) |>
              tidyr::pivot_wider(names_from = x, values_from = habitat) |>
              dplyr::select(-y) %>%
              as.matrix()
          )
        )

      species_distributions <- critter_habitats$habitat |>
        purrr::set_names(critter_habitats$critter)
    }


    out <- list(
      critter_distributions = species_distributions,
      critter_correlations = final_species_cores,
      wtf = critter_correlations
    )
  }

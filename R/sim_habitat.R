#' Simulate Spatially Correlated Species Habitat
#'
#' @description
#' Generates a set of spatially correlated habitat-quality surfaces for one or
#' more species using a Matern covariance function (following Thorson &
#' Barnett, 2017). Habitat values can be positively or negatively correlated
#' across species, reflecting ecological affinities or competitive exclusion.
#' The resulting matrices are the primary input for the \code{habitat} argument
#' of \code{\link{create_critter}}.
#'
#' @details
#' Spatial correlation among patches follows a Matern covariance with
#' smoothness \eqn{n = 1}. The rate parameter \code{kp} controls the spatial
#' decorrelation length: a value of \code{0.1} produces broad, smooth habitat
#' gradients, while \code{1.0} produces patchy, rapidly decorrelating maps.
#'
#' Cross-species correlations are drawn uniformly from
#' \eqn{[-max\_abs\_cor, -min\_abs\_cor] \cup [min\_abs\_cor, max\_abs\_cor]} when
#' \code{critter_correlations = NA}. Supply a full \eqn{n \times n} correlation
#' matrix to fix species correlations exactly.
#'
#' When \code{output = "list"}, habitat matrices are returned in \code{[ny, nx]}
#' format with column names \code{1:nx} and row names \code{1:ny}, matching
#' the convention expected by \code{\link{create_critter}}.
#'
#' @param critters Character vector of species names (length \eqn{n}).
#'   Names are carried through to the output list.
#' @param kp Positive numeric. Matern spatial range parameter. Lower values
#'   produce smoother, more correlated spatial distributions; higher values
#'   produce patchier maps. Typical values: \code{0.05}--\code{0.5}.
#' @param critter_correlations \eqn{n \times n} numeric correlation matrix
#'   specifying inter-species correlations in habitat quality, or \code{NA}
#'   (default) to draw random correlations. Diagonal must be 1; off-diagonal
#'   entries must be in \code{(-1, 1)}.
#' @param resolution Integer scalar or length-2 integer vector \code{c(nx, ny)}
#'   giving grid dimensions. Matches the \code{resolution} argument of
#'   \code{\link{create_critter}}.
#' @param patch_area Numeric. Area of each patch (km^2). Used to compute
#'   inter-patch distances in km.
#' @param rescale_habitat Logical. If \code{TRUE} (default), rescales all
#'   habitat values to \code{[0, log(max_delta)]}. The log scale means
#'   that values drive movement taxis multiplicatively.
#' @param max_delta Positive numeric. Upper bound of habitat after rescaling;
#'   passed to \code{\link[scales]{rescale}}. Default \code{3}.
#' @param max_abs_cor Numeric in (0, 1]. Maximum absolute value of randomly
#'   generated inter-species correlations. Must exceed \code{min_abs_cor}.
#' @param min_abs_cor Numeric in [0, 1). Minimum absolute value of randomly
#'   generated inter-species correlations. Default \code{0}.
#' @param output Character. One of \code{"df"} (default) or \code{"list"}.
#'   \describe{
#'     \item{\code{"df"}}{Returns a data frame with columns
#'       \code{critter}, \code{patch}, \code{x}, \code{y}, and
#'       \code{habitat}.}
#'     \item{\code{"list"}}{Returns a named list (one entry per species) of
#'       \code{[ny, nx]} habitat matrices suitable for direct use as the
#'       \code{habitat} argument to \code{\link{create_critter}}.}
#'   }
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{critter_distributions}}{Either a long data frame (when
#'     \code{output = "df"}) or a named list of \code{[ny, nx]} habitat
#'     matrices (when \code{output = "list"}).}
#'   \item{\code{critter_correlations}}{An \eqn{n \times n} correlation
#'     matrix of the realised inter-species habitat correlations.}
#'   \item{\code{wtf}}{The input or generated correlation matrix (before
#'     nearest positive-definite projection).}
#' }
#'
#' @references
#' Thorson, J.T. & Barnett, L.A.K. (2017). Comparing estimates of abundance
#' trends and distribution shifts using single- and multispecies models of
#' fishes and biogenic habitat. *ICES Journal of Marine Science*, 74(5),
#' 1311--1321. \doi{10.1093/icesjms/fsw193}
#'
#' @seealso \code{\link{create_critter}}, \code{\link{simmar}}
#'
#' @export
#'
#' @examples
#' # Simulate habitat for 3 species on a 10x10 grid
#' hab <- sim_habitat(
#'   critters   = c("tuna", "grouper", "snapper"),
#'   kp         = 0.1,
#'   resolution = c(10, 10),
#'   patch_area = 4,
#'   output     = "list"
#' )
#'
#' # hab$critter_distributions$tuna  # [10, 10] matrix
sim_habitat <-
  function(critters,
           kp,
           critter_correlations = NA,
           resolution,
           patch_area,
           rescale_habitat = TRUE,
           max_delta = 3,
           max_abs_cor = 1,
           min_abs_cor = 0,
           output = "df") {

    if (max_abs_cor <= min_abs_cor){
      stop("max_abs_cor must be greater than min_abs_cor")
    }

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
        runif(n_species_cores, min = min_abs_cor, max = max_abs_cor) * sample(c(-1,1), n_species_cores, replace =  TRUE)


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

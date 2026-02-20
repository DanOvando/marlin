#' Prepare Instantaneous Movement (Generator) Matrix
#'
#' @description
#' Converts a sparse adjacency-weighted matrix (output of the habitat-taxis
#' step in \code{\link{simmar}}) into an instantaneous transition rate matrix
#' (generator matrix) suitable for matrix exponentiation via
#' \code{expm::expm}. The diagonal is set to minus the column sums so that
#' each column sums to zero, as required for a continuous-time Markov chain.
#'
#' @details
#' Movement in marlin is modelled as a diffusion-taxis process following
#' Thorson & Barnett (2017). The habitat taxis component modifies the
#' off-diagonal elements of the base diffusion matrix to reflect drift towards
#' high-quality habitat. \code{prep_movement} completes the generator matrix
#' by filling in the diagonal, after which \code{expm::expm} is called to
#' obtain the discrete-time transition matrix over one time step.
#'
#' @param multiplier Sparse matrix (from \code{Matrix::sparseMatrix}).
#'   Off-diagonal elements represent the instantaneous transition rates from
#'   patch j (column) to patch i (row), derived by element-wise multiplication
#'   of the base diffusion matrix and the habitat-taxis multiplier matrix.
#' @param resolution Integer vector of length 2. Grid dimensions \code{c(nx,
#'   ny)}. Currently unused inside the function body but retained for
#'   compatibility with earlier code.
#'
#' @return A sparse matrix with the same dimensions as \code{multiplier},
#'   with diagonal elements set to \eqn{-\sum_{i \neq j} m_{ij}} (negative
#'   column sums of the off-diagonal).
#'
#' @references
#' Thorson, J.T. & Barnett, L.A.K. (2017). Comparing estimates of abundance
#' trends and distribution shifts using single- and multispecies models of
#' fishes and biogenic habitat. *ICES Journal of Marine Science*, 74(5),
#' 1311--1321. \doi{10.1093/icesjms/fsw193}
#'
#' @seealso \code{\link{simmar}}, \code{\link{sparsify_movement}}
#'
#' @export
#' Prepare movement matrix
#'
#' @param multiplier multiplier for adjacency matrix
#' @param time_step time step in question
#' @param resolution spatial resolution
#'
#' @return a prepared movement matrix
#' @export
#'
prep_movement <-
  function(multiplier,
           resolution) {
    # reminder,
    #
    # .2^2 == exp(2 * log(.2)), hence notation in Thorson et al. going back and forth between different structures


    # # set up spatial grid
    # adjacent <-
    #   tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) %>%
    #   dist() %>%
    #   as.matrix()
    #
    # # Mark adjacent cells
    # adjacent[adjacent != 1] <- 0
    #
    # # mark off water
    #
    # water <- !is.na(multiplier)
    #
    # multiplier[is.na(multiplier)] <- 0
    #
    # adjacent <- adjacent * water # set as non-adjacent any patches that contain land
    #
    # # generate base movement matrix
    # move_mat <- adjacent * multiplier

    # fill in diagonal
    move_mat <- multiplier
    diag(move_mat) <- -1 * Matrix::colSums(move_mat)

    # ensure matrix class
    # move_mat <- as.matrix(move_mat)

    return(move_mat)
  } # close calc_move_mat

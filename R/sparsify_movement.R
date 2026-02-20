#' Sparsify a Dense Transition Matrix
#'
#' @description
#' Converts a dense, column-stochastic transition matrix (output of
#' \code{expm::expm}) into a sparse matrix by retaining only the largest
#' destination probabilities for each origin patch, until the cumulative
#' retained mass reaches the \code{retain} threshold. The kept values are
#' renormalised so each column still sums to exactly 1.
#'
#' @details
#' This operation dramatically reduces memory use and speeds up subsequent
#' matrix-vector multiplications for large grids. By default (\code{retain =
#' 0.999}), at most 0.1\% of probability mass is dropped per patch, so
#' movement dynamics are essentially unchanged.
#'
#' Called internally by \code{\link{simmar}} after computing the
#' discrete-time transition matrix via matrix exponentiation.
#'
#' @param T_dense Numeric matrix (patches x patches). Dense column-stochastic
#'   transition matrix; each column gives the distribution of destinations
#'   for individuals starting in that patch.
#' @param retain Numeric in (0, 1]. Minimum cumulative probability mass to
#'   retain per origin patch. Higher values preserve more of the tails of
#'   the movement distribution at the cost of denser storage. Default
#'   \code{0.999}.
#'
#' @return A sparse matrix (class \code{"dgCMatrix"}) with the same
#'   dimensions as \code{T_dense}, with near-zero off-diagonal entries
#'   dropped and columns renormalised to sum to 1.
#'
#' @seealso \code{\link{prep_movement}}, \code{\link{simmar}}
#'
#' @keywords internal
sparsify_transition <- function(trans_mat, retain = 0.999) {

  P <- nrow(trans_mat)

  ii_list <- vector("list", P)
  jj_list <- vector("list", P)
  xx_list <- vector("list", P)

  for (from in seq_len(P)) {

    col_vals <- trans_mat[, from]
    nz <- which(col_vals > 0)

    if (length(nz) == 0L) next

    vals <- col_vals[nz]

    # Match old behavior exactly (including tie-breaking)
    ord <- order(vals, decreasing = TRUE)
    idx_sorted  <- nz[ord]
    vals_sorted <- vals[ord]

    cum_mass <- cumsum(vals_sorted)
    k <- which(cum_mass >= retain)[1]
    if (is.na(k)) k <- length(vals_sorted)

    keep_rows <- idx_sorted[seq_len(k)]
    keep_vals <- col_vals[keep_rows]
    keep_vals <- keep_vals / sum(keep_vals)

    ii_list[[from]] <- keep_rows
    jj_list[[from]] <- rep.int(from, length(keep_rows))
    xx_list[[from]] <- keep_vals
  }

  ii <- unlist(ii_list, use.names = FALSE)
  jj <- unlist(jj_list, use.names = FALSE)
  xx <- unlist(xx_list, use.names = FALSE)

  Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(P, P))
}

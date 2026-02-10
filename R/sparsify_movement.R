# -----------------------------------------------------------------------------
# sparsify_transition()
#
# Purpose:
# Convert a dense transition matrix T (patch × patch, column-stochastic)
# into a sparse matrix while preserving ≥ retain fraction of outgoing
# probability mass from each origin patch.
#
# Conceptually:
#   Each column = distribution of where individuals starting in patch "from"
#   move after one time step.
#
# Algorithm:
#   For each origin patch (column):
#     1. Sort destination probabilities from largest → smallest
#     2. Keep the minimum number of destinations needed to reach cumulative
#        probability ≥ retain (e.g. 0.999)
#     3. Renormalize kept values so column sums exactly to 1
#     4. Store only those entries → sparse matrix
#
# Guarantees:
#   - No probability mass leakage (columns always sum to 1)
#   - Dominant dispersal kernel structure preserved
#   - Long-distance tiny tails removed (→ sparsity)
#
# Assumptions:
#   - T is non-negative
#   - Columns of T approximately sum to 1 (post expm movement matrix)
# -----------------------------------------------------------------------------

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

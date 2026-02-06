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

sparsify_transition <- function(T, retain = 0.999) {

  # Number of spatial patches (matrix is P × P)
  P <- nrow(T)

  # Storage for sparse matrix triplet form
  #   ii = destination patch index (row)
  #   jj = origin patch index (column)
  #   xx = transition probability value
  ii <- integer(0)
  jj <- integer(0)
  xx <- numeric(0)

  # Loop over origin patches (columns = "from")
  for (from in seq_len(P)) {

    # Probability distribution of destinations from this origin patch
    col_vals <- T[, from]

    # Ignore exact zeros (no movement probability → no need to store)
    nz <- which(col_vals > 0)

    # If no movement from this patch, skip
    if (length(nz) == 0) next

    # Extract nonzero values
    vals <- col_vals[nz]

    # Sort by probability mass (largest first)
    ord <- order(vals, decreasing = TRUE)

    idx_sorted  <- nz[ord]     # destination patch indices
    vals_sorted <- vals[ord]   # corresponding probabilities

    # Cumulative probability mass retained as we add destinations
    cum_mass <- cumsum(vals_sorted)

    # Minimum number of destinations needed to retain target mass
    k <- which(cum_mass >= retain)[1]

    # If retain never reached (rare; numerical or degenerate cases),
    # keep everything
    if (is.na(k)) k <- length(vals_sorted)

    # Indices and values to keep
    keep_rows <- idx_sorted[seq_len(k)]
    keep_vals <- col_vals[keep_rows]

    # Renormalize so column sums exactly to 1 (mass conservation)
    keep_vals <- keep_vals / sum(keep_vals)

    # Append to sparse storage vectors
    ii <- c(ii, keep_rows)
    jj <- c(jj, rep.int(from, length(keep_rows)))
    xx <- c(xx, keep_vals)
  }

  # Build sparse transition matrix
  #   rows = destination patch
  #   cols = origin patch
  Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(P, P))
}

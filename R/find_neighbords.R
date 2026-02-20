#' Find 4-Neighbour (Rook) Adjacency Matrix for a Spatial Grid
#'
#' @description
#' Builds a sparse adjacency matrix \eqn{A} (patches x patches) where
#' \eqn{A[p, q] = 1} if patches \eqn{p} and \eqn{q} share an edge (rook
#' adjacency), and 0 otherwise. An optional water mask restricts edges so that
#' land patches are excluded from the adjacency structure.
#'
#' @details
#' Patch indexing follows the marlin convention:
#' \deqn{p = (x - 1) \times Y + y \quad \text{(y varies fastest)}}
#' which matches the ordering produced by
#' \code{tidyr::expand_grid(x = 1:X, y = 1:Y) |> mutate(patch = 1:n())}.
#'
#' Used internally during movement-matrix construction in
#' \code{\link{simmar}} to identify which patches are adjacent when
#' updating habitat taxis.
#'
#' @param resolution Integer vector of length 2: \code{c(X, Y)}, where
#'   \code{X} is the number of columns (x-direction) and \code{Y} the number
#'   of rows (y-direction).
#' @param water Logical vector of length \code{X * Y} in patch order
#'   (\code{TRUE} = water/open, \code{FALSE} = land/closed). Pass
#'   \code{rep(TRUE, X * Y)} for a fully open grid.
#'
#' @return A symmetric sparse matrix (class \code{"dgCMatrix"}) of dimensions
#'   \eqn{P \times P} (where \eqn{P = X \times Y}) with dimnames
#'   \code{"1"}, \code{"2"}, ..., \code{"P"}.
#'
#' @seealso \code{\link{prep_movement}}, \code{\link{simmar}}
#'
#' @export
find_neighbors <- function(resolution, water_mask) {

  # --- Basic input checks ------------------------------------------------------
  stopifnot(length(resolution) == 2)

  n_cols_x <- as.integer(resolution[1])  # X in your notation (number of columns)
  n_rows_y <- as.integer(resolution[2])  # Y in your notation (number of rows)
  stopifnot(n_cols_x >= 1L, n_rows_y >= 1L)

  n_patches <- n_cols_x * n_rows_y       # P = total number of patches in the grid

  stopifnot(length(water_mask) == n_patches)
  water_mask <- as.logical(water_mask)

  # Patch ids are 1..P, in the same order as expand_grid(x=1:X, y=1:Y)
  patch_id <- seq_len(n_patches)

  # --- Convert patch id -> (x, y) coordinates --------------------------------
  # Inverse mapping for:
  #   p = (x - 1)*Y + y   where y varies fastest
  #
  # Given p in 1..P:
  #   x = floor((p-1)/Y) + 1
  #   y = ((p-1) mod Y) + 1
  x_coord <- ((patch_id - 1L) %/% n_rows_y) + 1L
  y_coord <- ((patch_id - 1L) %%  n_rows_y) + 1L

  # --- Build the (directed) edge list for 4-neighbors -------------------------
  # We only generate "right" and "up" edges, then mirror them to make symmetry.
  #
  # Right neighbor means (x+1, y):
  #   Because y varies fastest, stepping one column to the right adds +Y to the patch id.
  #
  # Up neighbor means (x, y+1):
  #   Stepping one row up adds +1 to the patch id (within the same column).
  from_patch <- integer(0)  # i indices
  to_patch   <- integer(0)  # j indices

  # (1) Right edges: for patches not in the last column (x < X)
  idx_has_right <- which(x_coord < n_cols_x)
  from_patch <- c(from_patch, patch_id[idx_has_right])
  to_patch   <- c(to_patch,   patch_id[idx_has_right] + n_rows_y)

  # (2) Up edges: for patches not in the top row (y < Y)
  idx_has_up <- which(y_coord < n_rows_y)
  from_patch <- c(from_patch, patch_id[idx_has_up])
  to_patch   <- c(to_patch,   patch_id[idx_has_up] + 1L)

  # --- Apply water mask --------------------------------------------------------
  # We only keep an adjacency edge if both endpoints are water.
  # (Land patches remain in the matrix as isolated nodes, which keeps dimensions P x P.)
  keep_edge <- water_mask[from_patch] & water_mask[to_patch]
  from_patch <- from_patch[keep_edge]
  to_patch   <- to_patch[keep_edge]

  # --- Make the adjacency undirected (symmetric) ------------------------------
  # The dense dist() baseline yields a symmetric adjacency matrix.
  # We create both directions for every kept edge.
  i_index <- c(from_patch, to_patch)
  j_index <- c(to_patch,   from_patch)

  # Give dimnames so identical() can match dist()'s as.matrix output.
  patch_names <- as.character(seq_len(n_patches))

  # --- Construct sparse adjacency matrix --------------------------------------
  # Entries are 1 for adjacent pairs, 0 otherwise.
  Matrix::sparseMatrix(
    i = i_index,
    j = j_index,
    x = 1,
    dims = c(n_patches, n_patches),
    dimnames = list(patch_names, patch_names),
    giveCsparse = TRUE
  )
}

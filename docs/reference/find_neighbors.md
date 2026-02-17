# Find 4-Neighbour (Rook) Adjacency Matrix for a Spatial Grid

Builds a sparse adjacency matrix \\A\\ (patches x patches) where \\A\[p,
q\] = 1\\ if patches \\p\\ and \\q\\ share an edge (rook adjacency), and
0 otherwise. An optional water mask restricts edges so that land patches
are excluded from the adjacency structure.

## Usage

``` r
find_neighbors(resolution, water_mask)
```

## Arguments

- resolution:

  Integer vector of length 2: `c(X, Y)`, where `X` is the number of
  columns (x-direction) and `Y` the number of rows (y-direction).

- water:

  Logical vector of length `X * Y` in patch order (`TRUE` = water/open,
  `FALSE` = land/closed). Pass `rep(TRUE, X * Y)` for a fully open grid.

## Value

A symmetric sparse matrix (class `"dgCMatrix"`) of dimensions \\P \times
P\\ (where \\P = X \times Y\\) with dimnames `"1"`, `"2"`, ..., `"P"`.

## Details

Patch indexing follows the marlin convention: \$\$p = (x - 1) \times Y +
y \quad \text{(y varies fastest)}\$\$ which matches the ordering
produced by
`tidyr::expand_grid(x = 1:X, y = 1:Y) |> mutate(patch = 1:n())`.

Used internally during movement-matrix construction in
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md) to
identify which patches are adjacent when updating habitat taxis.

## See also

[`prep_movement`](https://danovando.github.io/marlin/reference/prep_movement.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)

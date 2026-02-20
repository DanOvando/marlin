# Sparsify a Dense Transition Matrix

Converts a dense, column-stochastic transition matrix (output of
[`expm::expm`](https://rdrr.io/pkg/expm/man/expm.html)) into a sparse
matrix by retaining only the largest destination probabilities for each
origin patch, until the cumulative retained mass reaches the `retain`
threshold. The kept values are renormalised so each column still sums to
exactly 1.

## Usage

``` r
sparsify_transition(trans_mat, retain = 0.999)
```

## Arguments

- retain:

  Numeric in (0, 1\]. Minimum cumulative probability mass to retain per
  origin patch. Higher values preserve more of the tails of the movement
  distribution at the cost of denser storage. Default `0.999`.

- T_dense:

  Numeric matrix (patches x patches). Dense column-stochastic transition
  matrix; each column gives the distribution of destinations for
  individuals starting in that patch.

## Value

A sparse matrix (class `"dgCMatrix"`) with the same dimensions as
`T_dense`, with near-zero off-diagonal entries dropped and columns
renormalised to sum to 1.

## Details

This operation dramatically reduces memory use and speeds up subsequent
matrix-vector multiplications for large grids. By default
(`retain = 0.999`), at most 0.1\\ movement dynamics are essentially
unchanged.

Called internally by
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md) after
computing the discrete-time transition matrix via matrix exponentiation.

## See also

[`prep_movement`](https://danovando.github.io/marlin/reference/prep_movement.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)

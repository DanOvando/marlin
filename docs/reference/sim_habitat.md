# Simulate Spatially Correlated Species Habitat

Generates a set of spatially correlated habitat-quality surfaces for one
or more species using a Matern covariance function (following Thorson &
Barnett, 2017). Habitat values can be positively or negatively
correlated across species, reflecting ecological affinities or
competitive exclusion. The resulting matrices are the primary input for
the `habitat` argument of
[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md).

## Usage

``` r
sim_habitat(
  critters,
  kp,
  critter_correlations = NA,
  resolution,
  patch_area,
  rescale_habitat = TRUE,
  max_delta = 3,
  max_abs_cor = 1,
  min_abs_cor = 0,
  output = "df"
)
```

## Arguments

- critters:

  Character vector of species names (length \\n\\). Names are carried
  through to the output list.

- kp:

  Positive numeric. Matern spatial range parameter. Lower values produce
  smoother, more correlated spatial distributions; higher values produce
  patchier maps. Typical values: `0.05`–`0.5`.

- critter_correlations:

  \\n \times n\\ numeric correlation matrix specifying inter-species
  correlations in habitat quality, or `NA` (default) to draw random
  correlations. Diagonal must be 1; off-diagonal entries must be in
  `(-1, 1)`.

- resolution:

  Integer scalar or length-2 integer vector `c(nx, ny)` giving grid
  dimensions. Matches the `resolution` argument of
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md).

- patch_area:

  Numeric. Area of each patch (km^2). Used to compute inter-patch
  distances in km.

- rescale_habitat:

  Logical. If `TRUE` (default), rescales all habitat values to
  `[0, log(max_delta)]`. The log scale means that values drive movement
  taxis multiplicatively.

- max_delta:

  Positive numeric. Upper bound of habitat after rescaling; passed to
  [`rescale`](https://scales.r-lib.org/reference/rescale.html). Default
  `3`.

- max_abs_cor:

  Numeric in (0, 1\]. Maximum absolute value of randomly generated
  inter-species correlations. Must exceed `min_abs_cor`.

- min_abs_cor:

  Numeric in \[0, 1). Minimum absolute value of randomly generated
  inter-species correlations. Default `0`.

- output:

  Character. One of `"df"` (default) or `"list"`.

  `"df"`

  :   Returns a data frame with columns `critter`, `patch`, `x`, `y`,
      and `habitat`.

  `"list"`

  :   Returns a named list (one entry per species) of `[ny, nx]` habitat
      matrices suitable for direct use as the `habitat` argument to
      [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md).

## Value

A list with three elements:

- `critter_distributions`:

  Either a long data frame (when `output = "df"`) or a named list of
  `[ny, nx]` habitat matrices (when `output = "list"`).

- `critter_correlations`:

  An \\n \times n\\ correlation matrix of the realised inter-species
  habitat correlations.

- `wtf`:

  The input or generated correlation matrix (before nearest
  positive-definite projection).

## Details

Spatial correlation among patches follows a Matern covariance with
smoothness \\n = 1\\. The rate parameter `kp` controls the spatial
decorrelation length: a value of `0.1` produces broad, smooth habitat
gradients, while `1.0` produces patchy, rapidly decorrelating maps.

Cross-species correlations are drawn uniformly from
`[-max_abs_cor, -min_abs_cor] [min_abs_cor, max_abs_cor]` when
`critter_correlations = NA`. Supply a full \\n \times n\\ correlation
matrix to fix species correlations exactly.

When `output = "list"`, habitat matrices are returned in `[ny, nx]`
format with column names `1:nx` and row names `1:ny`, matching the
convention expected by
[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md).

## References

Thorson, J.T. & Barnett, L.A.K. (2017). Comparing estimates of abundance
trends and distribution shifts using single- and multispecies models of
fishes and biogenic habitat. *ICES Journal of Marine Science*, 74(5),
1311–1321.
[doi:10.1093/icesjms/fsw193](https://doi.org/10.1093/icesjms/fsw193)

## See also

[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)

## Examples

``` r
# Simulate habitat for 3 species on a 10x10 grid
hab <- sim_habitat(
  critters   = c("tuna", "grouper", "snapper"),
  kp         = 0.1,
  resolution = c(10, 10),
  patch_area = 4,
  output     = "list"
)

# hab$critter_distributions$tuna  # [10, 10] matrix
```

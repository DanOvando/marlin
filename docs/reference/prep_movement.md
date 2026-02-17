# Prepare Instantaneous Movement (Generator) Matrix

Converts a sparse adjacency-weighted matrix (output of the habitat-taxis
step in
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)) into
an instantaneous transition rate matrix (generator matrix) suitable for
matrix exponentiation via
[`expm::expm`](https://rdrr.io/pkg/expm/man/expm.html). The diagonal is
set to minus the column sums so that each column sums to zero, as
required for a continuous-time Markov chain.

## Usage

``` r
prep_movement(multiplier, resolution)
```

## Arguments

- multiplier:

  multiplier for adjacency matrix

- resolution:

  spatial resolution

- time_step:

  time step in question

## Value

A sparse matrix with the same dimensions as `multiplier`, with diagonal
elements set to \\-\sum\_{i \neq j} m\_{ij}\\ (negative column sums of
the off-diagonal).

a prepared movement matrix

## Details

Movement in marlin is modelled as a diffusion-taxis process following
Thorson & Barnett (2017). The habitat taxis component modifies the
off-diagonal elements of the base diffusion matrix to reflect drift
towards high-quality habitat. `prep_movement` completes the generator
matrix by filling in the diagonal, after which
[`expm::expm`](https://rdrr.io/pkg/expm/man/expm.html) is called to
obtain the discrete-time transition matrix over one time step.

## References

Thorson, J.T. & Barnett, L.A.K. (2017). Comparing estimates of abundance
trends and distribution shifts using single- and multispecies models of
fishes and biogenic habitat. *ICES Journal of Marine Science*, 74(5),
1311–1321.
[doi:10.1093/icesjms/fsw193](https://doi.org/10.1093/icesjms/fsw193)

## See also

[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
`sparsify_movement`

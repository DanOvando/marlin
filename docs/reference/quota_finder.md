# Objective Function for Quota-Constrained Fishing Mortality

Computes the squared difference between realised total catch and a
target quota under a uniformly scaled fishing mortality. Used as the
objective for `optim` inside
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)
whenever a species-level catch quota is active and would be exceeded.

## Usage

``` r
quota_finder(
  fmult,
  quota,
  fauna,
  current_season,
  movement,
  f_p_a,
  last_n_p_a,
  f_p_a_fl,
  critter,
  rec_devs,
  patches,
  ages,
  fleets
)
```

## Arguments

- fmult:

  Numeric in 0, 1. Effort multiplier applied uniformly to `f_p_a` before
  simulating the population.

- quota:

  Numeric. Target total catch (in numbers) that must not be exceeded.

- fauna:

  Named list of fauna objects.

- current_season:

  Integer. Current season index.

- movement:

  List of movement matrices (one per season block).

- f_p_a:

  Numeric matrix `[patches, ages]` of baseline fishing mortality before
  quota reduction.

- last_n_p_a:

  Numeric matrix `[patches, ages]` of numbers at the start of the
  current step.

- f_p_a_fl:

  3-D array `[patches, ages, fleets]` of proportional fishing mortality
  share by fleet.

- critter:

  Character. Name of the target species.

- rec_devs:

  Numeric vector of recruitment deviates (length = patches).

- patches:

  Integer. Number of patches.

- ages:

  Integer. Number of age classes.

- fleets:

  Named list of fleet objects.

## Value

Numeric scalar: squared deviation of total catch from `quota`. Minimised
by `optim` to find the binding effort multiplier.

## Details

When a quota is binding,
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md) calls
`optim` with this function to find the scalar `fmult` in 0, 1 that
reduces the fishing mortality matrix `f_p_a` (and hence catch) to
exactly the quota level. All fleets experience the same proportional
reduction.

## See also

[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)

# Assign MSY-Based Reference Points to Each Species

Finds the effort level that maximises yield for each species in `fauna`
by calling `optim` on
[`find_msy`](https://danovando.github.io/marlin/reference/find_msy.md),
then stores the resulting MSY, SSBmsy, Bmsy, and Umsy in each critter's
`ref_points` field.

## Usage

``` r
assign_ref_points(fauna, fleets)
```

## Arguments

- fauna:

  a list of critters

- fleets:

  a list of cleets

## Value

A copy of `fauna` with the `ref_points` field populated for each
species.

a fauna object but with MSY based reference points included for each
critter

## Details

The optimisation searches over a scalar effort multiplier applied to the
first fleet's `base_effort` (range 0.001–10). MSY is estimated for each
species independently (i.e. holding all else equal; no joint
multi-species MSY is computed). Results are stored in each critter's
`ref_points` data frame with columns `ssb_msy`, `b_msy`, `n_msy`, `msy`,
`u_msy`, `base_e_msy_mult`, and `base_e_msy`.

## See also

[`find_msy`](https://danovando.github.io/marlin/reference/find_msy.md),
[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md),
[`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md)

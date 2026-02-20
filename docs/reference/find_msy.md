# Find MSY Referemce Points

Deprecated and likely not to come back, only really possible with one
fleet, one critter, and no space, so in other words not what `marlin` is
intended for

Estimates maximum sustainable yield (MSY) conditions for one species by
running a simulation under a scaled effort level and returning either
the negative total yield (for minimisation by `optim`) or the full
[`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md)
output at that effort level.

Mostly an internal helper for
[`assign_ref_points`](https://danovando.github.io/marlin/reference/assign_ref_points.md),
but exported for use in custom MSY calculations.

## Usage

``` r
find_msy(effort_mult, fauna, fleets, opt = TRUE, target_critter)

find_msy(effort_mult, fauna, fleets, opt = TRUE, target_critter)
```

## Arguments

- effort_mult:

  Numeric. Scalar multiplier applied to the first fleet's `base_effort`
  in `fleets`. The optimiser searches over this value to find the effort
  that maximises yield.

- fauna:

  Named list of fauna objects from
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md).

- fleets:

  Named list of fleet objects from
  [`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md),
  already tuned with
  [`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md).

- opt:

  Logical. If `TRUE` (default), returns negative total yield (scalar)
  for use as an `optim` objective. If `FALSE`, returns the full
  [`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md)
  output at the given effort level.

- target_critter:

  Character. Name of the species whose yield is maximised. Must match a
  name in `fauna`.

## Value

nothing at the moment

When `opt = TRUE`: a negative numeric scalar (negative total yield for
the target critter, for minimisation). When `opt = FALSE`: a
[`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md)
output list at the given effort level.

## See also

[`assign_ref_points`](https://danovando.github.io/marlin/reference/assign_ref_points.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md)

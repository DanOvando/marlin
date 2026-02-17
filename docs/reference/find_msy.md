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

  the multiplier on the base effort in fleet (assumes fleet model is
  already tuned)

- fauna:

  a list of critters

- fleets:

  a list of fleets

- opt:

  TRUE = optimize, FALSE = return MSY conditions

- target_critter:

  the name of the critter to find MSY for

## Value

nothing at the moment

When `opt = TRUE`: a negative numeric scalar (negative total yield for
the target critter, for minimisation). When `opt = FALSE`: a
[`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md)
output list at the given effort level.

MSY conditions

## See also

[`assign_ref_points`](https://danovando.github.io/marlin/reference/assign_ref_points.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md)

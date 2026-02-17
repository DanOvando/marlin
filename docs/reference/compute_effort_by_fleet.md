# Compute effective effort per fleet

For each fleet, calculates the share of `base_effort` that overlaps with
viable habitat (`b0_p > 0`) within the fleet's fishing grounds. Effort
is a fleet-level quantity (not species-specific): a fleet exerts one
pool of effort that affects each species differently through
catchability.

## Usage

``` r
compute_effort_by_fleet(fauna_species, fleets)
```

## Arguments

- fauna_species:

  A single fauna element (used only for `b0_p`).

- fleets:

  A fleet list (output of
  [`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md)).

## Value

Named numeric vector of effective effort per fleet.

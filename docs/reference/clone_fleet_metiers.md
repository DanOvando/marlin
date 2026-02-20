# Deep-clone all metiers in a fleet list

Creates independent copies of all metier R6 objects so that
modifications during tuning do not mutate the caller's fleet objects.

## Usage

``` r
clone_fleet_metiers(fleets)
```

## Arguments

- fleets:

  A fleet list (output of
  [`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md)).

## Value

A deep-cloned copy of `fleets` with independent metier objects.

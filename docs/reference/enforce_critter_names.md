# Enforce Critter Names on Storage Steps

Ensures that a storage step list is named and ordered consistently with
the fauna names used elsewhere in the simulation.

## Usage

``` r
enforce_critter_names(x, fauni)
```

## Arguments

- x:

  List representing one time-step of simulation storage.

- fauni:

  Character vector of expected fauna names.

## Value

The input list `x`, named and reordered to match `fauni`.

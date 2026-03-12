# Precompute Baranov Catch Equation Inputs

Assembles the per-species selectivity, mortality, biomass, and price
matrices needed by the IFD effort-allocation solver.

## Usage

``` r
precompute_baranov_inputs(storage, fauna, fleets, target_fleet, E_exo, P)
```

## Arguments

- storage:

  List of per-species population state (one element per species), each
  containing a `b_p_a` matrix (patches x ages).

- fauna:

  List of fauna objects (e.g. from
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md)).

- fleets:

  List of fleet objects (e.g. from
  [`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md)).

- target_fleet:

  Integer index of the fleet being optimised.

- E_exo:

  List of effort vectors for all other fleets (one per fleet), or `NULL`
  if no exogenous effort.

- P:

  Integer number of spatial patches.

## Value

A list with elements `alpha_mats`, `other_mort_mats`, `biomass_mats`
(each a list of matrices, one per species) and `price_s` (numeric vector
of prices).

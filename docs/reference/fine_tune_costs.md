# Fine-Tune Cost Parameters for Open-Access Fleets

Objective function used by
[`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md)
to calibrate `cost_per_unit_effort` for open-access fleets so that the
simulated depletion or exploitation rate matches a user-supplied target.

## Usage

``` r
fine_tune_costs(log_cost_mult, target, fauna, fleets, years = 25, tune_type)
```

## Arguments

- log_cost_mult:

  Numeric vector. Log-scale cost-per-unit-effort multipliers, one per
  open-access fleet.

- target:

  Data frame with columns `critter` and `target`, giving the desired
  depletion level for each species.

- fauna:

  List of fauna objects (from
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md)).

- fleets:

  List of fleet objects (from
  [`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md)).

- years:

  Integer. Number of years to simulate during optimisation (default 25).

- tune_type:

  Character. One of `"depletion"` or `"explt"`, indicating what metric
  is being matched to `target`.

## Value

Numeric scalar: sum of squared differences between simulated and target
values (used as the objective for optimisation).

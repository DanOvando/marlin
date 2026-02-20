# Tune fleet parameters to match targets

Convenience wrapper to adjust one or more fleet parameters so simulated
quantities (for example catch, effort, or revenue summaries) better
match user-specified targets.

## Usage

``` r
tune_fleets(
  fauna,
  fleets,
  years = 50,
  tune_type = "f",
  tune_costs = TRUE,
  depl_tol = 0.025
)
```

## Arguments

- fleets:

  Named list of fleet objects to tune.

- ...:

  Additional arguments controlling what to tune, target values, and how
  the tuning is performed (see function body for supported options).

## Value

A list containing the updated fleet objects and any tuning diagnostics
produced by the routine.

## Details

This function is intended for quick calibration workflows. It evaluates
model output under the current fleet settings, updates the requested
parameters using the chosen tuning rule, and returns the updated fleet
objects along with diagnostics.

## See also

[`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)

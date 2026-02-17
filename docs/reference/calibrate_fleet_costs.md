# Calibrate cost parameters from equilibrium simulation output

Extracts revenue and effort from an equilibrium time step, then solves
for `cost_per_unit_effort` and `effort_reference` for each fleet so that
costs match the target cost-to-revenue ratio.

## Usage

``` r
calibrate_fleet_costs(eq, fleets)
```

## Arguments

- eq:

  A single equilibrium time step from
  [`simmar`](https://danovando.github.io/marlin/reference/simmar.md)
  output (i.e. `sim[[length(sim)]]`).

- fleets:

  A fleet list (output of
  [`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md)).

## Value

A data frame with columns: `fleet`, `cost_per_unit_effort`,
`effort_reference`, `revenue`, `total_effort`, `n_open_patches`,
`cr_ratio`, `travel_weight`, `implied_cr`.

## Details

The calibration formula is: \$\$c_0 = \frac{\text{cr\\ratio} \times
R}{E^{ref} \times P \times (1 + \theta)}\$\$

The effort cost exponent \\\gamma\\ drops out because tuning runs under
constant effort, distributing effort roughly uniformly (\\E_l / E^{ref}
\approx 1\\, so \\1^\gamma = 1\\).

# Tune Fleet Catchability and Costs to Target Initial Conditions

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

- fauna:

  A named list of fauna objects from
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md).

- fleets:

  A named list of fleet objects from
  [`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md).

- years:

  Integer. Number of years to simulate during tuning runs. Longer values
  ensure equilibrium is reached. Default `50`.

- tune_type:

  Character. One of `"f"` / `"explt"` (tune to fishing mortality rate)
  or `"depletion"` (tune to B/B0). See Details.

- tune_costs:

  Logical. If `TRUE` (default), calibrate `cost_per_unit_effort` for
  each fleet so that equilibrium costs match the fleet's `cr_ratio`.

- depl_tol:

  Numeric. Relative tolerance for the depletion-tuning convergence
  check. A warning is issued if any species' achieved depletion differs
  from the target by more than this fraction. Default `0.025` (2.5\\

A tuned copy of `fleets` with updated `catchability`,
`spatial_catchability`, `vul_p_a`, and (if `tune_costs = TRUE`)
`cost_per_unit_effort` and `effort_reference` fields for each metier.
Adjusts the catchability coefficient of each fleet's metiers so that the
simulated fishery reaches either a target fishing mortality rate or a
target depletion level (B/B0). Optionally tunes `cost_per_unit_effort`
for each fleet to match a target cost-to-revenue ratio (`cr_ratio`).This
function should be called after
[`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md)
and before passing fleets to
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md); the
tuned fleet list is returned and should replace the original.

### Tuning type

- `"f"` (or `"explt"`):

  Sets catchability directly from each critter's `init_explt` and the
  metier's `p_explt` share. Analytical; no numerical solver required.
  Fast but does not guarantee a precise equilibrium depletion.

- `"depletion"`:

  Uses a two-phase Broyden solver via `nleqslv` to find catchabilities
  that produce the target depletion specified in each critter's
  `depletion` field. A logistic link function keeps catchabilities in
  (0, 1) throughout optimisation. Post-solve validation warns if
  achieved depletion differs from the target by more than `depl_tol`.

When `tune_costs = TRUE` and `tune_type = "depletion"`, a two-pass cost
calibration is used: an initial heuristic pass before the depletion
solver, then a refined pass using the actual equilibrium. This ensures
costs accurately reflect final equilibrium conditions.

Note that tuning is approximate: post-tuning values will not perfectly
match targets because some calibration steps depend on earlier steps.

[`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md)

# Tune Fleet Catchability and Costs to Match Targets

Calibrates fleet catchability coefficients (and optionally costs) so
that the simulated fishery matches user-specified exploitation rates or
depletion levels. This should be called after
[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md)
and
[`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md),
and before
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md).

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

  Named list of critter objects from
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md).

- fleets:

  Named list of fleet objects from
  [`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md).

- years:

  Integer. Number of years to simulate when running to equilibrium for
  calibration. Higher values give more accurate equilibria but take
  longer. Default `50`.

- tune_type:

  Character. Tuning mode: `"f"` or `"explt"` (equivalent) to match
  exploitation rates, or `"depletion"` to match SSB/SSB0 targets.
  Default `"f"`.

- tune_costs:

  Logical. If `TRUE` (default), calibrates `cost_per_unit_effort` for
  each fleet to match its `cr_ratio`. Set to `FALSE` to skip cost
  calibration (e.g. for `constant_effort` fleets where costs are
  irrelevant).

- depl_tol:

  Numeric. Relative tolerance for depletion-mode validation. A warning
  is issued if the achieved depletion for any species differs from its
  target by more than this fraction. Default `0.025` (2.5 percent).

## Value

A named list of tuned fleet objects (same structure as the input
`fleets`), with updated catchability, spatial catchability, and (if
`tune_costs = TRUE`) cost parameters. Original fleet models are
preserved.

## Details

### Tuning modes

- `"f"` / `"explt"`:

  Sets each fleet's catchability so that the total instantaneous fishing
  mortality matches each species' `init_explt` (from
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md)).
  Catchability is computed analytically from `init_explt`, `p_explt`
  (each metier's share), and effective effort. Fast (one simulation for
  cost calibration only).

- `"depletion"`:

  Uses a numerical solver
  ([`nleqslv`](https://rdrr.io/pkg/nleqslv/man/nleqslv.html), Broyden
  method) to find fishing mortalities that drive each species to its
  target `depletion` (SSB/SSB0, set in
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md)).
  Slower (multiple simulations). Warns if achieved depletion differs
  from the target by more than `depl_tol`.

### What gets calibrated

- **Catchability**: Each metier's `catchability` and
  `spatial_catchability` are set so that fishing mortality on each
  species matches the target. The share of total F allocated to each
  fleet is controlled by the metier's `p_explt` values (normalised to
  sum to 1 across fleets for each species).

- **Costs** (when `tune_costs = TRUE`): Runs an equilibrium simulation
  and calibrates `cost_per_unit_effort` and `effort_reference` for each
  fleet so that the cost-to-revenue ratio matches the fleet's
  `cr_ratio`. Relevant for `open_access` and `sole_owner` fleet models.

### Workflow

1.  Deep-clones fleet metiers to avoid mutating the caller's objects.

2.  Temporarily sets all fleet models to `"constant_effort"`.

3.  Normalises `p_explt` across fleets for each species.

4.  Sets initial catchability estimates.

5.  Runs equilibrium simulation(s) and calibrates costs.

6.  (Depletion mode only) Solves for catchability via `nleqslv` and
    validates achieved depletion.

7.  Restores original fleet models and returns tuned fleets.

## See also

[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md),
[`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`Metier`](https://danovando.github.io/marlin/reference/Metier.md)

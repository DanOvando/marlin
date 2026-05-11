# Create a Fleet Object

Builds a fleet object specifying the fishing behaviour, spatial
dynamics, and economic parameters for one group of vessels targeting one
or more species. Always run
[`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md)
after creating fleets to calibrate catchability and costs to the desired
initial conditions.

## Usage

``` r
create_fleet(
  metiers,
  mpa_response = "stay",
  fleet_model = "constant_effort",
  oa_max_growth_per_year = 0.5,
  oa_max_decline_per_year = 0.5,
  oa_signal_half = 0.3,
  cost_per_unit_effort = 1,
  spatial_allocation = "rpue",
  effort_cost_exponent = 1.2,
  travel_fraction = 0,
  ports = NULL,
  cost_per_distance = 1,
  cr_ratio = 1,
  resolution,
  patch_area = 1,
  base_effort = NULL,
  fishing_grounds = NULL,
  responsiveness = 0.025,
  memory_halflife = 0
)
```

## Arguments

- metiers:

  Named list of
  [`Metier`](https://danovando.github.io/marlin/reference/Metier.md) R6
  objects, one per species in `fauna`. Each metier specifies price,
  selectivity, and relative catchability for that fleet-species
  combination. Names must match species names in `fauna`.

- mpa_response:

  Character. Vessel response to MPA closures: `"stay"` (vessels stay in
  the closed area) or `"leave"` (vessels redistribute to open patches,
  concentrating effort).

- fleet_model:

  Character. Effort dynamics model; see Details. One of
  `"constant_effort"` (default), `"open_access"`, `"sole_owner"`, or
  `"manual"`.

- oa_max_growth_per_year:

  Numeric. Maximum fractional increase in total effort per year under
  highly profitable conditions for open-access and sole-owner fleets
  (e.g. `0.5` = +50% per year). Converted to a per-step multiplier
  internally. Must be \> 0.

- oa_max_decline_per_year:

  Numeric in (0, 1). Maximum fractional decrease in total effort per
  year under highly unprofitable conditions (e.g. `0.5` = -50% per
  year). Must be in (0, 1).

- oa_signal_half:

  Numeric in (0, 1). Profitability signal value at which effort
  adjustment speed reaches half its maximum. Lower values make fleets
  more sensitive. Typical range `0.15`–`0.4`. Default `0.3`.

- cost_per_unit_effort:

  Numeric. Base cost per unit of effort. Overridden by
  [`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md)
  when `tune_costs = TRUE`; rarely needs manual adjustment.

- spatial_allocation:

  Character. Spatial effort allocation strategy; see Details. Default
  `"rpue"`.

- effort_cost_exponent:

  Numeric. Exponent \\\gamma\\ in the effort cost function, controlling
  congestion / convexity. Values \> 1 make additional units of effort
  progressively more expensive. Default `1.2`.

- travel_fraction:

  Numeric in \[0, 1). Fraction of total costs attributable to travel at
  equilibrium. `0` means no spatial cost heterogeneity (all patches
  equally costly). Controls how strongly port proximity shapes the
  spatial cost surface.

- ports:

  Data frame with columns `x` and `y` giving port patch coordinates.
  Minimum distances from each patch to the nearest port are used to
  compute the travel-cost component when `travel_fraction > 0`.

- cost_per_distance:

  Numeric. Deprecated; use `travel_fraction` instead.

- cr_ratio:

  Numeric. Target cost-to-revenue ratio at equilibrium. `1` implies zero
  profits (open-access equilibrium). Used by
  [`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md)
  to calibrate `cost_per_unit_effort`.

- resolution:

  Integer scalar or length-2 integer vector `c(nx, ny)`. Must match the
  resolution of the `fauna` objects.

- patch_area:

  Numeric. Area of each patch (km^2). Used to compute port distances.

- base_effort:

  Numeric. Total effort units available to the fleet. Defaults to
  `prod(resolution)` (one unit per patch). Catchability is calibrated
  relative to this value by
  [`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md).

- fishing_grounds:

  Data frame with columns `x`, `y`, and `fishing_ground` (logical or
  numeric). Restricts where effort can be deployed; `TRUE`/`1` = open,
  `FALSE`/`0` = closed. When `spatial_allocation = "manual"`, numeric
  values in `fishing_ground` are used as effort weights. Defaults to all
  patches open.

- responsiveness:

  Numeric. Per-step responsiveness of patch effort to the objective
  signal in
  [`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md)
  (the \\\eta\\ parameter of the multiplicative update). Larger values
  move effort more aggressively toward high-objective patches each step;
  if too large, can drive period-2 sawtooth oscillation. Default
  `0.025`.

- memory_halflife:

  Non-negative numeric. Half-life **in years** of the EWMA applied to
  this fleet's spatial objective surface inside
  [`simmar`](https://danovando.github.io/marlin/reference/simmar.md).
  `0` (default) disables smoothing — the fleet sees only the previous
  step's objective, matching legacy behavior. The parameter is
  season-agnostic: simmar converts it internally to time steps
  (`halflife_steps = halflife * steps_per_year`) so a given value
  produces the same calendar-time smoothing regardless of how many
  seasons per year the model uses. Larger values dampen high-frequency
  feedback oscillations by blending in past objective surfaces; the
  per-step weight on the current surface is \\\alpha = 1 -
  0.5^{1/halflife\_{steps}}\\. Smoothing updates only patches that are
  currently open; closed patches retain their last-seen smoothed value
  ("freeze and resume"). Early post-bootstrap steps use a Welford-style
  ramp (effective \\\alpha\_\mathrm{eff} = \max(\alpha, 1/n)\\) so the
  smoothed surface is not anchored to the first observed buffet column.

  Practical guidance: values around `0.5`–`1.5` years are the typical
  sweet spot. Larger halflives introduce phase lag of roughly \\1.44
  \times \mathrm{halflife}\\ years between a true change in patch
  marginal value and the fleet's perceived value, which can produce
  low-frequency overshoot/undershoot — a different pathology from the
  high-frequency sawtooth that motivates the parameter. If catch
  trajectories under a given halflife show slow swings that aren't
  present at `halflife = 0`, reduce it.

## Value

A named list (fleet object) with all parameters needed by
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md) and
[`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md),
including computed travel weights, normalised cost-per-patch, and (for
open-access / sole-owner fleets) the derived annual effort-adjustment
parameters.

## Details

### Fleet models

- `"constant_effort"`:

  Total effort is fixed at `base_effort` each time step. Use for
  scenarios where fishing pressure is prescribed externally.

- `"open_access"`:

  Total effort adjusts each step based on a normalised
  average-profitability signal. Equilibrates where total profits = 0.
  Entry/exit speed is controlled by `oa_max_growth_per_year`,
  `oa_max_decline_per_year`, and `oa_signal_half`.

- `"sole_owner"`:

  Identical dynamics to `"open_access"` but uses the marginal (not
  average) profit signal. Equilibrates at MEY where marginal profit = 0.
  Requires
  [`calc_marginal_value()`](https://danovando.github.io/marlin/reference/calc_marginal_value.md)
  to be computed each step;
  [`simmar`](https://danovando.github.io/marlin/reference/simmar.md)
  handles this automatically.

- `"manual"`:

  Total effort each step is taken directly from a user-supplied vector;
  see the `manager` argument of
  [`simmar`](https://danovando.github.io/marlin/reference/simmar.md).

### Spatial allocation

The `spatial_allocation` argument determines how total fleet effort is
distributed among patches each step. Options:

- `"rpue"`:

  Revenue per unit effort (default). Effort concentrates in high-revenue
  patches.

- `"revenue"`:

  Total revenue. Similar to `"rpue"` but favours larger patches.

- `"ppue"`:

  Profit per unit effort (cost-aware).

- `"profit"`:

  Total profit (cost-aware).

- `"cpue"` / `"catch"`:

  Catch-based variants.

- `"marginal_revenue"` / `"marginal_profit"`:

  Uses finite- difference marginal returns from
  [`calc_marginal_value`](https://danovando.github.io/marlin/reference/calc_marginal_value.md).
  Required for `"sole_owner"` fleets. Requires
  `fleet_model = "sole_owner"` or explicit pre-computation.

- `"manual"`:

  Effort distributed proportionally to continuous weights in
  `fishing_grounds$fishing_ground` (0–1 valued).

- `"uniform"`:

  Effort spread equally across all open patches.

### Costs

Total cost per fleet is: \$\$C = c_0 \\ E^{ref} \sum_l
\left\[\left(\frac{E_l}{E^{ref}}\right)^\gamma + \theta \\ \tilde{d}\_l
\frac{E_l}{E^{ref}}\right\]\$\$ where \\c_0\\ is `cost_per_unit_effort`,
\\E^{ref}\\ is the reference effort per patch, \\\gamma\\ is
`effort_cost_exponent`, \\\theta\\ is `travel_weight` (derived from
`travel_fraction`), and \\\tilde{d}\_l\\ is the normalised distance from
patch \\l\\ to the nearest port.

## See also

[`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md),
[`Metier`](https://danovando.github.io/marlin/reference/Metier.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a metier for a single species
met <- Metier$new(
  critter     = fauna[["tuna"]],
  price       = 10,
  sel_form    = "logistic",
  sel_start   = 0.3,
  sel_delta   = 0.1,
  catchability = 0.01,
  p_explt     = 1
)

# Constant-effort fleet
fleet <- create_fleet(
  metiers    = list(tuna = met),
  resolution = c(10, 10)
)

# Open-access fleet with port-based travel costs
ports <- data.frame(x = 1, y = 1)
oa_fleet <- create_fleet(
  metiers              = list(tuna = met),
  fleet_model          = "open_access",
  spatial_allocation   = "ppue",
  travel_fraction      = 0.3,
  ports                = ports,
  cr_ratio             = 0.9,
  resolution           = c(10, 10)
)

fleets <- list(fleet = fleet)
fleets <- tune_fleets(fauna, fleets, tune_type = "depletion")
} # }
```

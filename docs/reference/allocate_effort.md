# Allocate Fishing Effort Across Patches Using a Multiplicative Velocity-Field Update

## Usage

``` r
allocate_effort(
  effort_by_patch,
  total_effort_by_fleet = NULL,
  buffet,
  fleets,
  open_patch,
  clip_z = 5,
  scale = c("mad", "iqr", "sd"),
  eps_mix = 0,
  floor_frac = 0,
  flatness_tol = 0.001,
  min_scale_abs = 1e-10,
  adaptive_floor_pct = 0.01
)
```

## Arguments

- effort_by_patch:

  Numeric matrix of effort by patch (rows) and fleet (columns). Column
  names should match fleet names.

- total_effort_by_fleet:

  Numeric vector of total effort per fleet. When `NULL` (default),
  computed as `colSums(effort_by_patch)`. Provides a way to update total
  effort independently of the current spatial distribution.

- buffet:

  Named list returned by
  [`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md)
  with `output_format = "matrix"`. Must contain: `r_p_fl`, `c_p_fl`,
  `prof_p_fl`, `rpue_p_fl`, `cpue_p_fl`, `ppue_p_fl`. For
  `"marginal_profit"` / `"marginal_revenue"` fleets, must also contain
  `mp_p_fl` / `mr_p_fl` from
  [`calc_marginal_value`](https://danovando.github.io/marlin/reference/calc_marginal_value.md).

- fleets:

  Named list of fleet objects. Each fleet must have
  `$spatial_allocation` set to one of: `"revenue"`, `"catch"`,
  `"profit"`, `"rpue"`, `"cpue"`, `"ppue"`, `"marginal_revenue"`,
  `"marginal_profit"`, `"manual"`, or `"uniform"`.

- open_patch:

  Logical matrix (patches x fleets) or logical vector (shared across all
  fleets) indicating which patches are currently open. Closed patches
  receive zero effort. Allows fleet-specific closures (e.g. different
  fishing grounds per fleet).

- clip_z:

  Numeric. Maximum absolute standardised objective value; prevents
  extreme updates from outliers. Default `5`.

- scale:

  Character. Method to scale the objective before computing the velocity
  field: `"mad"` (median absolute deviation; default, most robust to
  outliers), `"iqr"`, or `"sd"`.

- eps_mix:

  Numeric in 0, 1. Exploration mixing fraction. When \> 0, blends the
  optimal allocation with a uniform distribution to allow re-entry into
  temporarily abandoned patches. Default `0`.

- floor_frac:

  Numeric. Fractional effort floor relative to the mean open-patch
  effort. Prevents complete abandonment of any patch. Default `0` (no
  floor).

- flatness_tol:

  Numeric. Range-based CV threshold for detecting a flat objective
  (range / \|median\|). Below this threshold, effort is distributed
  uniformly. Default `1e-3` (0.1\\

  min_scale_absNumeric. Absolute minimum scale value to prevent division
  by near-zero numbers. Default `1e-10`.

  adaptive_floor_pctNumeric. Percentage of the median objective used as
  an adaptive minimum scale floor. Default `0.01` (1\\ A named list:

  `effort_new`

  :   Numeric matrix of updated effort (patches x fleets), with column
      names preserved.

  `velocity_by_patch`

  :   Numeric matrix of velocity fields (patches x fleets).

  `z_by_patch`

  :   Numeric matrix of standardised objective values (patches x
      fleets).

  `flat_objective`

  :   Named logical vector (one per fleet): `TRUE` if the objective was
      flat for that fleet (effort distributed uniformly). `NA` for
      manual/uniform fleets.

  `cv`

  :   Named numeric vector: range-based CV of the objective across open
      patches (used for flatness detection).

  `obj_range`

  :   Named numeric vector: range of the objective across open patches.

  Updates patch-level effort for one or more fleets using a
  replicator-style gradient-flow heuristic derived from an objective
  signal (the "buffet") returned by
  [`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md).
  Each fleet's objective is selected automatically from the buffet based
  on its `spatial_allocation` setting.Called automatically by
  [`simmar`](https://danovando.github.io/marlin/reference/simmar.md)
  each time step; can also be used directly for custom allocation
  scenarios.

  ### Algorithm

  For each fleet:

  1.  **Flatness detection**: If the objective is nearly uniform across
      patches (range-based CV \< `flatness_tol`), effort is distributed
      uniformly across open patches. This avoids spurious edge
      concentration when there is no meaningful gradient.

  2.  **Standardise**: The objective is centred (median-subtracted) and
      scaled (by MAD, IQR, or SD) to produce a \\z\\ score, clipped to
      `[-clip_z, clip_z]`.

  3.  **Velocity field**: The centered \\z\\ score becomes a velocity
      \\v = z - \bar{z}\\, where \\\bar{z}\\ is the mean across open
      patches.

  4.  **Multiplicative update**: \$\$e\_{new,p} \propto e\_{old,p}
      \exp(\eta \\ v_p)\$\$ Closed patches receive zero effort. Total
      effort is conserved.

  ### Special cases

  `spatial_allocation = "manual"`

  :   Effort is distributed proportionally to the continuous weights in
      `fishing_grounds$fishing_ground`. The optimizer is bypassed
      entirely; `buffet` values are not used for this fleet.

  `spatial_allocation = "uniform"`

  :   Effort is spread equally across all open patches regardless of the
      objective or weights.

  MPA closures

  :   Patches flagged as closed in `open_patch` receive zero effort.
      When a fleet's `mpa_response = "leave"`, the concentrator vector
      in
      [`simmar`](https://danovando.github.io/marlin/reference/simmar.md)
      zeros out MPA patches before calling this function.

  [`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md),
  [`calc_marginal_value`](https://danovando.github.io/marlin/reference/calc_marginal_value.md),
  [`simmar`](https://danovando.github.io/marlin/reference/simmar.md)

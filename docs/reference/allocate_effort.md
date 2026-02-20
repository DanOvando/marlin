# Allocate fishing effort across patches

Updates patch-level effort for one or more fleets based on an objective
signal contained in the `buffet` output from
[`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md).

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

  Numeric matrix with patches in rows and fleets in columns.

- total_effort_by_fleet:

  Optional numeric vector of total effort per fleet. If `NULL`, totals
  are computed from `effort_by_patch`.

- buffet:

  Named list from
  [`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md)
  with matrix outputs (for example `r_p_fl`, `c_p_fl`, `prof_p_fl`, and
  per-unit variants).

- fleets:

  Named list of fleet objects (each with `$spatial_allocation`).

- open_patch:

  Logical matrix (patches x fleets) indicating where each fleet may
  fish.

- clip_z:

  Numeric scalar controlling how strongly extreme objective values
  affect the update.

- scale:

  Character; method used to scale the objective before updating effort.

- eps_mix:

  Numeric between 0 and 1; mix-in fraction for uniform exploration.

- floor_frac:

  Numeric; optional effort floor applied across open patches.

- flatness_tol:

  Numeric; threshold for treating the objective as flat across patches.

- min_scale_abs:

  Numeric; minimum scale to avoid division by near-zero values.

- adaptive_floor_pct:

  Numeric; fraction of the median objective used as an adaptive scale
  floor.

## Value

A named list containing the updated effort matrix (`effort_new`) and
additional diagnostic matrices/vectors used to track the update (for
example standardised objectives and flatness indicators).

## Details

This helper is called by
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md) each
time step. For each fleet, it converts the selected objective (e.g.
revenue, catch, profit, or a per-unit variant) into a smooth
multiplicative update so that higher-objective patches receive more
effort, subject to `open_patch`.

## See also

[`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md),
[`calc_marginal_value`](https://danovando.github.io/marlin/reference/calc_marginal_value.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)

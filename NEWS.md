# marlin 1.1.4

## Bug Fixes

- Catch is now weighed at the age the fish were caught at. Since the 1.0.0 rewrite of `sim_fish()`, the Baranov catch of age *a* fish (numbers, *F* and *Z* all at age *a*) was converted to biomass with `weight_at_age[a + 1]`, the next age's weight, while the pre-1.0 code used the same age. This inflated catch biomass, revenue and profit by `weight_at_age[a + 1] / weight_at_age[a]`, most at the youngest ages, and through revenue and profit it also affected spatial effort allocation, open-access and sole-owner effort dynamics, and `tune_fleets()`'s cost calibration. Survival, numbers, biomass and spawning biomass were not affected directly. The plus group was always weighed correctly. A new test reproduces one step's catch by hand.

## New Features

- Catch is now also returned in numbers. `sim_fish()` returns `cn_p_a` alongside `c_p_a` (weight); `allocate_yields()` splits it across fleets as `cn_p_a_fl`; `simmar()` stores both, shifted to the step the catch came from just like `c_p_a`; and `process_marlin()` adds a `cn` column to `$fauna` and a `catch_numbers` column to `$fleets`. Useful for simulating sampled fish. Existing outputs are unchanged.

## Documentation

- `c_p_a`, `c_p_a_fl` and the processed `c` / `catch` columns are catch in weight (biomass). The docs for `simmar()`, `process_marlin()`, `allocate_yields()` and `plot_marlin()` had described them as catch in numbers, which invites weighting catch a second time.

# marlin 1.1.3

## Bug Fixes

- Fleet costs for multi-species fleets with a convex effort cost (`effort_cost_exponent` > 1) now match the documented cost function. `allocate_yields()` runs once per species and the per-species costs are summed, and through 1.1.2 each species' share was computed by dividing patch effort by the number of species *before* the congestion exponent. That sums back to the fleet's cost only when the exponent is 1; for an exponent of gamma it scaled the congestion term by n_species^(1 - gamma) (0.66 for four species at gamma = 1.3), so multi-species fleets paid less congestion cost than the formula and `tune_fleets()`'s cost calibration assumed. Each species now carries an equal share of the cost computed on the fleet's total effort. Single-species fleets and `effort_cost_exponent = 1` are unaffected; multi-species simulations with a convex cost will see higher costs and, for open-access and sole-owner fleets, lower equilibrium effort.

# marlin 1.1.1

## Behavior Changes

- The default youngest age (`min_age`) tracked by `create_critter()` / `Fish$new()` is now one time step (`1 / seasons`) rather than `0`. Because recruits are deposited at the youngest age *after* the population is aged, the youngest cohort is already one time step removed from the spawning pool, so age 0 ("at birth") was never really represented; starting at one time step is more consistent with that, keeps the extrapolated (below the Lorenzen fitting range) age-0 natural mortality out of the model, and makes per-recruit quantities comparable with age-1-indexed assessment frameworks. Pass `min_age = 0` to recover the previous behavior. Existing simulations will see a shifted age vector and, with default Lorenzen M, somewhat different absolute per-recruit quantities (this does not move SPR- or F_msy-type reference points).

- Lorenzen natural mortality (`lorenzen_m = TRUE`, the default) is now anchored so that the **mean** M over the adult age range equals the supplied `m`, rather than pinning `m` at the asymptotic (largest/oldest) size. A single reported `m` is therefore interpreted as an average adult rate — consistent with how `tmax`- and growth-based empirical M estimators are derived, and analogous to the age-range Lorenzen option in Stock Synthesis (Powers 2014, ICES J. Mar. Sci. 71:1629). The new `lorenzen_ref_age` argument to `create_critter()` / `Fish$new()` sets the lower age bound of that range (defaults to the age at maturity). Existing simulations using the default Lorenzen M will produce slightly different M-at-age and, in turn, unfished equilibria and tuned catchabilities.

## Bug Fixes

- `Fish$plot_movement()` no longer errors on its first call in a fresh session (e.g. on a critter restored via `readRDS()`). `Matrix` is now a declared dependency, so its S4 methods are registered when marlin loads, and the plot-time multiply densifies the (small) movement matrices to avoid relying on those methods at all.

# marlin 1.1.0

# marlin 1.0.1

# marlin 1.0.0

## Breaking Changes

- Travel costs are now set by `travel_fraction`, which species the proportion of the total costs that come from travel. cost_per_distance is now deprecated and may cause breaks

- Old open access parameters no longer work. see `?create_fleet` for new parameters

## Improvements

- Movement is now dealt with via sparse matrices, vastly increasing speed and memory efficiency in most cases

- `process_marlin` now works vis `data.table`, is orders of magnitude faster

- fleet spatial allocation strategies are now much more stable and consistent

## New Capabilities

- Added calculation of marginal revenues and marginal profits in space and time by fleet, allowing for spatial allocation
models `marginal_profit` and `marginal_revenue`, as well as fleet_model `sole_owner`

- Added flexible patch areas

# marlin 0.8.0

# marlin 0.7.0

# marlin 0.6.1

# marlin 0.6.0

# marlin 0.5.0

# marlin 0.2.0

# marlin 1.0.3

# marlin 1.0.2

# marlin 1.0.1

# marlin 1.0.0

# marlin 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.

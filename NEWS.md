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

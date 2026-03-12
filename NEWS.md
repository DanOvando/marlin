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

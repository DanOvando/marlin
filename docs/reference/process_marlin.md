# Tidy Simulation Output from simmar

Converts the nested list output of
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md) into
two tidy tibbles: one for population state (`$fauna`) and one for fleet
outcomes (`$fleets`). This is the primary post-processing step before
plotting or analysis.

## Usage

``` r
process_marlin(sim, steps_to_keep = NULL, time_step = NULL, keep_age = TRUE)
```

## Arguments

- sim:

  Named list returned by
  [`simmar`](https://danovando.github.io/marlin/reference/simmar.md).
  Step names should follow the `"year_season"` convention (with optional
  `"step_"` prefix, which is automatically stripped).

- steps_to_keep:

  Character or integer vector of step names / indices to include in the
  output. Useful for reducing memory use when only the final years are
  needed. Default `NULL` keeps all steps.

- time_step:

  Numeric. Fraction of a year per step (e.g. `0.25` for quarterly
  seasons). When `NULL`, inferred from the step names in `sim` (requires
  at least two steps).

- keep_age:

  Logical. If `TRUE` (default), return age-structured output. If
  `FALSE`, aggregate across all ages before returning.

## Value

A named list with two elements:

- `fauna`:

  Tibble of population state; see Details.

- `fleets`:

  Tibble of fleet outcomes; see Details.

## Details

### Output columns

**`$fauna`** — one row per critter × patch × age × step:

- `critter`: species name

- `patch`: patch index

- `x`, `y`: spatial coordinates

- `age`: age class (numeric)

- `mean_length`: mean length at age (from life history)

- `step`: decimal year (e.g. 5.75 = year 5, season 4 of 4)

- `n`: numbers at age

- `b`: biomass at age

- `ssb`: spawning stock biomass at age

- `c`: catch in numbers at age

**`$fleets`** — one row per critter × patch × age × fleet × step:

- `critter`, `patch`, `x`, `y`, `age`, `fleet`, `step`, `mean_length`:
  as above

- `catch`: catch in numbers

- `revenue`: revenue

- `effort`: effort units

- `cpue`: catch per unit effort (`NA` where effort = 0)

When `keep_age = FALSE`, age classes are summed and a single
`age = "all"` row is returned per spatial unit per step.

## See also

[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`plot_marlin`](https://danovando.github.io/marlin/reference/plot_marlin.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simmar(fauna = fauna, fleets = fleets, years = 50)

# Full output (all steps, all ages)
proc <- process_marlin(sim, time_step = 1)

# Last 10 years only, aggregated across ages
last_steps <- tail(names(sim), 10)
proc_last  <- process_marlin(sim,
                             steps_to_keep = last_steps,
                             time_step     = 1,
                             keep_age      = FALSE)

# Plot SSB over time
plot_marlin(proc, plot_var = "ssb", plot_type = "time")
} # }
```

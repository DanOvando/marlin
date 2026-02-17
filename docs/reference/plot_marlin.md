# Plot marlin Simulation Results

Generates time series, spatial, length composition, or age composition
plots from one or more
[`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md)
outputs. Multiple scenarios can be overlaid by passing them as named
arguments.

## Usage

``` r
plot_marlin(
  ...,
  steps_to_plot = NA,
  plot_var = "ssb",
  plot_type = "time",
  fauna = NULL,
  drop_recruits = TRUE,
  plots = "fauna",
  max_scale = TRUE
)
```

## Arguments

- ...:

  One or more named
  [`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md)
  outputs. Names appear in the plot legend. If unnamed, letters a, b, c,
  ... are used. Each argument should be the full list returned by
  `process_marlin` (i.e. with `$fauna` and `$fleets` elements).

- steps_to_plot:

  Numeric or character vector. Which steps to include. Default `NA` uses
  all available steps.

- plot_var:

  Character. Quantity to plot. One of:

  `"ssb"`

  :   Spawning stock biomass (default)

  `"b"`

  :   Total biomass

  `"c"`

  :   Catch in numbers

  `"n"`

  :   Abundance in numbers

- plot_type:

  Character. Type of plot. One of `"time"` (default), `"space"`,
  `"length"`, or `"age"`.

- fauna:

  A `fauna` list from
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md).
  Required only when `plot_type = "length"` (used to access the
  length-at-age key).

- drop_recruits:

  Logical. If `TRUE` (default), drops the youngest 10\\ length classes
  that can dominate the axis scale.

- plots:

  Character. One of `"fauna"` (default; plots population state from
  `$fauna`) or `"fleets"` (plots fleet outcomes from `$fleets`).

- max_scale:

  Logical. If `TRUE` (default for `"time"` and `"space"`), normalises
  each series by its maximum so that values are expressed as a
  proportion of the maximum observed value.

## Value

A `ggplot2` object. Print it to render the plot, or save with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Details

### Plot types

- `"time"`:

  Line chart of the chosen `plot_var` summed across patches and ages
  over time, faceted by species. When `max_scale = TRUE`, each series is
  normalised by its own maximum so that all fits share a 0, 1 scale and
  cross-species comparisons are straightforward.

- `"space"`:

  Tile map of `plot_var` at the last step in `steps_to_plot`, summed
  across ages, faceted by species and fit.

- `"length"`:

  Density plot of length composition at specified steps, derived by
  projecting age-specific abundances through the critter's length-at-age
  key. Requires `fauna`. A maximum of 10 steps can be plotted
  simultaneously.

- `"age"`:

  Density plot of age composition at the last step in `steps_to_plot`.

## See also

[`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`theme_marlin`](https://danovando.github.io/marlin/reference/theme_marlin.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare two scenarios
proc_base <- process_marlin(sim_base)
proc_mpa  <- process_marlin(sim_mpa)

# Time series of SSB
plot_marlin(`No MPA` = proc_base, `With MPA` = proc_mpa,
            plot_var = "ssb", plot_type = "time")

# Spatial biomass map at final time step
plot_marlin(proc_mpa, plot_var = "b", plot_type = "space")

# Length composition (requires fauna)
plot_marlin(proc_base, plot_type = "length", fauna = fauna)
} # }
```

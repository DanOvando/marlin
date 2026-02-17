# Exploratory Fishing Step (No Population Effect)

Applies a given matrix of effort by patch and fleet to determine the
potential returns from fishing without modifying the underlying
population. This is an exploratory step used internally by
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md) to
build the "buffet" that
[`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md)
uses for spatial effort reallocation.

## Usage

``` r
go_fish(
  e_p_fl,
  fauna,
  n_p_a,
  fleets,
  groupers = c("fleet", "patch"),
  output_format = c("tidy", "matrix")
)
```

## Arguments

- e_p_fl:

  Numeric matrix of effort by patch (rows) and fleet (columns), with
  column names matching fleet names. Typically the current step's effort
  matrix from
  [`simmar`](https://danovando.github.io/marlin/reference/simmar.md).

- fauna:

  Named list of fauna objects from
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md).

- n_p_a:

  Named list (one element per species in `fauna`) of numbers-at-age
  matrices `[patch, age]`, representing the current population state to
  fish against. Typically extracted from `sim[[step]][[critter]]$n_p_a`.

- fleets:

  Named list of fleet objects from
  [`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md).

- groupers:

  Character vector. Grouping columns for tidy output format (default:
  `c("fleet", "patch")`). Ignored when `output_format = "matrix"`.

- output_format:

  Character. Output format:

  `"tidy"`

  :   (default) Returns a list of tidy data frames with species summed
      across and effort joined.

  `"matrix"`

  :   Returns a list of patches x fleets matrices (species summed).
      Designed for direct use with
      [`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md).
      Substantially faster.

## Value

A named list with six elements, each summed across species:

- `c_p_fl`:

  Catch by patch and fleet

- `r_p_fl`:

  Revenue by patch and fleet

- `prof_p_fl`:

  Profit by patch and fleet

- `cpue_p_fl`:

  Catch per unit effort (`NA` where effort = 0)

- `rpue_p_fl`:

  Revenue per unit effort (`NA` where effort = 0)

- `ppue_p_fl`:

  Profit per unit effort (`NA` where effort = 0)

Format of each element depends on `output_format`.

## Details

`go_fish` runs a single, movement-free time step (`move_fish = 0`) for
each critter in `fauna` using the supplied population state `n_p_a`,
then aggregates catch, revenue, and profit across species via
[`aggregate_yields`](https://danovando.github.io/marlin/reference/aggregate_yields.md).

The function is also useful for "what-if" analyses: you can probe the
revenue surface under different effort distributions before committing
to a simulation run.

## See also

[`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md),
[`calc_marginal_value`](https://danovando.github.io/marlin/reference/calc_marginal_value.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`aggregate_yields`](https://danovando.github.io/marlin/reference/aggregate_yields.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Probe the revenue surface under current effort
n_p_a_now <- lapply(fauna, function(cr) cr$n_p_a_0)
buffet <- go_fish(e_p_fl, fauna, n_p_a_now, fleets, output_format = "matrix")

# Use with allocate_effort to redistribute effort spatially
result <- allocate_effort(
  effort_by_patch = e_p_fl,
  buffet          = buffet,
  fleets          = fleets,
  open_patch      = open_patches
)
} # }
```

# Aggregate Per-Species Yield Outputs Into Fleet-Level Totals

Sums catch, revenue, and profit across species from the outputs of
[`allocate_yields`](https://danovando.github.io/marlin/reference/allocate_yields.md)
to produce a fleet-level "buffet" of spatial yield metrics. Also
computes per-unit-effort (PUE) counterparts. This is the final step that
builds the buffet consumed by
[`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md).

## Usage

``` r
aggregate_yields(
  yields,
  e_p_fl,
  output_format = c("matrix", "tidy"),
  groupers = c("fleet", "patch")
)
```

## Arguments

- yields:

  Named list of per-species yield outputs (one element per species).
  Each element must be the output of
  [`allocate_yields`](https://danovando.github.io/marlin/reference/allocate_yields.md)
  and contain patches x fleets matrices: `r_p_fl`, `c_p_fl`, and
  `prof_p_fl`. Names should match species names in `fauna`.

- e_p_fl:

  Numeric matrix of effort by patch (rows) and fleet (columns), with
  fleet names as column names. Used to compute per-unit-effort metrics.

- output_format:

  Character. Output format:

  `"matrix"`

  :   (default) Returns list of patches x fleets matrices with fleet
      names as column names. Fast; designed for internal use in
      [`simmar`](https://danovando.github.io/marlin/reference/simmar.md).

  `"tidy"`

  :   Returns list of tidy data frames with species summed and effort
      joined. Slower; used inside
      [`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md)
      for user-facing output.

- groupers:

  Character vector. Grouping columns for tidy format (default:
  `c("fleet", "patch")`). Ignored when `output_format = "matrix"`.

## Value

A named list with six elements (summed across all species):

- `r_p_fl`:

  Revenue by patch and fleet

- `c_p_fl`:

  Catch by patch and fleet

- `prof_p_fl`:

  Profit by patch and fleet

- `rpue_p_fl`:

  Revenue per unit effort (`NA` where effort = 0)

- `cpue_p_fl`:

  Catch per unit effort (`NA` where effort = 0)

- `ppue_p_fl`:

  Profit per unit effort (`NA` where effort = 0)

## Details

Used in two contexts:

1.  Inside
    [`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md),
    to aggregate exploratory fishing results across species before
    returning the buffet.

2.  Inside the
    [`simmar`](https://danovando.github.io/marlin/reference/simmar.md)
    step loop, after
    [`allocate_yields`](https://danovando.github.io/marlin/reference/allocate_yields.md)
    has been called for each species.

## See also

[`allocate_yields`](https://danovando.github.io/marlin/reference/allocate_yields.md),
[`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Inside simmar after the species loop:
buffet <- aggregate_yields(yields_this_step, updated_e_p_f)

# Inside go_fish with tidy output:
buffet <- aggregate_yields(yields, e_p_fl, output_format = "tidy")
} # }
```

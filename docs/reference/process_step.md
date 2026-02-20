# Parse a Step Name into Year, Season, and Decimal Year

Converts a step name of the form `"step_YEAR_SEASON"` or `"YEAR_SEASON"`
(as returned by
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md)) into
a named list containing the clean year-season string, the integer year,
the integer season, and the decimal year used in
[`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md)
output.

## Usage

``` r
process_step(step)
```

## Arguments

- step:

  Character. Step name, e.g. `"step_5_3"` or `"5_3"`.

## Value

A named list with:

- `year_season`:

  Character. Clean `"year_season"` string (any `"step_"` prefix
  removed).

- `year`:

  Integer. The year component.

- `season`:

  Integer. The season component.

## See also

[`clean_steps`](https://danovando.github.io/marlin/reference/clean_steps.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md)

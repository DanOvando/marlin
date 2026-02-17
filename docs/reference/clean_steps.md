# Clean Step Names from simmar Output

Strips any `"step_"` prefix from step names in the output of
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
returning clean `"year_season"` format strings (e.g. `"5_3"` for year 5,
season 3).

## Usage

``` r
clean_steps(step)
```

## Arguments

- step:

  Character. Step name to clean, e.g. `"step_5_3"` or `"5_3"`.

## Value

Character string in `"year_season"` format.

## See also

[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md)

## Examples

``` r
clean_steps("step_1_2")  # Returns "1_2"
clean_steps("10_4")      # Returns "10_4" (no-op)
```

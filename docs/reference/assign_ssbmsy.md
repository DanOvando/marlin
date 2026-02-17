# Find Baseline SSBmsy

And assigns it to each critter in fauna object

## Usage

``` r
assign_ssbmsy(fauna, sel_start = 0.01, years = 50)
```

## Arguments

- fauna:

  a fauna object

- sel_start:

  selectivity as a multiplier of length at maturity

- years:

  the number of years to run things out for

- mult:

  multiplier to effort

- use:

  one of "graphs" or "optim"

## Value

results of simulation

## Examples

``` r
if (FALSE) { # \dontrun{
years <- 50

seasons <- 4

fauna <-
  list(
    "bigeye" = create_critter(
      scientific_name = "thunnus obesus",
      adult_diffusion = 10,
      density_dependence = "post_dispersal",
      seasons = seasons,
      resolution = resolution,
      age_mature = 1,
      steepness = 0.9,
      ssb0 = 1000
    )
  )

baseline_ssbmsy <- find_ssbmsy(fauna = fauna)
} # }
```

# Tune Diffusion

Tunes the diffusion parameters to achieve a desired "home range". The
home range is set by finding the diffusion rates that results in 95% or
less of animals being within "home range" distance of their starting
patch after one year.

## Usage

``` r
tune_diffusion(home_range, mode = "opt")
```

## Arguments

- home_range:

  the desired home range, defined as the linear distance such that 95%
  of animals travel that distance or less from their starting location

- mode:

  one of "opt" to optimize diffusion rate or "plot" to return diffusion
  matrix

## Value

The tuned diffusion rate

## Examples

``` r
home_range = 42
diffusion_rate <- tune_diffusion(home_range)
diffusion_rate
#> [1] 92.9285
```

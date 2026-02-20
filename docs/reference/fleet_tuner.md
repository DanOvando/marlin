# Objective function for depletion-based fleet tuning

Given a vector of log fishing mortalities (one per species), sets fleet
catchabilities accordingly (via the logistic link function to keep q in
(0, 1)), runs a forward simulation, and returns log-ratio residuals
between achieved and target depletion.

## Usage

``` r
fleet_tuner(log_fs, fauna, fleets, e_fl, years = 50)
```

## Arguments

- log_fs:

  numeric vector of log instantaneous fishing mortality, one element per
  species in `fauna`

- fauna:

  a fauna object

- fleets:

  a fleet object (already set to constant effort)

- e_fl:

  numeric vector of effective effort per fleet

- years:

  number of years to simulate

## Value

numeric vector of residuals (log achieved depletion - log target
depletion), one per species

## Details

Used as the objective function for
[`nleqslv`](https://rdrr.io/pkg/nleqslv/man/nleqslv.html) inside
[`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md).

# Generate Length-at-Age Key

Produces an age-by-length-bin matrix giving the probability of being in
each length bin at each age, based on the specified growth model.

## Usage

``` r
generate_length_at_age_key(
  min_age,
  max_age,
  length_bin_width = 1,
  growth_params = list(),
  growth_model = "von_bertalanffy",
  cv,
  time_step = 1,
  linf_buffer = 10
)
```

## Arguments

- min_age:

  Numeric. Minimum age tracked in the model. Best left at 0, as the
  model does not explicitly track recruitment delays.

- max_age:

  Numeric. Maximum age tracked by the model (individuals this age or
  older are in the plus group).

- length_bin_width:

  Numeric. Width of each length bin in the key (default 1).

- growth_params:

  Named list of growth parameters. Contents depend on `growth_model`:
  for `"von_bertalanffy"`: `linf`, `vbk`, `t0`; for `"power"`:
  `length_a`, `length_b`, `t0`; for `"growth_cessation"`: `l0`, `rmax`,
  `k`, `t50`.

- growth_model:

  Character. Growth model to use. One of `"von_bertalanffy"`, `"power"`,
  or `"growth_cessation"`.

- cv:

  Numeric. Coefficient of variation of length-at-age (log-space).

- time_step:

  Numeric. Time step the model is running on (1 / seasons).

- linf_buffer:

  Numeric. Multiplier around Linf used to set the upper bound of the
  length key, accounting for fish larger than Linf (default 10).

## Value

A tibble with columns `age`, `length_bin`, `mean_length_at_age`,
`sigma_at_age`, `next_length_bin`, and `p_bin` (probability of being in
each length bin at each age).

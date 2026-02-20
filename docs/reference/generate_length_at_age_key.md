# generate_length_to_age_key

produces an age by length bins matrix with probability of being in
length bin at age

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

  minimum age tracked in the model. Best to leave at 0, as the model
  does not explicitly track delays for recruitment

- max_age:

  maximum age tracked by the model (individuals this age or older are
  tracked in the plus group)

- cv:

  the coefficient of variation of the length-at-age relationship
  (log-space)

- time_step:

  the time step the model is running on (1 / seasons)

- linf_buffer:

  multiplier around linf to create length at age key, taking into
  account that some fish will be larger than Linf

- k:

  the

- linf:

  asymptotic length of the species in a von Bertalanffy growth function

- t0:

  hypothetical age at which the fish would have length 0 (e.g. -0.5
  years)

## Value

a length-at-age key

# Iteratively Optimise an MPA Network

Uses a stochastic importance-resampling (SIR) algorithm to identify an
MPA network that maximises a weighted combination of biodiversity and
economic objectives. At each iteration the function evaluates the
marginal value of toggling individual patches in or out of the MPA, then
expands or contracts the network by the patch(es) with the highest
marginal objective value.

## Usage

``` r
optimize_mpa(
  fauna,
  fleets,
  starting_conditions = NA,
  alpha = 0.33,
  max_prop_mpa = 1,
  resolution,
  prop_sampled = 0.2,
  max_delta = 2,
  workers = 6,
  bio_objective = "max_ssb",
  econ_objective = "yield",
  work_backwards = TRUE,
  patches_at_a_time = 1
)
```

## Arguments

- fauna:

  Named list of fauna objects from
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md).

- fleets:

  Named list of fleet objects from
  [`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md),
  tuned with
  [`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md).

- starting_conditions:

  Initial population state, typically `sim[[length(sim)]]` from a prior
  [`simmar`](https://danovando.github.io/marlin/reference/simmar.md)
  run. Passed as `initial_conditions` to each inner `simmar` call.

- alpha:

  Numeric in 0, 1. Weight on biodiversity (SSB/SSB0) in the composite
  objective. `alpha = 1` ignores economics; `alpha = 0` ignores biology.
  Default `0.33`.

- max_prop_mpa:

  Numeric in (0, 1\]. Maximum proportion of patches that can be
  designated as MPA. Setting this below 1 forces
  `work_backwards = FALSE` automatically. Default `1`.

- resolution:

  Integer scalar or length-2 vector `c(nx, ny)`. Must match the
  resolution of `fauna` and `fleets`.

- prop_sampled:

  Numeric in (0, 1\]. Proportion of candidate patches to evaluate at
  each iteration. Lower values speed up computation but increase
  stochasticity. Default `0.2`.

- max_delta:

  Numeric. Unused; kept for backwards compatibility.

- workers:

  Integer. Number of parallel workers for
  [`future::multisession`](https://future.futureverse.org/reference/multisession.html).
  Default `6`.

- bio_objective:

  Character. Biodiversity sub-objective:

  `"max_ssb"`

  :   (default) Sum of SSB/SSB0 across species; favours absolute
      biodiversity gains.

  `"min_loss"`

  :   Count of species that do not decline relative to baseline; favours
      preventing any losses.

- econ_objective:

  Character. Economic sub-objective used to evaluate each candidate
  patch:

  `"yield"`

  :   (default) Total catch across species.

  `"revenue"`

  :   Total revenue.

  `"profits"`

  :   Total profit.

- work_backwards:

  Logical. If `TRUE` (default and recommended), start at full MPA
  coverage and iteratively remove patches. Set to `FALSE` to grow the
  network from empty.

- patches_at_a_time:

  Integer. Number of patches to toggle per iteration. Default `1`.

## Value

A named list:

- `outcomes`:

  Tibble of biodiversity and economic outcomes at each MPA size level.

- `objective`:

  Tibble of the composite objective value (`obj`) at each MPA size, plus
  component scores (`ssb`, `econ`, `yield`).

- `mpa_network`:

  Tibble mapping each iteration (proportion protected) to the
  corresponding MPA location data frame.

## Details

### Algorithm

When `work_backwards = TRUE` (recommended), the algorithm starts with
100\\ improves the weighted objective. This "reverse greedy" approach
tends to produce more stable networks than the forward greedy (start
empty, add patches) because it avoids early lock-in to suboptimal
configurations.

The composite objective at each iteration is: \$\$obj = \alpha \cdot
\tilde{B} + (1 - \alpha) \cdot \tilde{E}\$\$ where \\\tilde{B}\\ and
\\\tilde{E}\\ are percent-rank-scaled biodiversity and economic
outcomes, respectively. Scaling by percent ranks means \\\alpha\\
strictly controls the weight given to biodiversity relative to economic
rank, not their absolute magnitudes.

Parallel evaluation of candidate patches uses
[`furrr::future_map`](https://furrr.futureverse.org/reference/future_map.html);
set `workers` to control the number of parallel sessions.

## See also

[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md),
[`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md)

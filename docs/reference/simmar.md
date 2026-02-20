# Run the marlin Spatially Explicit Fisheries Simulation

`simmar` is the main simulation engine for the marlin package. Given
initialised `fauna` and `fleets` objects, it advances the
age-structured, spatially explicit population forward in time for a
specified number of years, applying fishing mortality, natural
mortality, movement, and recruitment each time step.

## Usage

``` r
simmar(
  fauna = list(),
  fleets = list(),
  manager = list(),
  habitat = list(),
  years = 100,
  steps = NA,
  starting_season = NA,
  initial_conditions = NA,
  starting_step = NA,
  keep_starting_step = TRUE,
  log_rec_devs = NULL,
  cor_rec = diag(length(fauna))
)
```

## Arguments

- fauna:

  Named list of critter objects from
  [`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md).
  All critters must share the same `time_step` (i.e. the same number of
  `seasons`).

- fleets:

  Named list of fleet objects from
  [`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md),
  tuned with
  [`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md).

- manager:

  Named list of management actions. Supported elements:

  `mpas`

  :   List with `$locations` (a data frame with columns `x`, `y`, `mpa`;
      `mpa = TRUE` = closed) and `$mpa_year` (integer year when closures
      are activated).

  `quotas`

  :   Named numeric vector of annual catch caps per species (names match
      `fauna`). Effort is scaled down by an optimised multiplier when
      the cap would be exceeded.

  `effort_cap`

  :   Named numeric vector capping total effort per fleet in open-access
      / sole-owner models.

  `closed_seasons`

  :   Named list of integer vectors specifying which seasons are closed
      to fishing for each critter.

- habitat:

  Named list (one entry per critter) of time-varying habitat. Each entry
  is a list of `[ny, nx]` matrices (one per year or per time step). When
  provided, the movement matrix is updated each step based on the
  difference in habitat quality between adjacent patches.

- years:

  Numeric. Number of years to simulate. Ignored when `steps` is
  provided.

- steps:

  Integer. Alternative to `years`: number of time steps to simulate.
  Takes precedence over `years` when not `NA`.

- starting_season:

  Integer. Starting season within a year. Rarely needed; used when
  chaining partial-year runs.

- initial_conditions:

  Initial population state, typically `sim[[length(sim)]]` from a
  previous `simmar` call. Defaults to the unfished equilibrium embedded
  in each critter object.

- starting_step:

  Character. Step name of the form `"year_season"` to start from.
  Controls step labelling when chaining runs.

- keep_starting_step:

  Logical. If `TRUE` (default), the starting (initial-conditions) step
  is included in the output.

- log_rec_devs:

  Numeric matrix of pre-generated log recruitment deviates with rows =
  steps and columns = critters (names must match `fauna`). When `NULL`
  (default), deviates are generated internally using each critter's
  `sigma_rec` and `ac_rec`.

- cor_rec:

  Numeric \\n \times n\\ correlation matrix for recruitment deviates
  across species. Default `diag(length(fauna))` (independent
  recruitment).

## Value

A named list of length `years / time_step` (or `steps`), where each
element is named `"year_season"` and contains a named sub-list (one per
critter) with the population state at that step. Each critter's state
includes:

- `n_p_a`:

  Numbers by patch and age (matrix).

- `b_p_a`:

  Biomass by patch and age (matrix).

- `ssb_p_a`:

  Spawning stock biomass by patch and age (matrix).

- `c_p_a`:

  Catch in numbers by patch and age (matrix).

- `c_p_fl`:

  Catch by patch and fleet (matrix).

- `r_p_fl`:

  Revenue by patch and fleet (matrix).

- `prof_p_fl`:

  Profit by patch and fleet (matrix).

- `e_p_fl`:

  Effort by patch and fleet (data frame).

Pass output to
[`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md)
for tidying.

## Details

### Simulation loop (per time step)

1.  **Effort dynamics**: For each fleet, total effort for this step is
    determined by the fleet model (`"constant_effort"`, `"open_access"`,
    `"sole_owner"`, or `"manual"`). Open-access and sole-owner fleets
    use a tanh-based normalised profitability signal derived from the
    previous step's revenues and costs.

2.  **Spatial allocation**:
    [`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md)
    redistributes each fleet's total effort across patches using a
    multiplicative velocity-field update driven by the chosen
    `spatial_allocation` metric (`"rpue"`, `"profit"`, etc.) from the
    previous step's buffet.

3.  **Fishing & population dynamics**: For each critter, age-structured
    fishing mortality \\F\_{p,a}\\ is applied, followed by natural
    mortality, movement, and Beverton-Holt recruitment.

4.  **Yield accounting**:
    [`allocate_yields`](https://danovando.github.io/marlin/reference/allocate_yields.md)
    and
    [`aggregate_yields`](https://danovando.github.io/marlin/reference/aggregate_yields.md)
    build the per-fleet buffet of catch, revenue, profit, and CPUE by
    patch.

5.  **Marginal values (if needed)**: When any fleet uses
    `spatial_allocation = "marginal_profit"` or `"marginal_revenue"`, or
    `fleet_model = "sole_owner"`,
    [`calc_marginal_value`](https://danovando.github.io/marlin/reference/calc_marginal_value.md)
    appends patch-level marginal returns to the buffet.

6.  **Quota enforcement**: If `manager$quotas` is set for a species, a
    scalar effort multiplier is solved via `optim` to keep total catch
    at or below the quota.

### Step naming

Steps are named `"year_season"` (e.g. `"5_3"` = year 5, season 3). Use
[`clean_steps`](https://danovando.github.io/marlin/reference/clean_steps.md)
to strip any `"step_"` prefix.

### Initial conditions

By default the simulation starts from the unfished equilibrium embedded
in each critter object. Pass `initial_conditions = sim[[length(sim)]]`
from a prior run to chain simulations (e.g. a baseline followed by an
MPA scenario).

## See also

[`process_marlin`](https://danovando.github.io/marlin/reference/process_marlin.md),
[`plot_marlin`](https://danovando.github.io/marlin/reference/plot_marlin.md),
[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md),
[`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md),
[`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md),
[`go_fish`](https://danovando.github.io/marlin/reference/go_fish.md),
[`allocate_effort`](https://danovando.github.io/marlin/reference/allocate_effort.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic run
sim <- simmar(fauna = fauna, fleets = fleets, years = 50)

# With MPA starting in year 10
mpa_locs <- expand_grid(x = 1:10, y = 1:10) |>
  mutate(mpa = x <= 5)
sim_mpa <- simmar(
  fauna   = fauna,
  fleets  = fleets,
  years   = 50,
  manager = list(
    mpas = list(locations = mpa_locs, mpa_year = 10)
  )
)

# Chain: baseline + policy scenario
sim_base    <- simmar(fauna = fauna, fleets = fleets, years = 50)
sim_policy  <- simmar(fauna = fauna, fleets = fleets, years = 20,
                      initial_conditions = sim_base[[length(sim_base)]],
                      manager = list(quotas = list(tuna = 500)))
} # }
```

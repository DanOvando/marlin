# R6 Class Representing a Fish Population

The core population object in marlin. Stores all life-history
parameters, spatial structure, and state variables for an age-structured
fish population on a 2-D grid. Instances are created via
[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md);
users rarely need to call `Fish$new()` directly.

## Details

### Key public fields (after initialisation)

- `n_p_a_0`:

  Matrix `[patches, ages]` of numbers at unfished equilibrium. Used as
  `initial_conditions` in
  [`simmar`](https://danovando.github.io/marlin/reference/simmar.md).

- `ssb0`:

  Unfished spawning stock biomass (summed across patches and ages).
  Reference point for depletion calculations.

- `b0`:

  Unfished total biomass.

- `b0_p`:

  Vector of unfished biomass per patch.

- `movement_matrix`:

  List of sparse transition matrices (one per season block), each
  `[patches, patches]`, column-stochastic. Built from the
  diffusion-taxis process during initialisation.

- `length_at_age`:

  Numeric vector of mean length at each age class, computed from the
  growth model.

- `sel_at_age`:

  Default selectivity at age (used only as a reference; actual
  selectivity is set per-metier).

- `ages`:

  Numeric vector of age classes.

- `resolution`:

  Integer vector `c(nx, ny)`.

- `patches`:

  Total number of patches (`nx * ny`).

- `time_step`:

  Fraction of a year per time step (`1 / seasons`).

- `depletion`:

  Target B/B0 at the start of the simulation.

- `init_explt`:

  Initial exploitation rate.

- `sigma_rec`:

  Standard deviation of log recruitment deviates.

- `ac_rec`:

  First-order autocorrelation of log recruitment deviates.

- `grid`:

  Data frame with `x`, `y`, and `patch_area` columns defining the
  spatial grid.

### `swim()` method

The main per-step update method. Given a fishing mortality matrix
`f_p_a`, applies Baranov-style mortality, natural mortality, movement,
and Beverton-Holt recruitment, returning updated population arrays.

## See also

[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md),
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`Metier`](https://danovando.github.io/marlin/reference/Metier.md)

## Methods

### Public methods

- [`Fish$new()`](#method-Fish-new)

- [`Fish$plot()`](#method-Fish-plot)

- [`Fish$plot_movement()`](#method-Fish-plot_movement)

- [`Fish$swim()`](#method-Fish-swim)

- [`Fish$clone()`](#method-Fish-clone)

------------------------------------------------------------------------

### Method [`new()`](https://rdrr.io/r/methods/new.html)

Used to store parameters describing the behavior of an object of class
Fish

#### Usage

    Fish$new(
      common_name = NULL,
      scientific_name = NULL,
      linf = NA,
      vbk = NA,
      t0 = -0.5,
      cv_len = 0.1,
      l0 = 18.85,
      rmax = 37.24,
      k = 0.89,
      t50 = 4.57,
      length_a = 0.1,
      length_b = 3,
      length_units = "cm",
      length_bin_width = 1,
      growth_model = "von_bertalanffy",
      min_age = 0,
      max_age = NA,
      weight_a = NA,
      weight_b = NA,
      pups = 10,
      weight_units = "kg",
      fec_form = "weight",
      fec_expo = 1,
      fec_at_age = NULL,
      length_50_mature = NA,
      length_95_mature = NULL,
      delta_mature = 0.1,
      age_50_mature = NULL,
      age_95_mature = NULL,
      age_mature = NA,
      semelparous = FALSE,
      m = NA,
      m_at_age = NULL,
      steepness = 0.8,
      r0 = 10000,
      ssb0 = NA,
      b0 = NA,
      popsize_measure = "b0",
      density_dependence = "global_habitat",
      adult_diffusion = 4,
      recruit_diffusion = 10,
      adult_home_range = NULL,
      recruit_home_range = NULL,
      query_fishlife = TRUE,
      sigma_rec = 0,
      ac_rec = 0,
      cores = 4,
      mat_mode = "age",
      default_wb = 2.8,
      tune_weight = FALSE,
      linf_buffer = 1.2,
      resolution = NULL,
      patch_area = 1,
      habitat = list(),
      season_blocks = list(),
      recruit_habitat = NULL,
      depletion = 1,
      burn_years = 50,
      seasons = 1,
      explt_type = "f",
      init_explt = 0.1,
      get_common_name = FALSE,
      spawning_seasons = NULL,
      max_hab_mult = 2,
      lorenzen_m = TRUE,
      lorenzen_c = -1
    )

#### Arguments

- `common_name`:

  common name of the species (can be used to lookup information using
  `taxize` and `fishlife`)

- `scientific_name`:

  common name of the species (can be used to lookup information in
  `fishlife`)

- `linf`:

  asymptotic length of the species in a von Bertalanffy growth function

- `vbk`:

  growth parameter *k* of the species in a von Bertalanffy growth
  function

- `t0`:

  hypothetical age at which the fish would have length 0 (e.g. -0.5
  years)

- `cv_len`:

  coefficient of variation around the length at age relationship (in log
  space)

- `l0`:

  length at age 0 in growth cessation model

- `rmax`:

  related to maximum growth rate in growth cessation model

- `k`:

  steepness of the logistic function in the growth cessation model \>=0

- `t50`:

  age at the logistic function midpoint in the growth cessation model

- `length_a`:

  a coefficient in power function growth

- `length_b`:

  b power coefficient in power function growth

- `length_units`:

  units of the length at age function (arbitrary)

- `length_bin_width`:

  the width of the length bins in the length-at-age key, defaults to 1cm

- `growth_model`:

  one of "von_bertalanffy" or "power"

- `min_age`:

  minimum age tracked in the model. Best to leave at 0, as the model
  does not explicitly track delays for recruitment

- `max_age`:

  maximum age tracked by the model (individuals this age or older are
  tracked in the plus group)

- `weight_a`:

  alpha parameter in the allometric weight function alpha x length ^
  beta

- `weight_b`:

  beta parameter in the allometric weight function alpha x length ^ beta

- `pups`:

  number of pups per individual for animals with pups rather than larvae

- `weight_units`:

  units of the allometric weight function (defaults to kg)

- `fec_form`:

  one of of "weight" (default) or "pups". When "weight", fecundity is a
  function of weight. When "pups", constant number of pups per
  individual produced

- `fec_expo`:

  exponent for fecundity at weight relationship, 1 = isometric \> 1
  hyperallometric

- `fec_at_age`:

  manual vector of fecundity at age

- `length_50_mature`:

  length at 50% maturity in a logistic maturity ogive

- `length_95_mature`:

  length at 95% maturity in a logistic maturity ogive

- `delta_mature`:

  as an alternative, the different in units of length between
  length_50_mature and length_95_mature

- `age_50_mature`:

  age at 50% maturity in a logistic maturity ogive if maturity is age
  based

- `age_95_mature`:

  age at 95% maturity in a logistic maturity ogive if maturity is age
  based

- `age_mature`:

  an alternative option to just set one age mature, which ends up as the
  age_50_mature

- `semelparous`:

  TRUE or FALSE. When FALSE (default), individuals can reproduce
  multiple times. When TRUE, individuals can only spawn once, so
  mortality at increase increases as a function of maturity at age

- `m`:

  instantaneous natural mortality rate. When `lorenzen_m = TRUE`, this
  is the average natural mortality across all ages

- `m_at_age`:

  a vector of natural mortality age at age in case manually specified

- `steepness`:

  steepness parameter (h) in a Beverton-Holt spawner-recruit function

- `r0`:

  asymptotic number of recruits under unfished conditions

- `ssb0`:

  asymptotic spawning stock biomass of recruits under unfished
  conditions. Tunnes r0 to achieve

- `b0`:

  desired level of unfished biomass

- `popsize_measure`:

  whether unfished population size should be set based on biomass or
  spawning stock biomass

- `density_dependence`:

  timing and nature of density dependence in the Beverton-Holt spawner
  recruit function, one of one of
  'global_habitat','local_habitat','pre_dispersal','post_dispersal','global_ssb'.
  When 'global_habitat' density dependence is a function of the total
  spawning biomass across all patches, and recruits are then distributed
  proportional to recruit habitat quality. When 'local_habitat' density
  dependence is a function of the spawning biomass in each patch, and
  recruits are then distributed proportional to recruit habitat quality.
  When 'pre_dispersal', density dependence is a function of the spawning
  biomass in each patch, and recruits are then distributed from their
  home patch based on the `recruit_home_range`. When 'post_dispersal',
  larvae are distributed from their home patch based on
  `recruit_home_range`, and then density dependence happens based on
  spawning biomass in the settling patches. When 'global_ssb', density
  dependence is a function of the total spawning biomass and recruits
  are then distributed in space proportional to spawning biomass.

- `adult_diffusion`:

  diffusion parameter *D* in the CTMC movement function for "adults"
  (not recruits). Deprecated, best to use `adult_home_range`

- `recruit_diffusion`:

  diffusion parameter *D* in the CTMC movement function for recruits.
  Deprecated, best to use `recruit_home_range`

- `adult_home_range`:

  the desired home range of adults. Overrides adult_diffusion if set

- `recruit_home_range`:

  the desired home range of recruits Overrides recruit_diffusion if set

- `query_fishlife`:

  TRUE or FALSE to query `Fishlife` for missing life history values.
  When set to FALSE all required life history values must be supplied by
  the user

- `sigma_rec`:

  the standard deviation of recruitment deviates in log-normal space

- `ac_rec`:

  the autocorrelation of recruitment deviates

- `cores`:

  the number of cores used to tun the weight relationship if used
  (deprecated)

- `mat_mode`:

  specifies whether maturity is a function of age (default) or length

- `default_wb`:

  deprecated

- `tune_weight`:

  deprecated

- `linf_buffer`:

  multiplier around linf to create length at age key, taking into
  account that some fish will be larger than Linf

- `resolution`:

  a vector of length two with number of patches in X and Y dimensions

- `patch_area`:

  area of each patch (e.g. KM2). Either a scaler, a matrix, or a vector.
  If vector must be in correct patch order (arrange(x,y))

- `habitat`:

  a matrix with dimensions X and Y specifying quality of adult
  (non-recruit) habitat. Determines taxis matrix

- `season_blocks`:

  list with elements indicating blocks of seasons. For example, if there
  are four seasons, setting `season_blocks = list(c(1,2),c(3,4))`
  indicates that seasons 1 and 2 are one block, 3 and 4 another. Allows
  for the model to be run at fine time scales while allowing some
  processes like movement or spawning to operate at coarser scales

- `recruit_habitat`:

  a matrix with dimensions X and Y with quality of recruit habitat
  (scales r0 in space as well as recruit diffusion under applicable
  forms of density dependence)

- `depletion`:

  depletion (SSB/SSB0) under initial fished conditions

- `burn_years`:

  number of years used to burn in the population to tune parameters
  without analytical solutions like SSB0

- `seasons`:

  number of seasons per year. 4 would indicate quarterly time steps, 12
  monthly, 365 daily.

- `explt_type`:

  deprecated plot object of class fish

- `init_explt`:

  instantaneous fishing mortality rate under initial fished conditions

- `get_common_name`:

  TRUE or FALSE to lookup common name from scientific name. Requires
  internet connection

- `spawning_seasons`:

  seasons in which spawning occurs

- `max_hab_mult`:

  maximum value of that habitat matrix multiplier (to prevent some
  habitats from being \>\>\> good, defaults to 2)

- `lorenzen_m`:

  TRUE or FALSE to use the Lorenzen function to calculate natural
  mortality at age

- `lorenzen_c`:

  the rate of the Lorenzen curve. Defaults to -1, larger values will
  make the difference between natural mortality at young vs old ages
  less pronounced

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

#### Usage

    Fish$plot(type = 2)

#### Arguments

- `type`:

#### Returns

a plot of the life history values plot diffusion

------------------------------------------------------------------------

### Method `plot_movement()`

#### Usage

    Fish$plot_movement()

#### Arguments

- `type`:

#### Returns

a plot of the diffusion after one year Swim

Swim advances the population one time step

------------------------------------------------------------------------

### Method `swim()`

#### Usage

    Fish$swim(
      burn_steps = 0,
      season = 1,
      f_p_a = NULL,
      last_n_p_a = NULL,
      adult_movement = NULL,
      tune_unfished = 0,
      rec_devs = NA,
      move_fish = 1
    )

#### Arguments

- `burn_steps`:

  number of steps for burn in period if applicable

- `season`:

  the current season

- `f_p_a`:

  matrix of fishing mortality by patch and age

- `last_n_p_a`:

  matrix of initial numbers by patch and age

- `adult_movement`:

  the adult movement matrix

- `tune_unfished`:

  boolean indicating whether to tune unfished

- `rec_devs`:

  externally supplied recruitment deviates

- `move_fish`:

  1 means fish move, 0 means they do not. Only set to 0 when running
  exploratory fishing

#### Returns

the population in the next time step

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Fish$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

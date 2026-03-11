# Create a Critter (Species) Object

Constructs a critter (species population) object for use in
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md). Life
history parameters (growth, maturity, natural mortality, weight-at-age)
are looked up from FishLife when `scientific_name` is supplied, or can
be passed manually via `...`. The returned object is an R6 `Fish`
instance initialised to unfished equilibrium and burned in for
`burn_years` years.

## Usage

``` r
create_critter(
  common_name = NULL,
  scientific_name = NULL,
  get_common_name = FALSE,
  query_fishlife = TRUE,
  critter_type = "fish",
  habitat = list(),
  season_blocks = list(),
  recruit_habitat = NULL,
  seasons = 1,
  lorenzen_c = -1,
  fec_form = "weight",
  adult_home_range = 1,
  recruit_home_range = 2,
  adult_diffusion = NULL,
  recruit_diffusion = NULL,
  depletion = 0.4,
  fished_depletion = NULL,
  init_explt = 0.1,
  f = NULL,
  explt_type = "f",
  burn_years = 50,
  linf = NA,
  vbk = NA,
  t0 = -0.5,
  weight_a = NA,
  weight_b = NA,
  m = NA,
  max_age = NA,
  age_mature = NA,
  length_50_mature = NA,
  fec_expo = 1,
  resolution = c(10, 10),
  patch_area = 1,
  spawning_seasons = NULL,
  density_dependence = "global_habitat",
  steepness = 0.8,
  growth_model = "von_bertalanffy",
  ...
)
```

## Arguments

- common_name:

  Character or `NULL`. Common name of the species (e.g.
  `"yellowfin tuna"`). Used to look up life history from FishLife when
  `scientific_name` is not supplied.

- scientific_name:

  Character or `NULL`. Scientific name (e.g. `"Thunnus albacares"`).
  Preferred over `common_name` for FishLife lookup; case-insensitive.

- get_common_name:

  Logical. If `TRUE`, resolves the common name from `scientific_name`
  via an internet lookup. Default `FALSE`.

- query_fishlife:

  Logical. If `TRUE` (default), missing life history parameters are
  looked up from FishLife (requires internet for taxonomic resolution
  via taxize). If `FALSE`, no external lookups are performed; all
  required parameters must be supplied by the user. User-supplied values
  always take precedence over FishLife regardless of this setting.

- critter_type:

  Character. Currently only `"fish"` is supported. Placeholder for
  future non-fish implementations.

- habitat:

  List of habitat-quality matrices (one per season block). Each matrix
  must be `[ny, nx]`. Values represent relative habitat quality on a log
  scale; they drive movement taxis and, depending on
  `density_dependence`, recruit distribution. Generate with
  [`sim_habitat`](https://danovando.github.io/marlin/reference/sim_habitat.md).

- season_blocks:

  List of integer vectors specifying which seasons belong to each
  habitat block (e.g. `list(1:2, 3:4)` for a two-block, four-season
  year). Must align with `habitat`.

- recruit_habitat:

  Matrix or `NULL`. Habitat matrix used to distribute new recruits.
  Defaults to adult habitat of the first season block when `NULL`.

- seasons:

  Integer. Number of seasons per year. Sets `time_step = 1 / seasons`.
  All critters in a `fauna` list must share the same value.

- lorenzen_c:

  Numeric. Exponent in the Lorenzen (1996) size-dependent natural
  mortality function. Negative values (default `-1`) give higher
  mortality for smaller fish.

- fec_form:

  Character. Fecundity form: `"weight"` (fecundity proportional to body
  weight; default) or `"pups"` (fixed litter size, e.g. for
  elasmobranchs).

- adult_home_range:

  Numeric or `NULL`. Characteristic adult home range (km). When not
  `NULL`, this is used to derive the adult diffusion coefficient and
  takes precedence over `adult_diffusion`. Set to `NULL` to use
  `adult_diffusion` directly.

- recruit_home_range:

  Numeric or `NULL`. Characteristic recruit home range (km). When not
  `NULL`, this is used to derive the recruit diffusion coefficient and
  takes precedence over `recruit_diffusion`. Set to `NULL` to use
  `recruit_diffusion` directly.

- adult_diffusion:

  Numeric or `NULL`. Adult diffusion rate (km^2 per time step). Used
  only when `adult_home_range` is `NULL`.

- recruit_diffusion:

  Numeric or `NULL`. Recruit diffusion rate. Used only when
  `recruit_home_range` is `NULL`.

- depletion:

  Numeric in (0, 1\]. Target spawning stock biomass relative to unfished
  SSB (B/B0) at the start of the simulation.

- fished_depletion:

  Numeric. Alias for `depletion`; takes precedence when both are
  supplied.

- init_explt:

  Numeric. Initial annual exploitation rate (fraction of the exploitable
  population removed per year).

- f:

  Numeric. Alias for `init_explt`.

- explt_type:

  Character. One of `"f"` (instantaneous fishing mortality; default) or
  `"fmsy"` (as a multiple of Fmsy).

- burn_years:

  Integer. Number of years to burn in the unfished population to
  equilibrium before starting the simulation. Default 50.

- linf:

  Numeric or `NA`. Asymptotic length (L-infinity) in a von Bertalanffy
  growth model. Looked up from FishLife when `NA` and
  `query_fishlife = TRUE`. Required when
  `growth_model = "von_bertalanffy"`.

- vbk:

  Numeric or `NA`. Growth rate parameter (k) in a von Bertalanffy growth
  model. Looked up from FishLife when `NA`.

- t0:

  Numeric. Hypothetical age at length zero. Default `-0.5`.

- weight_a:

  Numeric or `NA`. Intercept alpha in the allometric weight-at-length
  relationship \\W = \alpha L^\beta\\. Looked up from FishLife when
  `NA`.

- weight_b:

  Numeric or `NA`. Exponent beta in the allometric weight-at-length
  relationship. Looked up from FishLife when `NA`.

- m:

  Numeric or `NA`. Instantaneous natural mortality rate. When
  `lorenzen_m = TRUE` (default in Fish), this is the asymptotic natural
  mortality. Looked up from FishLife when `NA`.

- max_age:

  Numeric or `NA`. Maximum age tracked by the model (plus group). Looked
  up from FishLife when `NA`; if still missing, computed as
  `-log(0.05) / m`.

- age_mature:

  Numeric or `NA`. Age at 50\\ for `age_50_mature`). Looked up from
  FishLife when `NA`. At least one maturity specification (`age_mature`,
  `age_50_mature + age_95_mature`, or `length_50_mature`) is required.

- length_50_mature:

  Numeric or `NA`. Length at 50\\ Looked up from FishLife when `NA`.
  Used when maturity is length-based.

- fec_expo:

  Numeric. Exponent for the fecundity-weight relationship. Values \> 1
  produce hyperallometric fecundity. Default `1`.

- resolution:

  Integer scalar or length-2 integer vector `c(nx, ny)` giving grid
  dimensions. A scalar is replicated to a square grid.

- patch_area:

  Numeric. Area of each patch (km^2). Used to scale diffusion and
  compute density-dependent quantities.

- spawning_seasons:

  Integer vector or `NULL`. Which seasons spawning occurs in. Defaults
  to all seasons when `NULL`.

- density_dependence:

  Character. Density dependence form; see Details. One of
  `"global_habitat"`, `"local_habitat"`, `"pre_dispersal"`,
  `"post_dispersal"`, `"global_ssb"`.

- steepness:

  Numeric in (0.2, 1). Beverton-Holt steepness parameter \\h\\ governing
  the shape of the spawner-recruit relationship.

- growth_model:

  Character. Growth model for length-at-age: `"von_bertalanffy"`
  (default), `"power"`, or `"growth_cessation"`.

- ...:

  Additional parameters forwarded to the
  [`Fish`](https://danovando.github.io/marlin/reference/Fish.md) R6
  class constructor. See
  [`?Fish`](https://danovando.github.io/marlin/reference/Fish.md) for
  the full list.

## Value

An R6 `Fish` object containing the age-structured population state at
unfished equilibrium, movement matrices, life-history schedules, and
reference points. Pass inside a named list as `fauna` to
[`simmar`](https://danovando.github.io/marlin/reference/simmar.md).

## Details

### Spatial setup

`resolution` determines the 2-D grid: `c(nx, ny)` produces `nx * ny`
patches indexed by `expand_grid(x, y)` (x varies slowest). `habitat`
should be a list (one matrix per season block) of `[ny, nx]`
habitat-quality matrices; these drive adult movement taxis and,
optionally, recruit distribution.

### Movement

Adult and recruit movement is modelled via diffusion-taxis. If
`adult_home_range` / `recruit_home_range` are provided (km), they are
used to derive diffusion coefficients automatically and take precedence
over `adult_diffusion` / `recruit_diffusion`. If a home range is `NULL`,
the corresponding diffusion value is used directly.

### Density dependence options

- `"global_habitat"`:

  Beverton-Holt spawner-recruit using global SSB; recruits distributed
  proportionally to habitat.

- `"local_habitat"`:

  Recruitment density dependence acts locally within each patch.

- `"pre_dispersal"`:

  Competition occurs before larval dispersal.

- `"post_dispersal"`:

  Competition occurs after larval dispersal.

- `"global_ssb"`:

  Beverton-Holt using global SSB with uniform spatial recruit
  distribution.

## See also

[`simmar`](https://danovando.github.io/marlin/reference/simmar.md),
[`create_fleet`](https://danovando.github.io/marlin/reference/create_fleet.md),
[`tune_fleets`](https://danovando.github.io/marlin/reference/tune_fleets.md),
[`sim_habitat`](https://danovando.github.io/marlin/reference/sim_habitat.md),
[`Fish`](https://danovando.github.io/marlin/reference/Fish.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal: look up life history from FishLife
tuna <- create_critter(
  scientific_name = "Thunnus albacares",
  resolution      = c(10, 10),
  seasons         = 4
)

# With explicit spatial structure
snapper <- create_critter(
  scientific_name    = "Lutjanus campechanus",
  adult_home_range   = 10,
  recruit_home_range = 20,
  resolution         = c(10, 10),
  seasons            = 1,
  density_dependence = "global_habitat",
  fished_depletion   = 0.4
)

fauna <- list(snapper = snapper)
} # }
```

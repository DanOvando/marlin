
<!-- README.md is generated from README.Rmd. Please edit that file -->

# marlin

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/mar)](https://CRAN.R-project.org/package=mar)
<!-- badges: end -->

marlin is a package for efficiently running simulations of marine fauna
and fisheries. Age-structured population model of different
(independent) animal types in a 2D system with multiple fishing fleets.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("DanOvando/marlin")
```

## Naviation

The core wrapper function is located in R/simmar.R. This funcion keeps
track of each of the populations and fleets.

The actual population models are found in src/fish\_model.cpp.
Additional modules will be put in there as they are developed

### Troubleshooting

Make sure you try the install with a fresh R session (go to
“Session\>Restart R” to make sure)

If you run into an error, first off try updating your R packages. From
there….

If your version of R is lower than 3.5, you might want to consider
updating R itself. Updating from 3.51 to 3.52 shouldn’t be any hassle.
BIG WARNING THOUGH, updating from say R 3.1 to 3.5 is a major update,
and you’ll lose all your installed packages in the process. I recommend
following the instructions
[here](https://www.datascienceriot.com/r/upgrade-R-packages/) to deal
with that, but even with that fix it can take a while, so I don’t
recommend doing a major R update if you’re on a deadline. There are also
packages to help you with this process, specifically
[`installR`](https://github.com/talgalili/installr/issues) for Windows
and [`updateR`](https://github.com/AndreaCirilloAC/updateR) for Mac.

From there…

  - On Windows, make sure you have the appropriate version of Rtools
    installed ([here](https://cran.r-project.org/bin/windows/Rtools/)),
    most likely Rtools35 if you have R version 3.3 or higher
      - Make sure that you select the box that says something about
        adding Rtools to the PATH variable
  - On macOS, there might be some issues with the your compiler,
    particularly if your version of R is less than 4.0.0.

If you get an error that says something like `clang: error: unsupported
option '-fopenmp'`, follow the instructions
[here](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos-before-r-4.0.0/)

Once you’ve tried those, restart your computer and try running

Below are a bunch of examples showing what `marlin` can do \#\# Simple
Example

Let’s start with a simple one-fleet one-critter example to illustrate
the various options in `marlin`

``` r
library(marlin)
library(tidyverse)
#> ── Attaching packages ────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──
#> ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
#> ✓ tibble  3.0.3     ✓ dplyr   1.0.1
#> ✓ tidyr   1.1.1     ✓ stringr 1.4.0
#> ✓ readr   1.3.1     ✓ forcats 0.5.0
#> ── Conflicts ───────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
options(dplyr.summarise.inform = FALSE)
theme_set(marlin::theme_marlin())

resolution <- 25 # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches 

years <- 20

seasons <- 1

time_step <- 1 / seasons

steps <- years * seasons

fauna <- 
  list(
    "bigeye" = create_critter(
      common_name = "bigeye tuna",
      adult_movement = 1,
      adult_movement_sigma = 10,
      rec_form = 1,
      seasons = seasons,
      fished_depletion = .5
    )
  )
#> ══  1 queries  ═══════════════
#> 
#> Retrieving data for taxon 'bigeye tuna'
#> ✔  Found:  bigeye+tuna[Common Name]
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0

# create a fleets object, which is a list of lists (of lists). Each fleet has one element, 
# with lists for each species inside there. Price specifies the price per unit weight of that 
# species for that fleet
# sel_form can be one of logistic or dome


fleets <- list(
  "longline" = create_fleet(
    list("bigeye" = Metier$new(
        critter = fauna$bigeye,
        price = 10,
        sel_form = "logistic",
        sel_start = 1,
        sel_delta = .01,
        catchability = 0,
        p_explt = 1
      )
    ),
    base_effort = resolution ^ 2
  )
)

a <- Sys.time()

fleets <- tune_fleets(fauna, fleets, tune_type = "depletion") 

Sys.time() - a
#> Time difference of 6.483664 secs


fauna$bigeye$plot()
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r

a <- Sys.time()

sim <- simmar(fauna = fauna,
                  fleets = fleets,
                  years = years)

Sys.time() - a
#> Time difference of 0.1297212 secs
```

we can then use `process_marlin` and `plot_marlin` to examine the
simulation

``` r

processed_marlin <- process_marlin(sim = sim, time_step = time_step)

plot_marlin(processed_marlin)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

``` r

plot_marlin(processed_marlin, plot_var = "c")
```

<img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" />

``` r

plot_marlin(processed_marlin, plot_var = "n", plot_type = "length", fauna = fauna)
#> Warning in plot_marlin(processed_marlin, plot_var = "n", plot_type = "length", :
#> trying to plot too many steps at once, cutting down to 10
#> dropping recruits from plot since drop_recruits = TRUE
```

<img src="man/figures/README-unnamed-chunk-3-3.png" width="100%" />

``` r

plot_marlin(processed_marlin, plot_var = "ssb", plot_type = "space")
#> Warning in plot_marlin(processed_marlin, plot_var = "ssb", plot_type = "space"):
#> Can only plot one time step for spatial plots, defaulting to last of the
#> supplied steps
```

<img src="man/figures/README-unnamed-chunk-3-4.png" width="100%" />

## Simple Example: seasonal habitat and movement and spatial q

Now let’s add in different adult habitats by season, different movement
rates by season, specified recruitment habitat, and a spatial dimension
to catchability, and quarterly time steps

``` r

seasons <- 4

time_step <-  1 / seasons

bigeye_habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  dnorm((x ^ 2 + y ^ 2), 300, 100)) %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()


bigeye_habitat2 <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  dnorm((x ^ .2 + y ^ .2), 100, 100)) %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()

bigeye_q <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  dplyr::mutate(habitat = rlnorm(resolution^2)) %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()


fauna <- 
  list(
    "bigeye" = create_critter(
      common_name = "bigeye tuna",
      seasonal_habitat = list(bigeye_habitat,bigeye_habitat2),
      season_blocks = list(c(1,2),c(3,4)),
     adult_movement = list(c(0,0),c(10,10)),# the mean number of patches moved by adults
      adult_movement_sigma = list(c(2,2), c(.1,.1)), # standard deviation of the number of patches moved by adults
      rec_form = 2,
      seasons = seasons,
      init_explt =  .1,
     explt_type = "f"
    )
  )
#> ══  1 queries  ═══════════════
#> 
#> Retrieving data for taxon 'bigeye tuna'
#> ✔  Found:  bigeye+tuna[Common Name]
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0

# create a fleets object, which is a list of lists (of lists). Each fleet has one element, 
# with lists for each species inside there. Price specifies the price per unit weight of that 
# species for that fleet
# sel_form can be one of logistic or dome


fleets <- list(
  "longline" = create_fleet(
    list("bigeye" = Metier$new(
        critter = fauna$bigeye,
        price = 10,
        sel_form = "logistic",
        sel_start = 1,
        sel_delta = .01,
        catchability = 1e-3,
        p_explt = 1,
        spatial_catchability = bigeye_q
      )
    ),
    base_effort = resolution ^ 2
  )
)

a <- Sys.time()

fleets <- tune_fleets(fauna, fleets) 

Sys.time() - a
#> Time difference of 3.820094 secs


fauna$bigeye$plot()
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r

a <- Sys.time()

sim2 <- simmar(fauna = fauna,
                  fleets = fleets,
                  years = years)

Sys.time() - a
#> Time difference of 1.374433 secs
  

processed_marlin <- process_marlin(sim = sim2, time_step = time_step)

plot_marlin(processed_marlin)
```

<img src="man/figures/README-unnamed-chunk-4-2.png" width="100%" />

``` r

plot_marlin(processed_marlin, plot_var = "c")
```

<img src="man/figures/README-unnamed-chunk-4-3.png" width="100%" />

``` r

plot_marlin(processed_marlin, plot_var = "n", plot_type = "length", fauna = fauna)
#> Warning in plot_marlin(processed_marlin, plot_var = "n", plot_type = "length", :
#> trying to plot too many steps at once, cutting down to 10
#> dropping recruits from plot since drop_recruits = TRUE
```

<img src="man/figures/README-unnamed-chunk-4-4.png" width="100%" />

``` r

plot_marlin(processed_marlin, plot_var = "ssb", plot_type = "space")
#> Warning in plot_marlin(processed_marlin, plot_var = "ssb", plot_type = "space"):
#> Can only plot one time step for spatial plots, defaulting to last of the
#> supplied steps
```

<img src="man/figures/README-unnamed-chunk-4-5.png" width="100%" />

## Two Species and two fleets with bells and whistles

``` r

seasons <- 4

steps <- years * seasons

time_step <- 1 / seasons
# for now make up some habitat


skipjack_habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  dplyr::mutate(habitat =  dnorm((x ^ 2 + y ^ 2), 20, 200)) %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()


bigeye_habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  dnorm((x ^ 2 + y ^ 2), 300, 100)) %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()


bigeye_habitat2 <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  dnorm((x ^ .2 + y ^ .2), 100, 100)) %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()

a <- Sys.time()

fauna <- 
  list(
    "skipjack" = create_critter(
      scientific_name = "Katsuwonus pelamis",
      seasonal_habitat = list(skipjack_habitat, skipjack_habitat), # pass habitat as lists
      season_blocks = list(c(1, 2), c(3, 4)), # seasons each habitat apply to
      recruit_habitat = skipjack_habitat,
      adult_movement = 2,# the mean number of patches moved by adults
      adult_movement_sigma = 2, # standard deviation of the number of patches moved by adults
      fished_depletion = .6, # desired equilibrium depletion with fishing (1 = unfished, 0 = extinct),
      rec_form = 1, # recruitment form, where 1 implies local recruitment
      seasons = seasons,
      init_explt = 0.2, 
      explt_type = "f"
      ),
    "bigeye" = create_critter(
      common_name = "bigeye tuna",
      seasonal_habitat = list(bigeye_habitat, bigeye_habitat2), # pass habitat as lists
      season_blocks = list(c(1, 2), c(3, 4)), # seasons each habitat apply to
      recruit_habitat = bigeye_habitat,
      adult_movement = 3,
      adult_movement_sigma = 1,
      fished_depletion = .1,
      rec_form = 1,
      seasons = seasons,
      init_explt = 0.3, 
      explt_type = "f"
    )
  )
#> ══  1 queries  ═══════════════
#> 
#> Retrieving data for taxon 'Katsuwonus pelamis'
#> ✔  Found:  Katsuwonus+pelamis
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0
#> ══  1 queries  ═══════════════
#> 
#> Retrieving data for taxon 'bigeye tuna'
#> ✔  Found:  bigeye+tuna[Common Name]
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0
Sys.time() - a
#> Time difference of 6.852328 secs

# create a fleets object, which is a list of lists (of lists). Each fleet has one element, 
# with lists for each species inside there. Price specifies the price per unit weight of that 
# species for that fleet
# sel_form can be one of logistic or dome


fleets <- list(
  "longline" = create_fleet(
    list(
       "skipjack" = Metier$new(
        critter = fauna$skipjack,
        price = 100,
        # price per unit weight
        sel_form = "logistic",
        # selectivity form, one of logistic or dome
        sel_start = .3,
        # percentage of length at maturity that selectivity starts
        sel_delta = .1,
        # additional percentage of sel_start where selectivity asymptotes
        catchability = .01,
        # overwritten by tune_fleet but can be set manually here
        p_explt = 1
      ),
       "bigeye" = Metier$new(
        critter = fauna$bigeye,
        price = 10,
        sel_form = "logistic",
        sel_start = .1,
        sel_delta = .01,
        catchability = 0,
        p_explt = 1
      )
    ),
    base_effort = resolution ^ 2
  ),
  "purseseine" = create_fleet(list(
    skipjack = Metier$new(
      critter = fauna$skipjack,
      price = 100,
      sel_form = "logistic",
      sel_start = 0.25,
      sel_delta = .1,
      catchability = .01,
      p_explt = 0.9
    ),
    bigeye = Metier$new(
      critter = fauna$bigeye,
      price = 100,
      sel_form = "logistic",
      sel_start = .25,
      sel_delta = .5,
      catchability = .01,
      p_explt = 1
    )
  ),
  base_effort = resolution ^ 2)
)

a <- Sys.time()

fleets <- tune_fleets(fauna, fleets) 

Sys.time() - a
#> Time difference of 6.629604 secs


# run simulations

# run the simulation using marlin::simmar
a <- Sys.time()

sim3 <- simmar(fauna = fauna,
                  fleets = fleets,
                  years = years)

Sys.time() - a
#> Time difference of 2.931982 secs
  
processed_marlin <- process_marlin(sim = sim3, time_step = time_step)

plot_marlin(processed_marlin)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r

plot_marlin(processed_marlin, plot_var = "c")
```

<img src="man/figures/README-example-2.png" width="100%" />

``` r

plot_marlin(processed_marlin, plot_var = "n", plot_type = "length", fauna = fauna)
#> Warning in plot_marlin(processed_marlin, plot_var = "n", plot_type = "length", :
#> trying to plot too many steps at once, cutting down to 10
#> dropping recruits from plot since drop_recruits = TRUE
```

<img src="man/figures/README-example-3.png" width="100%" />

``` r

plot_marlin(processed_marlin, plot_var = "ssb", plot_type = "space")
#> Warning in plot_marlin(processed_marlin, plot_var = "ssb", plot_type = "space"):
#> Can only plot one time step for spatial plots, defaulting to last of the
#> supplied steps
```

<img src="man/figures/README-example-4.png" width="100%" />

## Evaluating MPAs

Now let’s compare the effect of an MPA on two species

sharks live nearshore

``` r
library(marlin)
library(tidyverse)

theme_set(marlin::theme_marlin())

resolution <- 20 # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches 

seasons <- 1

years <- 50

tune_type <- "depletion"

steps <- years * seasons

# for now make up some habitat

yft_habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  .05 * x) %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()
 

mako_habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  dnorm(x,17,5)) %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()


# create a fauna object, which is a list of lists

fauna <- 
  list(
    "Yellowfin Tuna" = create_critter(
      scientific_name = "Thunnus albacares",
      seasonal_habitat = list(yft_habitat), # pass habitat as lists
      recruit_habitat = yft_habitat,
      adult_movement = 0,# the mean number of patches moved by adults
      adult_movement_sigma = 4, # standard deviation of the number of patches moved by adults
      fished_depletion = .4, # desired equilibrium depletion with fishing (1 = unfished, 0 = extinct),
      rec_form = 1, # recruitment form, where 1 implies local recruitment
      seasons = seasons,
      init_explt = 0.12, 
      explt_type = "f"
      ),
    "Shortfin Mako" = create_critter(
      scientific_name = "Isurus oxyrinchus",
      seasonal_habitat = list(mako_habitat), # pass habitat as lists
      recruit_habitat = mako_habitat,
      adult_movement = 5,
      adult_movement_sigma = 3,
      fished_depletion = .3,
      rec_form = 1,
      burn_years = 200,
      seasons = seasons,
      init_explt = .12, 
      explt_type = "f",
      fec_form = "linear",
      weight_a = 2 # average two offspring per shark
    )
  )
#> ══  1 queries  ═══════════════
#> 
#> Retrieving data for taxon 'Thunnus albacares'
#> ✔  Found:  Thunnus+albacares
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0
#> ══  1 queries  ═══════════════
#> 
#> Retrieving data for taxon 'Isurus oxyrinchus'
#> ✔  Found:  Isurus+oxyrinchus
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0

fauna$`Yellowfin Tuna`$plot()
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

``` r


# create a fleets object, which is a list of lists (of lists). Each fleet has one element, 
# with lists for each species inside there. Price specifies the price per unit weight of that 
# species for that fleet
# sel_form can be one of logistic or dome

fleets <- list("longline" = create_fleet(list(
  `Yellowfin Tuna` = Metier$new(
    critter = fauna$`Yellowfin Tuna`,
    price = 100, # price per unit weight
    sel_form = "logistic", # selectivity form, one of logistic or dome
    sel_start = .3, # percentage of length at maturity that selectivity starts
    sel_delta = .1, # additional percentage of sel_start where selectivity asymptotes
    catchability = .01, # overwritten by tune_fleet but can be set manually here
    p_explt = 1
    ),
  `Shortfin Mako` = Metier$new(
    critter = fauna$`Shortfin Mako`,
    price = 10,
    sel_form = "logistic",
    sel_start = .1,
    sel_delta = .01,
    catchability = 0,
    p_explt = 1
  )),
  mpa_response = "stay",
  base_effort = resolution^2
))

a <- Sys.time()

fleets <- tune_fleets(fauna, fleets, tune_type = tune_type) # tunes the catchability by fleet to achieve target depletion

Sys.time() - a
#> Time difference of 8.407858 secs

# run simulations

a <- Sys.time()

nearshore <- simmar(fauna = fauna,
                  fleets = fleets,
                  years = years)

Sys.time() - a
#> Time difference of 0.2205288 secs
  
proc_nearshore <- process_marlin(nearshore, time_step =  fauna[[1]]$time_step)
```

create an MPA

``` r
set.seed(42)
#specify some MPA locations
mpa_locations <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
mutate(mpa = x > 15 & y < 15)

mpa_locations %>% 
  ggplot(aes(x,y, fill = mpa)) + 
  geom_tile() + 
  scale_fill_brewer(palette = "Accent", direction  = -1, name = "MPA") + 
  scale_x_continuous(name = "Lat") + 
  scale_y_continuous(name = "Lon")
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

And now apply MPA

``` r
a <- Sys.time()

nearshore_mpa <- simmar(
  fauna = fauna,
  fleets = fleets,
  years = years,
  mpas = list(locations = mpa_locations,
              mpa_year = floor(years * .5))
)

Sys.time() - a
#> Time difference of 0.2250888 secs

proc_nearshore_mpa <- process_marlin(nearshore_mpa, time_step =  fauna[[1]]$time_step)
```

Now though, consider a different scenario. Here the tunas still slightly
prefer their same nearshore habitat, but now the shortfin mako
population primarily lives farther offshore. We will first simulate that
population without the MPA, and then assess the effects of the exact
same MPA on this new scenario.

``` r

yft_habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  .2 * x) %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()
 

mako_habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  dnorm(x, 9,5)) %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()


# create a fauna object, which is a list of lists

fauna <- 
  list(
    "Yellowfin Tuna" = create_critter(
      scientific_name = "Thunnus albacares",
      seasonal_habitat = list(yft_habitat), # pass habitat as lists
      recruit_habitat = yft_habitat,
      adult_movement = 0,# the mean number of patches moved by adults
      adult_movement_sigma = 20, # standard deviation of the number of patches moved by adults
      fished_depletion = .4, # desired equilibrium depletion with fishing (1 = unfished, 0 = extinct),
      rec_form = 1, # recruitment form, where 1 implies local recruitment
      seasons = seasons,
      init_explt = .12, 
      explt_type = "f"
      ),
    "Shortfin Mako" = create_critter(
      scientific_name = "Isurus oxyrinchus",
      seasonal_habitat = list(mako_habitat), # pass habitat as lists
      recruit_habitat = mako_habitat,
      adult_movement = 3,
      adult_movement_sigma = 1,
      fished_depletion = .3,
      rec_form = 1,
      burn_years = 200,
            seasons = seasons,
            init_explt = .12, 
      explt_type = "f",
          fec_form = "linear",
      weight_a = 2
    )
  )
#> ══  1 queries  ═══════════════
#> 
#> Retrieving data for taxon 'Thunnus albacares'
#> ✔  Found:  Thunnus+albacares
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0
#> ══  1 queries  ═══════════════
#> 
#> Retrieving data for taxon 'Isurus oxyrinchus'
#> ✔  Found:  Isurus+oxyrinchus
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0

fleets <- tune_fleets(fauna, fleets, tune_type = tune_type) # tunes the catchability by fleet to achieve target depletion

# run simulations

# run the simulation using marlin::simmar
a <- Sys.time()

offshore <- simmar(fauna = fauna,
                  fleets = fleets,
                  years = years)

Sys.time() - a
#> Time difference of 0.2119761 secs
  
proc_offshore <- process_marlin(offshore, time_step =  fauna[[1]]$time_step)

a <- Sys.time()

offshore_mpa_sim <- simmar(
  fauna = fauna,
  fleets = fleets,
  years = years,
  mpas = list(locations = mpa_locations,
              mpa_year = floor(years * .5))
)

Sys.time() - a
#> Time difference of 0.2097239 secs


proc_offshore_mpa <- process_marlin(offshore_mpa_sim, time_step =  fauna[[1]]$time_step)
```

``` r
plot_marlin(
  `MPA: Sharks Offshore` = proc_offshore_mpa,
  `No MPA` = proc_nearshore,
  `MPA: Sharks Nearshore` = proc_nearshore_mpa,
  steps_to_plot = NA,
  plot_var = "ssb"
)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

``` r
plot_marlin(
  `MPA: Sharks Offshore` = proc_offshore_mpa,
  `No MPA` = proc_nearshore,
  `MPA: Sharks Nearshore` = proc_nearshore_mpa,
  plot_var = "b",
  plot_type = "space",
  steps_to_plot = c(1,25,50))
#> Warning in plot_marlin(`MPA: Sharks Offshore` = proc_offshore_mpa, `No MPA` =
#> proc_nearshore, : Can only plot one time step for spatial plots, defaulting to
#> last of the supplied steps
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

## Defacto MPAs through bycatch penalties

Now consider a case where prices for shortfin mako are negative creating
defacto MPAs

``` r
years <- 100

tune_type <- "explt"

# make up some habitat

yft_habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  .05 * x) %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()


mako_habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  x > 12 & y >12) %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()


# create a fauna object, which is a list of lists

fauna <- 
  list(
    "Yellowfin Tuna" = create_critter(
      scientific_name = "Thunnus albacares",
      seasonal_habitat = list(yft_habitat), 
      recruit_habitat = yft_habitat,
      adult_movement = 0,
      adult_movement_sigma = 4, 
      fished_depletion = .4, 
      rec_form = 1, 
      seasons = seasons,
      init_explt = 0.12, 
      explt_type = "f"
    ),
    "Shortfin Mako" = create_critter(
      scientific_name = "Isurus oxyrinchus",
      seasonal_habitat = list(mako_habitat), 
      recruit_habitat = mako_habitat,
      adult_movement = 5,
      adult_movement_sigma = 3,
      fished_depletion = .3,
      rec_form = 1,
      burn_years = 200,
      seasons = seasons,
      init_explt = 0.1, 
      explt_type = "f"
    )
  )
#> ══  1 queries  ═══════════════
#> 
#> Retrieving data for taxon 'Thunnus albacares'
#> ✔  Found:  Thunnus+albacares
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0
#> ══  1 queries  ═══════════════
#> 
#> Retrieving data for taxon 'Isurus oxyrinchus'
#> ✔  Found:  Isurus+oxyrinchus
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0

# create a fleets object, accounting a negative price to shortfin makos

fleets <- list("longline" = create_fleet(list(
  `Yellowfin Tuna` = Metier$new(
    critter = fauna$`Yellowfin Tuna`,
    price = 100, # price per unit weight
    sel_form = "logistic", # selectivity form, one of logistic or dome
    sel_start = .3, # percentage of length at maturity that selectivity starts
    sel_delta = .1, # additional percentage of sel_start where selectivity asymptotes
    catchability = .01, # overwritten by tune_fleet but can be set manually here
    p_explt = 1
  ),
  `Shortfin Mako` = Metier$new(
    critter = fauna$`Shortfin Mako`,
    price = -20000,
    sel_form = "logistic",
    sel_start = .1,
    sel_delta = .01,
    catchability = 0,
    p_explt = 1
  )),
  mpa_response = "stay",
  base_effort = resolution^2
))

a <- Sys.time()

fleets <- tune_fleets(fauna, fleets, tune_type = tune_type) # tunes the catchability by fleet to achieve target depletion

Sys.time() - a
#> Time difference of 0.2139809 secs

# run simulations

# run the simulation using marlin::simmar
negative_prices <- simmar(fauna = fauna,
                          fleets = fleets,
                          years = years)

proc_negative_prices <- process_marlin(negative_prices, time_step =  fauna[[1]]$time_step)
```

``` r
plot_marlin(
  `De-Facto MPA` = proc_negative_prices,
  plot_var = "ssb",
  plot_type = "space",
  steps_to_plot = c(1,25,50))
#> Warning in plot_marlin(`De-Facto MPA` = proc_negative_prices, plot_var =
#> "ssb", : Can only plot one time step for spatial plots, defaulting to last of
#> the supplied steps
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

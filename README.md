
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

## Example

Create two critters, skipjack tuna and bigeye tuna, and simulate their
unfished conditions

``` r
library(marlin)
library(tidyverse)
#> ── Attaching packages ────────────────── tidyverse 1.3.0 ──
#> ✓ ggplot2 3.3.1     ✓ purrr   0.3.4
#> ✓ tibble  3.0.1     ✓ dplyr   1.0.0
#> ✓ tidyr   1.1.0     ✓ stringr 1.4.0
#> ✓ readr   1.3.1     ✓ forcats 0.5.0
#> ── Conflicts ───────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
options(dplyr.summarise.inform = FALSE)
resolution <- 25
habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  dnorm((x ^ 2 + y ^ 2), 20, 10))

habitat_mat <-
  matrix(
    rep(habitat$habitat, resolution),
    nrow = resolution ^ 2,
    ncol = resolution ^ 2,
    byrow = TRUE
  )

skj_hab <- habitat_mat / rowSums(habitat_mat)

habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  dnorm((x ^ 2 + y ^ 2), 600, 100))


# habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
#   mutate(habitat =  1)


habitat_mat <-
  matrix(
    rep(habitat$habitat, resolution),
    nrow = resolution ^ 2,
    ncol = resolution ^ 2,
    byrow = TRUE
  )

bet_hab <- habitat_mat / rowSums(habitat_mat)


fauna <-
  list(
    "skipjack" = create_critter(
      scientific_name = "Katsuwonus pelamis",
      habitat = skj_hab,
      adult_movement = 2,
      fished_depletion = .75
    ),
    "bigeye" = create_critter(
      common_name = "bigeye tuna",
      habitat = bet_hab,
      adult_movement = 10,
      fished_depletion = .3
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

# plot(fauna$bigeye$distance[10,],fauna$bigeye$move_mat[,10])

fauna$bigeye$move_mat %>% 
  as_tibble() %>% 
  mutate(x = 1:nrow(.)) %>% 
  pivot_longer(-x, names_to = "y", values_to = "movement") %>% 
  mutate(y = as.numeric(y)) %>% 
  ggplot(aes(x, y, fill = movement)) + 
  geom_tile()
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r


fleets <- list("longline" = list(
  skipjack = list(
    price = 0,
    sel_form = "logistic",
    sel_start = .9,
    sel_delta = .1,
    catchability = .1
  ),
  bigeye = list(
    price = 1000,
    sel_form = "logistic",
    sel_start = 1.25,
    sel_delta = .01,
    catchability = 0.1
  )
),
"purseseine" = list(
  skipjack = list(
    price = 100,
    sel_form = "logistic",
    sel_start = 0.25,
    sel_delta = .1,
    catchability = .9
  ),
  bigeye = list(
    price = 100,
    sel_form = "dome",
    sel_start = .5,
    sel_delta = .5,
    catchability = .1
  )
))



fleets <- create_fleet(fleets = fleets, fauna = fauna, base_effort = resolution^2)

fleets <- tune_fleets(fauna, fleets, steps = 100)

# run simulations
steps <- 100

a <- Sys.time()

storage <- simmar(fauna = fauna,
                  fleets = fleets,
                  steps = steps)

Sys.time() - a
#> Time difference of 0.8430099 secs
  

ssb_skj <- rowSums(storage[[steps]]$skipjack$ssb_p_a)

check <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(skj = ssb_skj)

ggplot(check, aes(x, y, fill = skj)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "skipjack")
```

<img src="man/figures/README-example-2.png" width="100%" />

``` r

ssb_bet <- rowSums(storage[[steps]]$bigeye$ssb_p_a)

check <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(bet = ssb_bet)

ggplot(check, aes(x, y, fill = bet)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "bigeye")
```

<img src="man/figures/README-example-3.png" width="100%" />

``` r


(sum(ssb_bet) / fauna$bigeye$ssb0) / fauna$bigeye$fished_depletion
#> [1] 1

(sum(ssb_skj) / fauna$skipjack$ssb0) / fauna$skipjack$fished_depletion
#> [1] 1
```

Now, add an MPA.

``` r

# mpa_locations <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
#   mutate(mpa =   (habitat$habitat < (qunif(.25, min = 0, max = max(habitat$habitat)))))


mpa_locations <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(mpa = between(x + y,10,15))

mpa_locations %>% 
  ggplot(aes(x,y, fill = mpa)) + 
  geom_tile()
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r


a <- Sys.time()

mpa_storage <- simmar(
  fauna = fauna,
  fleets = fleets,
  steps = steps,
  mpas = list(locations = mpa_locations,
              mpa_step = 50)
)

Sys.time() - a
#> Time difference of 0.5919321 secs

ssb_skj <- rowSums(mpa_storage[[steps]]$skipjack$ssb_p_a)

check <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(skj = ssb_skj)

ggplot(check, aes(x, y, fill = skj)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10") +
  labs(title = "skipjack")
```

<img src="man/figures/README-unnamed-chunk-2-2.png" width="100%" />

``` r

ssb_bet <- rowSums(mpa_storage[[steps]]$bigeye$ssb_p_a)

check <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(bet = ssb_bet)

# plot(check$bet[check$x == 1])

ggplot(check, aes(x, y, fill = bet)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10") +
  labs(title = "bigeye")
```

<img src="man/figures/README-unnamed-chunk-2-3.png" width="100%" />

``` r

# 
# (sum(ssb_bet) / fauna$bigeye$ssb0) / fauna$bigeye$fished_depletion
# 
# (sum(ssb_skj) / fauna$skipjack$ssb0) / fauna$skipjack$fished_depletion

bet_trajectory <- map_dbl(mpa_storage,~ sum(.x$bigeye$ssb_p_a))

plot(bet_trajectory)
```

<img src="man/figures/README-unnamed-chunk-2-4.png" width="100%" />

``` r

skj_trajectory <- map_dbl(mpa_storage,~ sum(.x$skipjack$ssb_p_a))

plot(skj_trajectory)
```

<img src="man/figures/README-unnamed-chunk-2-5.png" width="100%" />

Ah interesting, so need to think through the movement a bit more: the
problem is that movement is effectively 0 for the really good habitats:
critters stay put once they get there. Is that so bad? Problem is that
it doesn’t really allow for spillover.

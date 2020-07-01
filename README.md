
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
#> ── Attaching packages ───────────────────────────── tidyverse 1.3.0 ──
#> ✓ ggplot2 3.3.1     ✓ purrr   0.3.4
#> ✓ tibble  3.0.1     ✓ dplyr   1.0.0
#> ✓ tidyr   1.1.0     ✓ stringr 1.4.0
#> ✓ readr   1.3.1     ✓ forcats 0.5.0
#> ── Conflicts ──────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
options(dplyr.summarise.inform = FALSE)
  resolution <- 25
  habitat <- expand_grid(x =1:resolution, y = 1:resolution) %>%
    mutate(habitat =  dnorm((x^2 + y^2), 600,100))

  habitat_mat <-
    matrix(
      rep(habitat$habitat, resolution),
      nrow = resolution^2,
      ncol = resolution^2,
      byrow = TRUE
    )

  skj_hab <- habitat_mat / rowSums(habitat_mat)

  habitat <- expand_grid(x =1:resolution, y = 1:resolution) %>%
    mutate(habitat =  dnorm((x^2 + y^2), 100,2))

  habitat_mat <-
    matrix(
      rep(habitat$habitat, resolution),
      nrow = resolution^2,
      ncol = resolution^2,
      byrow = TRUE
    )

  bet_hab <- habitat_mat / rowSums(habitat_mat)


  fauna <-
    list("skipjack" = create_critter(
      scientific_name = "Katsuwonus pelamis",
      habitat = skj_hab,
      adult_movement = 2
    ),
    "bigeye" = create_critter(
      common_name = "bigeye tuna",
      habitat = bet_hab,
      adult_movement = 10
    ))
#> ══  1 queries  ═══════════════
#> 
#> Retrieving data for taxon 'Katsuwonus pelamis'
#> ✔  Found:  Katsuwonus+pelamis
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0


  # run simulations
  steps <- 100
  
  a <- Sys.time()
  
  storage <- marsim(fauna = fauna, 
         steps = steps)
  
  Sys.time() - a
#> Time difference of 0.609565 secs
  
  
  rec <- map(storage, ~.x[[2]]$n_p_a) %>%
    map_df(~tibble(rec = .x[,1]),.id = "i") %>%
    mutate(i = as.numeric(i)) %>%
    filter(i > 1) %>%
    group_by(i) %>%
    summarise(recs = sum(rec))

  ssb <- map(storage, ~.x[[2]]$ssb_p_a)%>%
    map_df(~tibble(ssb = rowSums(.x)),.id = "i") %>%
    mutate(i = as.numeric(i)) %>%
    filter(i > 1) %>%
    group_by(i) %>%
    summarise(ssb = sum(ssb))


  ssb_skj <- rowSums(storage[[steps]]$skipjack$ssb_p_a)

  check <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
    mutate(skj = ssb_skj)

   ggplot(check, aes(x,y, fill = skj)) +
    geom_tile() +
     scale_fill_viridis_c() + 
     labs(title = "skipjack")
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
   
  ssb_bet <- rowSums(storage[[steps]]$bigeye$ssb_p_a)

  check <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
    mutate(skj = ssb_bet)

   ggplot(check, aes(x,y, fill = ssb_bet)) +
    geom_tile() +
     scale_fill_viridis_c() + 
     labs(title = "bigeye")
```

<img src="man/figures/README-example-2.png" width="100%" />

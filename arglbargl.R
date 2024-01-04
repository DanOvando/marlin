library(marlin)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
theme_set(marlin::theme_marlin(base_size = 42))

resolution <- c(2,2) # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches

years <- 20

seasons <- 9

time_step <- 1 / seasons

steps <- years * seasons

fauna <-
  list(
    "bigeye" = create_critter(
      common_name = "bigeye tuna",
      adult_diffusion = 10,
      density_dependence = "post_dispersal",
      seasons = seasons,
      fished_depletion = .25,
      resolution = resolution,
      steepness = 0.6,
      ssb0 = 42,
      m = 0.4
    )
  )

fauna$bigeye$plot()

resolution <- fauna[[1]]$resolution

plot(fauna$bigeye$n_p_a_0[1,])

fleets <- list(
  "longline" = create_fleet(
    list("bigeye" = Metier$new(
      critter = fauna$bigeye,
      price = 10,
      sel_form = "logistic",
      sel_start = 1,
      sel_delta = .01,
      catchability = .1,
      p_explt = 1
    )
    ),
    base_effort = 0,
    resolution = resolution
  )
)

fleets$longline$metiers$bigeye$spatial_catchability

storage <- simmar(fauna = fauna,
                  fleets = fleets,
                  years = years)
a = process_marlin(storage)
plot_marlin(a)


# fleets <- tune_fleets(fauna, fleets, tune_type = "depletion")

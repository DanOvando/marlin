library(marlin)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
theme_set(marlin::theme_marlin(base_size = 42))

resolution <- 2 # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches 

years <- 50

seasons <- 1

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
      ssb0 = 1000
    )
  )

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
    base_effort = resolution ^ 2,
    resolution = resolution
  )
)

fauna$bigeye$plot()

fleets <- tune_fleets(fauna, fleets, tune_type = "depletion") 

start_time <- Sys.time()

for (i in 1:(10 * 80)){

example_sim <- simmar(fauna = fauna,
                      fleets = fleets,
                      years = years)
}
Sys.time() - start_time


library(marlin)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
theme_set(marlin::theme_marlin(base_size = 12))

resolution <- 28 # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches

years <- 50

seasons <- 1

time_step <- 1 / seasons

steps <- years * seasons

fauna <-
  list(
    "bigeye" = create_critter(
      common_name = "bigeye tuna",
      adult_home_range = 1,
      recruit_home_range = 4,
      density_dependence = "pre_dispersal",
      seasons = seasons,
      fished_depletion = .25,
      resolution = resolution,
      steepness = 0.6,
      ssb0 = 1000
    )
  )

fauna$bigeye$plot_movement()

fauna$bigeye$plot()

fishing_grounds <- expand.grid(x = 1:resolution, y = 1:resolution) |>
  mutate(fishing_ground = FALSE)

fishing_grounds$fishing_ground[1:2] <- TRUE

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
    )),
    base_effort = resolution^2,
    resolution = resolution,
    fishing_grounds = fishing_grounds
  )
)


fleets$longline$metiers$bigeye$plot_selectivity()


fleets <- tune_fleets(fauna, fleets, tune_type = "depletion")

fleets$longline$metiers$bigeye$plot_catchability()


start_time <- Sys.time()

for (i in 1:(1)) {
  example_sim <- simmar(
    fauna = fauna,
    fleets = fleets,
    years = years
  )
}
Sys.time() - start_time


proc_sim <- process_marlin(example_sim)

plot_marlin(proc_sim, plot_var = "c", plot_type = "space", max_scale = FALSE)

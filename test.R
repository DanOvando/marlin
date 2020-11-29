f <- .1

sel <- 1

library(marlin)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
theme_set(marlin::theme_marlin())

resolution <- 10 # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches 

years <- 20

seasons <- 1

time_step <- 1 / seasons

steps <- years * seasons


find_msy <- function(log_f, sel = 1) {

fauna <- 
  list(
    "bigeye" = create_critter(
      scientific_name = "Thunnus obesus",
      adult_movement = 1,
      adult_movement_sigma = 10,
      rec_form = 1,
      seasons = seasons,
      init_explt =  exp(log_f),
      resolution = resolution
    )
  )

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
      sel_start = sel,
      sel_delta = .01,
      catchability = .1,
      p_explt = 1 )
    ),
    base_effort = resolution ^ 2
  )
)

fleets <- tune_fleets(fauna, fleets, tune_type = "explt") 


sim <- simmar(fauna = fauna,
              fleets = fleets,
              years = years)


processed_marlin <-
  process_marlin(
    sim = sim,
    time_step = time_step,
    steps_to_keep = last(names(sim)),
    keep_age = FALSE
  )

# plot_marlin(processed_marlin, plot_var = "b")

# plot_marlin(processed_marlin, plot_var = "c")


yield <- -processed_marlin$fleets$catch %>% sum()
# yield
}

check <- tibble(log_f = log(seq(1e-6,2 , by = .05))) %>% 
  mutate(yields = map_dbl(log_f, find_msy))

check %>% 
  ggplot(aes(exp(log_f), yields)) + 
  geom_point()

library(marlin)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
theme_set(marlin::theme_marlin())

resolution <- 10 # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches 

years <- 200

seasons <- 1

time_step <- 1 / seasons

steps <- years * seasons

fauna <- 
  list(
    "bigeye" = create_critter(
      common_name = "bigeye tuna",
      adult_movement = 0,
      adult_movement_sigma = 4,
      rec_form = 3,
      seasons = seasons,
      fished_depletion = .2,
      resolution = resolution,
      steepness = 0.6,
      ssb0 = 1000
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
      sel_start = .25,
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


fauna$bigeye$plot()

a <- Sys.time()

sim <- simmar(fauna = fauna,
              fleets = fleets,
              years = years)

Sys.time() - a


max_delta <- .03

effort_cost_exponent <- fleets$longline$effort_cost_exponent

starting <- sim[[length(sim)]]

revenues <-  map_df(starting, ~.x$r_p_a_fl %>%
                      reshape2::melt() %>%
                      group_by(Var3) %>%
                      summarise(revenue = sum(value))) %>% 
  group_by(Var3) %>% 
  rename(fleet = Var3) %>% 
  summarise(revenue = sum(revenue))

effort <- map_df(starting[1], ~.x$e_p_fl) %>% 
  ungroup() %>% 
  mutate(patch = 1:nrow(.)) %>% 
  pivot_longer(-patch, names_to = "fleet", values_to = "effort") %>% 
  group_by(fleet) %>% 
  summarise(e2 = sum(((effort + 1)^effort_cost_exponent)),
            check = sum(effort))

# revenues$revenue <- 339613.9

profits <- revenues %>% 
  left_join(effort, by = "fleet") %>% 
  mutate(cost = revenue / e2) %>% 
  mutate(profit = revenue - cost * e2)

max_rev <- map_dbl(fauna, "ssb0")

prices = pluck(fleets, 1,1) %>% 
  map_dbl("price") 

max_p <- sum(max_rev * prices[names(max_rev)])

profits <- profits %>% 
  mutate(theta = log((1 + max_delta)) / (max_p))

og_fleet <- fleets

fleets$longline$cost_per_unit_effort <- profits$cost[profits$fleet == "longline"]

fleets$longline$profit_sensitivity <- profits$theta[profits$fleet == "longline"]

# fleets$longline$fleet_model <- "open access"

fleets$longline$spatial_allocation <- "revenue"


mpa_locations <- expand_grid(x = 1:resolution, y = 1:resolution) %>% 
  mutate(mpa = FALSE)

mpa_locations$mpa[1:75] <- TRUE

mpa_locations %>% 
  ggplot(aes(x,y, fill = mpa)) + 
  geom_tile() + 
  scale_fill_brewer(palette = "Accent", direction  = -1, name = "MPA") + 
  scale_x_continuous(name = "Lat") + 
  scale_y_continuous(name = "Lon")


a <- Sys.time()

sim <- simmar(fauna = fauna,
              fleets = fleets,
              years = years,
              mpas = list(locations = mpa_locations,
                          mpa_year = floor(years * .5)))

Sys.time() - a

processed_marlin <- process_marlin(sim = sim, time_step = time_step)


plot_marlin(processed_marlin, plot_var = "c", max_scale = FALSE)

plot_marlin(processed_marlin, plot_var = "ssb", max_scale = FALSE)




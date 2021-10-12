library(marlin)
library(tidyverse)

theme_set(marlin::theme_marlin())

resolution <- 20 # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches 

seasons <- 1

years <- 50

tune_type <- "depletion"

spatial_allocation <- "rpue"

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
      seasonal_habitat = yft_habitat, # pass habitat as lists
      recruit_habitat = yft_habitat,
      adult_movement = 0,# the mean number of patches moved by adults
      adult_movement_sigma = 10, # standard deviation of the number of patches moved by adults
      fished_depletion = .15, # desired equilibrium depletion with fishing (1 = unfished, 0 = extinct),
      rec_form = 1, # recruitment form, where 1 implies local recruitment
      seasons = seasons,
      init_explt = .3, 
      explt_type = "f"
    ),
    "Shortfin Mako" = create_critter(
      scientific_name = "Isurus oxyrinchus",
      seasonal_habitat = list(mako_habitat), # pass habitat as lists
      recruit_habitat = mako_habitat,
      adult_movement = 0,
      adult_movement_sigma = 10,
      fished_depletion = .3,
      rec_form = 1,
      burn_years = 200,
      seasons = seasons,
      init_explt = .2, 
      explt_type = "f",
      fec_form = "pups",
      pups = 2
    )
  )

fauna$`Shortfin Mako`$plot()
# create a fleets object, which is a list of lists (of lists). Each fleet has one element, 
# with lists for each species inside there. Price specifies the price per unit weight of that 
# species for that fleet
# sel_form can be one of logistic or dome

fleets <- list("longline" = create_fleet(list(
  `Yellowfin Tuna` = Metier$new(
    critter = fauna$`Yellowfin Tuna`,
    price = 100, # price per unit weight
    sel_form = "logistic", # selectivity form, one of logistic or dome
    sel_start = .1, # percentage of length at maturity that selectivity starts
    sel_delta = .1, # additional percentage of sel_start where selectivity asymptotes
    catchability = .01, # overwritten by tune_fleet but can be set manually here
    p_explt = 1
  ),
  `Shortfin Mako` = Metier$new(
    critter = fauna$`Shortfin Mako`,
    price = 0,
    sel_form = "logistic",
    sel_start = .1,
    sel_delta = .01,
    catchability = 0,
    p_explt = 1
  )),
  mpa_response = "stay",
  effort_cost_exponent = 1.3,
  base_effort = resolution^2,
  spatial_allocation = spatial_allocation
))

a <- Sys.time()

fleets <- tune_fleets(fauna, fleets, tune_type = tune_type) # tunes the catchability by fleet to achieve target depletion

Sys.time() - a

# run simulations

a <- Sys.time()

nearshore <- simmar(fauna = fauna,
                    fleets = fleets,
                    years = years)

Sys.time() - a


# assign costs


max_delta <- .03

effort_cost_exponent <- fleets$longline$effort_cost_exponent

starting <- nearshore[[length(nearshore)]]

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

fleets$longline$spatial_allocation <- "profit"


nearshore <- simmar(
  fauna = fauna,
  fleets = fleets,
  years = years)


proc_nearshore <- process_marlin(nearshore, time_step =  fauna[[1]]$time_step)

set.seed(42)
#specify some MPA locations
mpa_locations <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(mpa = x > 10 & y < 15)

mpa_locations %>% 
  ggplot(aes(x,y, fill = mpa)) + 
  geom_tile() + 
  scale_fill_brewer(palette = "Accent", direction  = -1, name = "MPA") + 
  scale_x_continuous(name = "Lat") + 
  scale_y_continuous(name = "Lon")

a <- Sys.time()

nearshore_mpa <- simmar(
  fauna = fauna,
  fleets = fleets,
  years = years,
  mpas = list(locations = mpa_locations,
              mpa_year = floor(years * .5))
)

Sys.time() - a

proc_nearshore_mpa <- process_marlin(nearshore_mpa, time_step =  fauna[[1]]$time_step)

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
      seasonal_habitat = yft_habitat, # pass habitat as lists
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
      fec_form = "pups",
      pups = 2
    )
  )

fauna$`Shortfin Mako`$plot()

fleets <- tune_fleets(fauna, fleets, tune_type = tune_type) # tunes the catchability by fleet to achieve target depletion

# run simulations

# run the simulation using marlin::simmar
a <- Sys.time()

offshore <- simmar(fauna = fauna,
                   fleets = fleets,
                   years = years)

Sys.time() - a

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


proc_offshore_mpa <- process_marlin(offshore_mpa_sim, time_step =  fauna[[1]]$time_step)


plot_marlin(
  `MPA: Sharks Offshore` = proc_offshore_mpa,
  `No MPA` = proc_nearshore,
  `MPA: Sharks Nearshore` = proc_nearshore_mpa,
  steps_to_plot = NA,
  plot_var = "c",
  max_scale = FALSE
)

plot_marlin(
  `MPA: Sharks Offshore` = proc_offshore_mpa,
  `No MPA` = proc_nearshore,
  `MPA: Sharks Nearshore` = proc_nearshore_mpa,
  steps_to_plot = NA,
  plot_var = "ssb",
  max_scale = FALSE
)


plot_marlin(
  `MPA: Sharks Offshore` = proc_offshore_mpa,
  `No MPA` = proc_nearshore,
  `MPA: Sharks Nearshore` = proc_nearshore_mpa,
  plot_var = "ssb",
  plot_type = "space",
  steps_to_plot = c(years-1))
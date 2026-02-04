

library(marlin)

library(tidyverse)

library(tictoc)

library(plot.matrix)

resolution <- c(20, 20) # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches

seasons <- 1

years <- 50

tune_type <- "depletion"

steps <- years * seasons

yft_home_range <- 6

yft_depletion <- 0.9

mako_depletion <- 0.4

mako_home_range <- 5

yft_b0 <- 1e9

mako_b0 <- 1e7

# for now make up some habitat

yft_habitat <- expand_grid(x = 1:resolution[1], y = 1:resolution[2]) %>%
  mutate(
    habitat = 0.9*x,
    habitat = habitat / max(habitat) * yft_home_range
  ) %>%
  pivot_wider(names_from = x, values_from = habitat) %>%
  select(-y) %>%
  as.matrix()

plot(yft_habitat)


mako_habitat <- expand_grid(x = 1:resolution[1], y = 1:resolution[2]) %>%
  mutate(
    habitat = dnorm(x, resolution, 100),
    habitat = habitat / max(habitat) * mako_home_range
  ) %>%
  pivot_wider(names_from = x, values_from = habitat) %>%
  select(-y) %>%
  as.matrix()


# create a fauna object, which is a list of lists
fauna <-
  list(
    "Yellowfin Tuna" = create_critter(
      scientific_name = "Thunnus albacares",
      habitat = yft_habitat, # pass habitat as lists
      recruit_habitat = yft_habitat,
      adult_home_range = yft_home_range, # cells per year
      depletion = yft_depletion, # desired equilibrium depletion with fishing (1 = unfished, 0 = extinct),
      density_dependence = "global_habitat", # recruitment form, where 1 implies local recruitment
      seasons = seasons,
      b0 = yft_b0
    ),
    "Shortfin Mako" = create_critter(
      scientific_name = "Isurus oxyrinchus",
      habitat = list(mako_habitat), # pass habitat as lists
      recruit_habitat = mako_habitat,
      adult_home_range = mako_home_range,
      depletion = mako_depletion,
      density_dependence = "global_habitat", # recruitment form, where 1 implies local recruitment
      burn_years = 200,
      seasons = seasons,
      fec_form = "pups",
      pups = 6,
      b0 = mako_b0,
      lorenzen_m = TRUE,
      lorenzen_c = -0.25
    )
  )

fauna$`Shortfin Mako`$plot()
# create a fleets object, which is a list of lists (of lists). Each fleet has one element,
# with lists for each species inside there. Price specifies the price per unit weight of that
# species for that fleet
# sel_form can be one of logistic or dome

fleets <- list("longline" = create_fleet(
  list(
    `Yellowfin Tuna` = Metier$new(
      critter = fauna$`Yellowfin Tuna`,
      price = 100, # price per unit weight
      sel_form = "logistic", # selectivity form, one of logistic or dome
      sel_start = .3, # percentage of lngth at maturity that selectivity starts
      sel_delta = .1, # additional percentage of sel_start where selectivity asymptotes
      catchability = .2, # overwritten by tune_fleet but can be set manually here
      p_explt = 1
    ),
    `Shortfin Mako` = Metier$new(
      critter = fauna$`Shortfin Mako`,
      price = 0,
      sel_form = "logistic",
      sel_start = .1,
      sel_delta = .01,
      catchability = 0.2,
      p_explt = 1
    )
  ),
  cost_per_unit_effort = 10000,
  effort_cost_exponent = 1.3,
  spatial_allocation = "marginal_profits",
  mpa_response = "stay",
  base_effort = prod(resolution),
  resolution = resolution
),
"longline2" = create_fleet(
  list(
    `Yellowfin Tuna` = Metier$new(
      critter = fauna$`Yellowfin Tuna`,
      price = 10, # price per unit weight
      sel_form = "logistic", # selectivity form, one of logistic or dome
      sel_start = .6, # percentage of length at maturity that selectivity starts
      sel_delta = .1, # additional percentage of sel_start where selectivity asymptotes
      catchability = .2, # overwritten by tune_fleet but can be set manually here
      p_explt = 1
    ),
    `Shortfin Mako` = Metier$new(
      critter = fauna$`Shortfin Mako`,
      price = 100,
      sel_form = "logistic",
      sel_start = .1,
      sel_delta = .01,
      catchability = 0.2,
      p_explt = 2
    )
  ),
  cost_per_unit_effort = 10000,
  effort_cost_exponent = 1.3,
  spatial_allocation = "marginal_profits",
  mpa_response = "stay",
  base_effort = 2*prod(resolution),
  resolution = resolution
))

a <- Sys.time()

# before <- fleets$longline$metiers$`Yellowfin Tuna`$catchability
#
# fleets$longline$metiers$`Yellowfin Tuna`$spatial_catchability

# tic()
# fleets1 <- tune_fleets(fauna, fleets, tune_type = "explt") # tunes the catchability by fleet to achieve target depletion
# toc()

# fleets <- fleets
#
# fleets$longline$spatial_allocation <- "ifdish"
# tic()
# fleets <- tune_fleets(fauna, fleets, tune_type = "depletion") # tunes the catchability by fleet to achieve target depletion
# toc()

tic()
fleets <- tune_fleets(fauna, fleets, tune_type = "depletion") # tunes the catchability by fleet to achieve target depletion
toc()
# fleets$longline$base_effort
# after =  fleets$longline$metiers$`Yellowfin Tuna`$catchability


# fleets$longline$metiers[[1]]$catchability
# fleets$longline$metiers[[1]]$catchability
# fleets$longline$metiers[[2]]$catchability

# new_fleets$longline$base_effort
# before == after

# fleets$longline$metiers$`Yellowfin Tuna`$spatial_catchability

# Sys.time() - a

# run simulations

# a <- Sys.time()

nearshore <- simmar(
  fauna = fauna,
  fleets = fleets,
  years = years
)

# Sys.time() - a

proc_nearshore <- process_marlin(nearshore, time_step = fauna[[1]]$time_step, keep_age = FALSE)

a = proc_nearshore$fleets  |>
  group_by(step, fleet) |>
  summarise(effort = sum(effort))

  proc_nearshore$fleets |>
  filter(step == min(step) | (step == max(step))) |>
  group_by(x,y,fleet,step) |>
  summarise(effort = unique(effort)) |>
  ggplot(aes(x,y,fill = effort)) +
  geom_tile() +
  facet_wrap(~step) +
  scale_fill_viridis_c()

plot_marlin(proc_nearshore, max_scale = FALSE, plot_var = "b")
plot_marlin(proc_nearshore, max_scale = FALSE, plot_var = "ssb")
plot_marlin(proc_nearshore, max_scale = FALSE, plot_var = "c")



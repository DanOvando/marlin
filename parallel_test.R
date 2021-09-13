
library(tidyverse)
library(marlin)
library(foreach)
library(doParallel)
library(furrr)
library(tictoc)

workers <- 4

seasons <- 4

n <- 250

years <- 10

resolution <- 20

steps <- years * seasons

time_step <- 1 / seasons

plan(strategy = multisession, workers = workers)

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
Sys.time() - a

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


# run simulations

# run the simulation using marlin::simmar
a <- Sys.time()

mpa_locations <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(mpa = x > 15 & y < 15)

task <- function(i,fauna, fleets, years, mpas){
  
  mpa_sim <- simmar(
    fauna = fauna,
    fleets = fleets,
    years = years,
    mpas = list(locations = mpas,
                mpa_year = 1)
  )
  

  biodiv <-
    (map_df(mpa_sim[[length(mpa_sim)]], ~ sum(.x$ssb_p_a) / .x$ssb0)) %>%
    pivot_longer(everything(), names_to = "critter",values_to = "biodiv")
  # calculate biodiversity component of objective function
  
  # econ <- sum(map_dbl(res, ~sum(.x$c_p_a))) #  calculate econ component of objective function
  #
  econ <-
    (map_df(mpa_sim[[length(mpa_sim)]], ~ sum(.x$r_p_a_fl, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "critter",values_to = "econ")
  #  calculate econ component of objective function, currently revenues across all fleets and species
  
  # out <- tibble(biodiv = biodiv, econ = econ)
  
  objective_outcomes <- biodiv %>%
    left_join(econ, by = "critter")
  
  outcomes <- list()
  
  outcomes$obj <- objective_outcomes
  
  outcomes$mpa <- mpas
  
  return(outcomes)
  
}


# test furrr --------------------------------------------------------------

tic()
furrr_test <- tibble(i = 1:n) %>%
  mutate(
    tmp = future_map(
      i,
      task,
      fauna = fauna,
      fleets = fleets,
      years = years,
      mpas = mpa_locations,
      .options = furrr_options(seed = TRUE),
      .progress = TRUE
    )
  )
toc()

# test dopar --------------------------------------------------------------

doParallel::registerDoParallel(cores = workers)

tic()
dopar_test <- foreach(i = 1:n) %dopar% {
  
  task(fauna = fauna, fleets = fleets, years = years, mpas = mpa_locations)
  
}
toc()

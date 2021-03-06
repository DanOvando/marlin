
## ----example-------------------------------------------------
library(marlin)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)


resolution <- 20 # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches 


# for now make up some habitat
skipjack_habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  dnorm((x ^ 2 + y ^ 2), 20, 200))

skipjack_habitat_mat <-
  matrix(
    rep(skipjack_habitat$habitat, resolution),
    nrow = resolution ^ 2,
    ncol = resolution ^ 2,
    byrow = TRUE
  )

skj_hab <- skipjack_habitat_mat / rowSums(skipjack_habitat_mat)

bigeye_habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  dnorm((x ^ 2 + y ^ 2), 300, 100))

# bigeye_habitat <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
#   mutate(habitat =  1)

bigeye_habitat_mat <-
  matrix(
    rep(bigeye_habitat$habitat, resolution),
    nrow = resolution ^ 2,
    ncol = resolution ^ 2,
    byrow = TRUE
  )

bet_hab <- bigeye_habitat_mat / rowSums(bigeye_habitat_mat)

# create a fauna object, which is a list of lists
# marlin::create_crutter will look up relvant life history information
# that you don't pass explicitly

fauna <-
  list(
    "skipjack" = create_critter(
      scientific_name = "Katsuwonus pelamis",
      habitat = skj_hab,
      adult_movement = 0,# the mean number of patches moved by adults
      adult_movement_sigma = 2, # standard deviation of the number of patches moved by adults
      fished_depletion = .3, # desired equilibrium depletion with fishing (1 = unfished, 0 = extinct)
      rec_form = 1, # recruitment form, where 1 implies local recruitment
      m = 0.5
    ),
    "bigeye" = create_critter(
      common_name = "bigeye tuna",
      habitat = bet_hab,
      adult_movement = 0,
      adult_movement_sigma = 2,
      fished_depletion = .4,
      m = 0.1,
      rec_form = 1
    )
  )

# plot(fauna$bigeye$distance[2,],fauna$bigeye$move_mat[,2])
# 
# fauna$bigeye$move_mat %>% 
#   as_tibble() %>% 
#   mutate(x = 1:nrow(.)) %>% 
#   pivot_longer(-x, names_to = "y", values_to = "movement") %>% 
#   mutate(y = as.numeric(y)) %>% 
#   ggplot(aes(x, y, fill = movement)) + 
#   geom_tile()


# create a fleets object, which is a list of lists (of lists). Each fleet has one element, 
# with lists for each species inside there. Price specifies the price per unit weight of that 
# species for that fleet
# sel_form can be one of logistic or dome


fleets <- list("longline" = list(
  skipjack = list(
    price = 100, # price per unit weight
    sel_form = "logistic", # selectivity form, one of logistic or dome
    sel_start = .99, # percentage of length at maturity that selectivity starts
    sel_delta = .01, # additional percentage of sel_start where selectivity asymptotes
    catchability = .1 # overwritten by tune_fleet but can be set manually here
  ),
  bigeye = list(
    price = 1000,
    sel_form = "logistic",
    sel_start = .99,
    sel_delta = .01,
    catchability = 0.1
  )
))



fleets <- create_fleet(fleets = fleets, fauna = fauna, base_effort = resolution^2) # creates fleet objects, basically adding in selectivity ogives

fleets <- tune_fleets(fauna, fleets, steps = 100) # tunes the catchability by fleet to achieve target depletion
## Note this will be a problem if there are more fleets than species, need to maybe assign proportion of catch that comes from 
## different fleets for each species?

# run simulations
steps <- 100 # number of time steps, in years for now


# run the simulation using marlin::simmar
a <- Sys.time()

storage <- simmar(fauna = fauna,
                  fleets = fleets,
                  steps = steps)

Sys.time() - a
  

# process results, will write some wrappers to automate this
ssb_skj <- rowSums(storage[[steps]]$skipjack$ssb_p_a)

check <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(skj = ssb_skj)

ggplot(check, aes(x, y, fill = skj)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "skipjack")

ssb_bet <- rowSums(storage[[steps]]$bigeye$ssb_p_a)

check <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(bet = ssb_bet)

ggplot(check, aes(x, y, fill = bet)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "bigeye")


# double check that target depletions are reached

(sum(ssb_bet) / fauna$bigeye$ssb0) / fauna$bigeye$fished_depletion

(sum(ssb_skj) / fauna$skipjack$ssb0) / fauna$skipjack$fished_depletion




## ------------------------------------------------------------


set.seed(42)

#specify some MPA locations
mpa_locations <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
 mutate(mpa = between(x,7,13) & between(y,7,13))

mpa_locations %>% 
  ggplot(aes(x,y, fill = mpa)) + 
  geom_tile()



# run the simulation, starting the MPAs at year 50 of the simulation
a <- Sys.time()

mpa_storage <- simmar(
  fauna = fauna,
  fleets = fleets,
  steps = steps,
  mpas = list(locations = mpa_locations,
              mpa_step = 50)
)

Sys.time() - a

ssb_skj <- rowSums(mpa_storage[[steps]]$skipjack$ssb_p_a)

check <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(skj = ssb_skj)

ggplot(check, aes(x, y, fill = skj)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "skipjack")

ssb_bet <- rowSums(mpa_storage[[steps]]$bigeye$ssb_p_a)

check <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(bet = ssb_bet)

# plot(check$bet[check$x == 1])

ggplot(check, aes(x, y, fill = bet)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "bigeye")

# 

mpa_storage[[77]]$bigeye$ssb_p_a -> a

plot(a[which(mpa_locations$mpa == FALSE)[10],])
# (sum(ssb_bet) / fauna$bigeye$ssb0) / fauna$bigeye$fished_depletion
# 
# (sum(ssb_skj) / fauna$skipjack$ssb0) / fauna$skipjack$fished_depletion

bet_outside_trajectory <- map_dbl(mpa_storage,~ sum(.x$bigeye$ssb_p_a[mpa_locations$mpa == FALSE,]))

plot(bet_outside_trajectory)


bet_inside_trajectory <- map_dbl(mpa_storage,~ sum(.x$bigeye$ssb_p_a[mpa_locations$mpa == TRUE,]))

plot(bet_inside_trajectory)

bet_outside_trajectory <- map_dbl(mpa_storage,~ sum(.x$bigeye$ssb_p_a[mpa_locations$mpa == FALSE,]))

plot(bet_outside_trajectory)


bet_trajectory <- map_dbl(mpa_storage,~ sum(.x$bigeye$ssb_p_a))

plot(bet_trajectory)

skj_trajectory <- map_dbl(mpa_storage,~ sum(.x$skipjack$ssb_p_a))

plot(skj_trajectory)



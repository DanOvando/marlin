library(marlin)
library(tidyverse)

years <- 20

resolution <-  10

seasons <- 4

steps <- years * seasons

time_step <-  1 / seasons

h1 <-  expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  .5 * x)

h2 <-  expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  -.5 * x)


h1 %>% 
  ggplot(aes(x,y,fill = habitat)) + 
  geom_tile()


bigeye_habitat <- h1 %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()


huh <- tidyr::pivot_longer(as.data.frame(bigeye_habitat), tidyr::everything())

h1$habitat2 <- as.numeric(huh$value)

h1 %>% 
  ggplot(aes(x,y,fill = habitat2)) + 
  geom_tile()


bigeye_habitat2 <- h2 %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()

bigeye_q <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
  dplyr::mutate(habitat = rlnorm(resolution^2)) %>% 
  pivot_wider(names_from = y, values_from = habitat) %>% 
  select(-x) %>% 
  as.matrix()


fauna <- 
  list(
    "bigeye" = create_critter(
      scientific_name = "thunnus obesus",
      base_habitat = list(bigeye_habitat,bigeye_habitat2),
      recruit_habitat = bigeye_habitat,
      season_blocks = list(c(1,2),c(3,4)),
      adult_diffusion = list(c(5,5), c(5,5)), # standard deviation of the number of patches moved by adults
      recruit_diffusion = 10,
      rec_form = 3,
      seasons = seasons,
      init_explt =  1,
      explt_type = "f"
    )
  )


fauna <- 
  list(
    "bigeye" = create_critter(
      scientific_name = "thunnus obesus",
      recruit_diffusion = 10,
      rec_form = 3,
      seasons = seasons,
      init_explt =  1,
      explt_type = "f"
    )
  )


fauna$bigeye$plot()
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
      sel_start = 1,
      sel_delta = .01,
      catchability = 1e-3,
      p_explt = 1,
      spatial_catchability = NA
    )
    ),
    base_effort = resolution ^ 2
  )
)

a <- Sys.time()

fleets <- tune_fleets(fauna, fleets) 

Sys.time() - a


a <- Sys.time()

sim2 <- simmar(fauna = fauna,
               fleets = fleets,
               years = years)

Sys.time() - a


processed_marlin <- process_marlin(sim = sim2, time_step = time_step)

processed_marlin$fauna %>% 
  filter(age == max(age)) %>% 
  group_by(step) %>% 
  mutate(n = n / sum(n)) %>% 
  ungroup() %>% 
  ggplot(aes(x,y,fill = n)) + 
  geom_tile() + 
  facet_wrap(~step) + 
  scale_fill_viridis_c()

plot_marlin(processed_marlin)

plot_marlin(processed_marlin, plot_var = "c")

plot_marlin(processed_marlin, plot_var = "n", plot_type = "length", fauna = fauna)

plot_marlin(processed_marlin, plot_var = "ssb", plot_type = "space", max_scale = TRUE)
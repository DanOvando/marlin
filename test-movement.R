library(marlin)
library(tidyverse)
library(gganimate)
years <- 20

resolution <-  10

seasons <- 4

steps <- years * seasons

time_step <-  1 / seasons

h1 <-  expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  dnorm(x, resolution / 2, .05 * resolution) *  dnorm(y, resolution / 2, .05 * resolution)) %>% 
  mutate(habitat = habitat * (x >= 4))

  

h2 <-  expand_grid(x = 1:resolution, y = 1:resolution) %>%
  mutate(habitat =  -.5 * x + 10) %>% 
  mutate(habitat = habitat * (x < 4))
  


h2 %>% 
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
      season_blocks = list(c(1,2),c(3,4)),
      adult_diffusion = list(c(5,5), c(5,5)), # standard deviation of the number of patches moved by adults
      recruit_diffusion = 10,
      rec_form = 0,
      seasons = seasons,
      init_explt =  1,
      explt_type = "f"
    )
  )


# fauna <- 
#   list(
#     "bigeye" = create_critter(
#       scientific_name = "thunnus obesus",
#       recruit_diffusion = 10,
#       rec_form = 3,
#       seasons = seasons,
#       init_explt =  1,
#       explt_type = "f"
#     )
#   )


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
    base_effort = resolution ^ 2,
    resolution = resolution
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



spawning_agg <- processed_marlin$fauna %>% 
  filter(age == max(age)) %>% 
  group_by(step) %>% 
  mutate(n = n / sum(n)) %>% 
  ungroup() %>% 
  ggplot(aes(x,y,fill = n)) + 
  geom_tile() + 
  transition_time(step) +
  ease_aes('linear') +
  scale_fill_viridis_c(name = "Tunas") + 
  scale_x_continuous(name = "longitude") + 
  scale_y_continuous(name = "latitude")


animate(spawning_agg, nframes = 100, fps=2)
anim_save("spawning_agg.gif")

plot_marlin(processed_marlin)

plot_marlin(processed_marlin, plot_var = "c")

plot_marlin(processed_marlin, plot_var = "n", plot_type = "length", fauna = fauna)

plot_marlin(processed_marlin, plot_var = "ssb", plot_type = "space", max_scale = TRUE)


# test habitat vector -----------------------------------------------------

range_shift <- vector(mode = "list", length = years)

for (i in 1:years){
  
  range_shift[[i]] <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
    mutate(habitat =  -((y - (1 + i))^2) / 10) %>% 
    pivot_wider(names_from = y, values_from = habitat) %>% 
    select(-x) %>% 
    as.matrix()
  
  
}



critter_habitat <- list(bigeye = range_shift)


image(critter_habitat$bigeye[[1]])

image(critter_habitat$bigeye[[years]])


a <- Sys.time()

sim_climate <- simmar(fauna = fauna,
               fleets = fleets,
               habitat = critter_habitat,
               years = years)

Sys.time() - a


processed_marlin <- process_marlin(sim = sim_climate, time_step = time_step, keep_age = FALSE)


range_shift <- processed_marlin$fauna %>% 
  group_by(step) %>% 
  mutate(n = n / max(n)) %>% 
  ungroup() %>% 
  ggplot(aes(x,y,fill = n)) + 
  geom_tile() + 
  transition_time(step) +
  ease_aes('linear') +
  scale_fill_viridis_c(name = "Tunas") + 
  scale_x_continuous(name = "longitude") + 
  scale_y_continuous(name = "latitude")


animate(range_shift, nframes = 100, fps=2)
anim_save("range_shift.gif")

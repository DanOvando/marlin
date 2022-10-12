library(marlin)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
theme_set(marlin::theme_marlin(base_size = 12))

resolution <- 20 # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches 

years <- 2 

seasons <- 1

time_step <- 1 / seasons

steps <- years * seasons


get_diffusion_frontier <- function(adult_diffusion, resolution = 10, years = 1, patch = 1){
  fauna <- 
    list(
      "bigeye" = create_critter(
        scientific_name = "thunnus obesus",
        adult_diffusion = (adult_diffusion),
        density_dependence = "post_dispersal",
        seasons = 1,
        resolution = resolution,
        steepness = 0.6,
        ssb0 = 1000,
        tune_diffusion = FALSE
      )
    )
  
  critter <- fauna[[1]]
  
  distance <-
    tidyr::expand_grid(x = 1:resolution, y = 1:resolution) %>%
    dist() %>%
    as.matrix()
  
  base_movement <-
    purrr::map2(critter$seasonal_diffusion,
                critter$base_habitat,
                ~ as.matrix(Matrix::expm((.x + .y) / seasons)))
  
  
  dist_org <- data.frame(distance = distance[patch,], p_move = base_movement[[1]][patch,] ) %>% 
    dplyr::arrange(distance)
  
  dist_org$cumulative_movement = cumsum(dist_org$p_move)
  
  move_frontier <- dist_org$distance[which.min((dist_org$cumulative_movement - 0.95)^2)[1]]

}


diffusion_frontier <- expand_grid(adult_diffusion = seq(0, 20, by = 0.5), resolution = seq(2,20, by = .2)) %>% 
  filter(adult_diffusion <= resolution) %>% 
  mutate(diffusion_frontier = map2_dbl(adult_diffusion, resolution, get_diffusion_frontier, years = 1)) 

diffusion_frontier %>% 
  ggplot(aes(diffusion_frontier, adult_diffusion, color = resolution, group = resolution)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_line(size = 2)


diffusion_frontier_model <- ranger(adult_diffusion ~ diffusion_frontier + resolution, data = diffusion_frontier %>% select(resolution,diffusion_frontier, adult_diffusion))

usethis::use_data(diffusion_frontier_model, overwrite = TRUE)

diffusion_frontier$pred <- predict(diffusion_frontier_model, data = diffusion_frontier)$predictions

diffusion_frontier %>% 
  ggplot(aes(diffusion_frontier, adult_diffusion, color = resolution, group = resolution)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_line(aes(diffusion_frontier, pred, color = resolution, group = resolution), alpha = 0.5) + 
  geom_point(size = 2)+
  scale_x_continuous(name = "True Max Cells Moved") + 
  scale_y_continuous(name = "Diffusion Parameter") + 
  scale_color_viridis_c(name = "Spatial Resolution") + 
  labs(caption = "Points are simulated values, lines are interopolation model")


diffusion_frontier %>% 
  filter(resolution == 10) %>% 
  ggplot(aes(diffusion_frontier, adult_diffusion, color = resolution, group = resolution)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point(size = 2)+
  geom_line(aes(diffusion_frontier, pred, color = resolution, group = resolution))

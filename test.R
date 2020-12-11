library(marlin)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
theme_set(marlin::theme_marlin())

resolution <- 10 # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches 

years <- 50

seasons <- 1

time_step <- 1 / seasons

steps <- years * seasons

# fauna <-
#   list(
#     "bigeye" = create_critter(
#       scientific_name = "gadus morhua",
#       age_95_mature = 3,
#       adult_movement = 0,
#       adult_movement_sigma = 50,
#       recruit_movement = 0,
#       recruit_movement_sigma = 50,
#       rec_form = 2,
#       seasons = seasons,
#       fished_depletion = .25,
#       resolution = resolution
#     )
#   )
# 
# fauna$bigeye$plot()


fauna <-
  list(
    "bigeye" = create_critter(
      common_name = "bigeye tuna",
      age_95_mature = 3,
      adult_movement = 0,
      adult_movement_sigma = 50,
      recruit_movement = 0,
      recruit_movement_sigma = 50,
      rec_form = 2,
      seasons = seasons,
      fished_depletion = .25,
      resolution = resolution
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
      catchability = 0,
      p_explt = 1
    )
    ),
    mpa_response = "stay",
    base_effort = resolution ^ 2
  )
)

a <- Sys.time()

test <- tune_fleets(fauna, fleets, tune_type = "explot") 

Sys.time() - a



a <- Sys.time()

sim <- simmar(fauna = fauna,
              fleets = fleets,
              years = years)

Sys.time() - a

fleets <- tune_fleets(fauna, fleets, tune_type = "explt") 

effort_mult = 1.5


find_msy <- function(effort_mult, fauna, fleets, opt = TRUE, target_critter) {
  
  
  tmp_fleets <- fleets
  
  for (f in seq_along(fleets)){
    
    tmp_fleets[[f]]$base_effort <- tmp_fleets[[f]]$base_effort * effort_mult
    
  }
  
  
  
  sim <- simmar(fauna = fauna,
                fleets = tmp_fleets,
                years = years)
  
  
  processed_marlin <-
    process_marlin(
      sim = sim,
      time_step = time_step,
      steps_to_keep = last(names(sim)),
      keep_age = FALSE
    )
  
  fauna_yield <- processed_marlin$fleets %>% 
    dplyr::filter(critter == target_critter) %>% 
    dplyr::summarise(yield = sum(catch))
  
  out <- -sum(fauna_yield$yield)
  
  if (opt == FALSE){
    
    out <- processed_marlin
  }
  # yield
  # 
  return(out)
}

assign_ref_points <- function(fauna, fleets){
  
  for (f in seq_along(fauna)){
    
    msy_mult <- nlminb(1e-3,find_msy, lower = 0, upper = 10, fauna = fauna, fleets = fleets, opt = TRUE,target_critter = names(fauna)[f])
    
    msy_state <- find_msy(msy_mult$par, fauna = fauna, fleets = fleets, opt = FALSE, target_critter =  names(fauna)[f])
    
    ref_points <- msy_state$fauna %>% 
      dplyr::filter(critter == names(fauna)[f]) %>% 
      dplyr::summarise(ssb_msy = sum(ssb), 
                       b_msy  = sum(b),
                       n_msy = sum(n),
                       msy = sum(c),
                       u_msy = sum(c) / sum(b))
    
    base_e_msy <- mean(purrr::map_dbl(fleets,"base_effort")) * msy_mult$par
    
    ref_points$base_e_msy_mult <- msy_mult$par
    
    ref_points$base_e_msy <- base_e_msy
    
    fauna[[f]]$ref_points <- ref_points
    
  }
  
  return(fauna)
}

fauna <- assign_ref_points(fauna = fauna, fleets = fleets)

# now, loop over MPA size and and fishing mortality rate and movement rates

mpa_test_grid <- expand_grid(mpa_size = seq(0,1, by = 0.05), f_fmsy = seq(0,3, by = 0.5))


mpa_test <- function(mpa_size, f_fmsy, fauna,fleets, random_mpas = FALSE){
  
  # f_fmsy <- 2
  # 
  # mpa_size <- 0.2
  # 
  base_e_msy <- purrr::map_dbl(fauna,c("ref_points","base_e_msy")) 
  
  msy <- purrr::map_dbl(fauna,c("ref_points","msy")) 
  
  
  tmp_fleets <- fleets
  
  tmp_fleets <- purrr::modify_in(tmp_fleets,list(1, "base_effort"), ~ base_e_msy * f_fmsy)
  
  if (random_mpas == TRUE) {
    mpa_locations <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
      mutate(patch = 1:nrow(.)) %>%
      mutate(mpa = patch %in% sample(1:nrow(.), round(mpa_size * nrow(.))))
  } else {
    mpa_locations <-
      expand_grid(x = 1:resolution, y = 1:resolution) %>%
      mutate(patch = 1:nrow(.)) %>%
      mutate(mpa = patch %in% (1:(round(nrow(
        .
      ) * mpa_size))))
  }
  
  # mpa_locations %>% 
  #   ggplot(aes(x,y, fill =   mpa)) + 
  #   geom_tile()
  
  mpa_sim <- simmar(
    fauna = fauna,
    fleets = tmp_fleets,
    years = years,
    mpas = list(locations = mpa_locations,
                mpa_year = floor(years * .5))
  )
  
  processed_marlin <-
    process_marlin(
      sim = mpa_sim,
      time_step = time_step,
      steps_to_keep = last(names(sim)),
      keep_age = FALSE
    )
  
  fauna_yield <- processed_marlin$fleets %>% 
    dplyr::ungroup() %>% 
    dplyr::summarise(yield = sum(catch, na.rm = TRUE))
  
  out <- fauna_yield$yield / msy
  
}

mpa_test_grid <- mpa_test_grid %>% 
  mutate(mpa_yields = map2_dbl(mpa_size, f_fmsy, mpa_test, fauna = fauna, fleets = fleets, random_mpas = FALSE))


 mpa_test_grid %>% 
  ggplot(aes(mpa_size, mpa_yields, color = factor(f_fmsy))) + 
  geom_hline(aes(yintercept = 0)) +
  geom_line(size = 2) + 
  theme_marlin() + 
  scale_x_continuous(labels = scales::percent, name = "% in MPA") + 
  scale_y_continuous(name = "Catch / MSY") + 
  scale_color_viridis_d()


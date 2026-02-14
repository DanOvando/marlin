library(marlin)

library(tidyverse)

resolution <- c(42,14)

sigma_r <-
max_shark_age <- 25

ages <- seq(0, 25, by = 1 / seasons)

length_at_age <- 216.4 * (1 - exp(-0.148 * ((ages) - -1.778991)))

pups_at_age <- pmax(0,-8.6 + 0.098 * length_at_age)

home_ranges = c("low" = 2, "medium" = 20, "high" = 200)

sigma_rec <- 0.2


critter <- marlin::create_critter(
  query_fishlife = FALSE,
  scientific_name = "Carcharhinus falciformis",
  common_name = "silky shark",
  m = 0.21,
  lorenzen_c = -1,
  linf = 216.4,
  vbk = 0.148,
  t0 = -1.778991, # 50cm at birth
  min_age = 0,
  max_age = 25,
  age_mature = 6,
  weight_a = 0.0000273,
  weight_b = 2.86,
  steepness = 0.3,
  b0 = 1e4,
  fec_at_age = pups_at_age,
  adult_home_range = home_ranges["high"],
  recruit_home_range = home_ranges["low"],
  resolution = resolution,
  patch_area = patch_area,
  seasons = seasons,
  habitat = shark_habitat,
  recruit_habitat = shark_recruit_habitat,
  density_dependence = "post_dispersal",
  fec_expo = 1,
  sigma_rec = sigma_rec,
  ac_rec = ac_rec)


critter <- marlin::create_critter(
  query_fishlife = FALSE,
  scientific_name = "Thunnus albacares",
  common_name = "yellowfin tuna",
  m = 0.3,
  growth_model = "growth_cessation",
  lorenzen_c = -1,
  l0 = 18.85,
  rmax = 37.24,
  k =  0.89,
  t50 = 4.57,
  min_age = 0,
  max_age = 15,
  age_mature = 3,
  weight_a = 0.00004,
  weight_b = 2.86,
  habitat = tuna_habitat,
  recruit_habitat = tuna_recruit_habitat,
  adult_home_range = home_ranges["high"],
  recruit_home_range = home_ranges["high"],
  density_dependence = "global_habitat",
  fec_expo = 1.3,
  steepness =  0.9,
  b0 = 1e6,
  resolution = resolution,
  patch_area = patch_area,
  sigma_rec = sigma_rec,
  ac_rec = ac_rec)


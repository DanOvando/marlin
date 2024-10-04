library(tidyverse)

library(marlin)

resolution <- 4

years <- 4

seasons <- 52

time_step <- 1 / seasons

fauna <-
  list(
    "squishy" = create_critter(
      common_name =  "squishy",
      growth_model = "power",
      length_a = .1,
      length_b = 3,
      cv_len = 0.4,
      m = 0.2,
      max_age = 5,
      age_mature = 4,
      delta_mature = .1,
      weight_a = 2,
      weight_b = 1,
      lorenzen_m = TRUE,
      adult_diffusion = .1,
      recruit_diffusion = 10,
      fished_depletion = .25,
      resolution = resolution,
      sigma_rec = 0.5,
      query_fishlife = FALSE,
      seasons = seasons
    )
  )


fleets <- list(
  "artisanal" = create_fleet(list(
    "squishy" = Metier$new(
      critter = fauna$squishy,
      p_explt = 0.5,
      sel_unit = "p_of_mat",
      sel_start = 0.5,
      sel_delta = 0.1
    )
  ), resolution = resolution),
  "commercial" = create_fleet(list(
    "squishy" = Metier$new(
      critter = fauna$squishy,
      p_explt = 0.5,
      sel_unit = "length",
      sel_start = 5,
      sel_delta = .01
    )
  ), resolution = resolution)
)


fauna$squishy$plot()

sels <- data.frame(
  artisinal = fleets$artisanal$metiers$squishy$sel_at_length
  ,
  commercial = fleets$commercial$metiers$squishy$sel_at_length,
  length = as.numeric(colnames(fauna$squishy$length_at_age_key))
)

sels |> 
  pivot_longer(-length, names_to = "fleet", values_to = "selectivity")  |> 
  ggplot(aes(length, selectivity, color = fleet)) + 
  geom_line()

fleets <- tune_fleets(fauna, fleets, tune_type = "depletion", years = 10)

squishy_sim <- simmar(fauna, fleets, years = years)

processed_squishy <- process_marlin(sim = squishy_sim, time_step = time_step)

plot_marlin(processed_squishy)

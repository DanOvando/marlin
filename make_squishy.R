library(tidyverse)

library(marlin)

library(ggridges)

library(tictoc)

library(gganimate)
# check on plus group issues
# add in basically mortality as mirror to fraction mature
# tweek selecitivy or q to allow for it to be time specific
# check if max maturity can be less than 1
# add effort as a management group
resolution <- 5

patch_area <- 1

kp <- .01

years <- 10

seasons <- 52 # weekly time step

time_step <- 1 / seasons

max_age <- 1.5

length_bin_width = 0.1

length_a <- 0.1

length_b <- 3

lorenzen_c <- -1

m <- 0.4

t0 <- -.5

ages <- seq(0,max_age, by = time_step)

length_at_age <- length_a * (ages - t0)^length_b

max_length <- max(length_at_age)

m_at_age <- m * (length_at_age / max(length_at_age))^lorenzen_c

critter_correlations <- matrix(c(1, -.8, -.8, 1), nrow = 2, byrow = TRUE)

species_distributions <- sim_habitat(critters = c("squishy", "squashy"), resolution = resolution, patch_area = patch_area, kp = kp, output = "list")

# m_at_age <- rep(0.2, length(ages))

fauna <-
  list(
    "squishy" = create_critter(
      spawning_seasons = c(1:20),
      habitat = species_distributions$critter_distributions$squishy,
      common_name =  "squishy",
      growth_model = "power",
      length_a = length_a,
      length_b = length_b,
      length_bin_width = length_bin_width,
      t0 = t0,
      cv_len = 0.6,
      m_at_age = m_at_age,
      max_age = max_age,
      age_mature = max_age - .5,
      semelparous = TRUE,
      delta_mature = .1,
      weight_a = 2,
      weight_b = 1,
      adult_diffusion = .1,
      recruit_diffusion = 10,
      init_explt = 1.2,
      resolution = resolution,
      sigma_rec = 0,
      ac_rec = .3,
      query_fishlife = FALSE,
      seasons = seasons
    ),
    "squashy" = create_critter(
      spawning_seasons = c(21:40),
      habitat = species_distributions$critter_distributions$squashy,
      common_name =  "squashy",
      growth_model = "power",
      length_a = length_a,
      length_b = 2,
      length_bin_width = length_bin_width,
      t0 = t0,
      cv_len = 0.6,
      m_at_age = m_at_age,
      max_age = max_age,
      age_mature = max_age  / 2,
      semelparous = FALSE,
      delta_mature = .1,
      weight_a = 2,
      weight_b = 1,
      adult_diffusion = .1,
      recruit_diffusion = 10,
      init_explt = .6,
      resolution = resolution,
      sigma_rec = .6,
      ac_rec = .3,
      query_fishlife = FALSE,
      seasons = seasons
    )
  )

fauna <-
  list(
    "squishy" = create_critter(
      spawning_seasons = c(1:20),
      habitat = species_distributions$critter_distributions$squishy,
      common_name =  "squishy",
      growth_model = "power",
      length_a = length_a,
      length_b = length_b,
      length_bin_width = length_bin_width,
      t0 = t0,
      cv_len = 0.6,
      m_at_age = m_at_age,
      max_age = max_age,
      age_mature = max_age - .5,
      semelparous = TRUE,
      delta_mature = .1,
      weight_a = 2,
      weight_b = 1,
      adult_diffusion = .1,
      recruit_diffusion = 10,
      init_explt = 1.2,
      resolution = resolution,
      sigma_rec = 0,
      ac_rec = .3,
      query_fishlife = FALSE,
      seasons = seasons
    ),
    "squashy" = create_critter(
      spawning_seasons = c(1:20),
      habitat = species_distributions$critter_distributions$squishy,
      common_name =  "squishy",
      growth_model = "power",
      length_a = length_a,
      length_b = length_b,
      length_bin_width = length_bin_width,
      t0 = t0,
      cv_len = 0.6,
      m_at_age = m_at_age,
      max_age = max_age,
      age_mature = max_age - .5,
      semelparous = FALSE,
      delta_mature = .1,
      weight_a = 2,
      weight_b = 1,
      adult_diffusion = .1,
      recruit_diffusion = 10,
      init_explt = 1.2,
      resolution = resolution,
      sigma_rec = 0,
      ac_rec = .3,
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
      sel_start = 0.1,
      sel_delta = 0.1
    ),
    "squashy" = Metier$new(
      critter = fauna$squashy,
      p_explt = 0.5,
      sel_unit = "p_of_mat",
      sel_start = 0.1,
      sel_delta = 0.1
    )
  ), resolution = resolution),
  "commercial" = create_fleet(list(
    "squishy" = Metier$new(
      critter = fauna$squishy,
      p_explt = 0.5,
      sel_unit = "length",
      sel_start = 0.25 * max_length,
      sel_delta = .01
    ),
    "squashy" = Metier$new(
      critter = fauna$squashy,
      p_explt = 0.5,
      sel_unit = "p_of_mat",
      sel_start = 0.5,
      sel_delta = 0.1
    )
  ), resolution = resolution)
)


fauna$squishy$plot()

fauna$squashy$plot()


sels <- data.frame(
  artisinal = fleets$artisanal$metiers$squishy$sel_at_length,
  commercial = fleets$commercial$metiers$squishy$sel_at_length,
  length = as.numeric(colnames(fauna$squishy$length_at_age_key))
)

sels |>
  pivot_longer(-length, names_to = "fleet", values_to = "selectivity")  |>
  ggplot(aes(length, selectivity, color = fleet)) +
  geom_line()

fleets <- tune_fleets(fauna, fleets, tune_type = "explt", years = 20)

tic()
squishy_sim <- simmar(fauna, fleets, years = years, cor_rec = critter_correlations)
toc()
processed_squishy <- process_marlin(sim = squishy_sim, time_step = time_step)

processed_squishy$fauna |>
  group_by(step, age, mean_length) |>
  summarise(number = sum(n)) |>
  group_by(step) |>
  mutate(pn = number / sum(number)) |>
  ungroup() |>
  ggplot(aes(mean_length,factor(step),height = pn)) +
  geom_ridgeline(stat = "identity", color = "transparent", fill = "red", scale = 10) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


processed_squishy$fauna |>
  filter(age == min(age)) |>
  ggplot(aes(step, n, color = critter)) +
  geom_line()

processed_squishy$fauna |>
  filter(step == max(step)) |>
  group_by(x,y,critter) |>
  summarise(n = sum(n)) |>
  group_by(critter) |>
  mutate(sn = as.numeric(scale(n))) |>
  ggplot(aes(x,y,fill = sn)) +
  geom_tile() +
  facet_wrap(~critter) +
  scale_fill_viridis_c()

plot_marlin(processed_squishy,fauna = fauna, plot_var = "n", max_scale = FALSE)

plot_marlin(processed_squishy,fauna = fauna, plot_var = "c", max_scale = FALSE)

plot_marlin(processed_squishy,fauna = fauna, plot_var = "b", max_scale = FALSE)

plot_marlin(processed_squishy,fauna = fauna, plot_var = "n", max_scale = FALSE, plot_type = "age", steps_to_plot = max(processed_squishy$fauna$step))


a = processed_squishy$fauna |>
  group_by(year, step, age, critter) |>
  summarise(n = sum(n)) |>
  ungroup() |>
  group_by(critter, year, step) |>
  mutate(sn = n / max(n)) |>
  group_by(critter, year) |>
  mutate(n = n / max(n)) |>
  ungroup()

b = processed_squishy$fauna |>
  group_by(year, age, critter) |>
  summarise(n = sum(n)) |>
  group_by(critter, year) |>
  mutate(yn = n / max(n)) |>
  ungroup()

mogive <- map(fauna, \(x) tibble(age = x$ages, maturity = x$maturity_at_age, semelparous = x$semelparous)) |>
  list_rbind(names_to = "critter")

b |>
  filter(year == max(year)) |>
  ggplot() +
  geom_col(aes(age, yn, fill = "Age Composition")) +
  geom_line(data = mogive,aes(age, maturity, linetype = semelparous,color = "Maturity At Age")) +
  facet_wrap(~critter) +
  scale_fill_manual(values = "blue") +
  scale_y_continuous()


a |>
  filter(year < 4) |>
  ggplot() +
  geom_density_ridges(aes(age, step, height = sn, group = step, fill = step),
                      stat = "identity",
                      color = "transparent",
                      scale = 4,
                      alpha = 0.5) +
  facet_wrap( ~ critter, scales = "free_x") +
  scale_fill_viridis_c()




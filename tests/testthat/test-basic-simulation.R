test_that("basic simulation runs and produces expected output structure", {
  resolution <- 5
  seasons <- 4
  time_step <- 1 / seasons
  years <- 3

  species_distributions <- sim_habitat(
    critters = c("critter_a", "critter_b"),
    resolution = resolution,
    patch_area = 1,
    kp = 0.5,
    output = "list"
  )

  fauna <- list(
    "critter_a" = create_critter(
      common_name = "critter_a",
      habitat = species_distributions$critter_distributions$critter_a,
      max_age = 5,
      linf = 50,
      vbk = 0.3,
      t0 = -0.5,
      m = 0.2,
      weight_a = 0.01,
      weight_b = 3,
      age_mature = 2,
      resolution = resolution,
      seasons = seasons,
      init_explt = 0.2,
      sigma_rec = 0,
      adult_diffusion = 1,
      recruit_diffusion = 5,
      query_fishlife = FALSE
    ),
    "critter_b" = create_critter(
      common_name = "critter_b",
      habitat = species_distributions$critter_distributions$critter_b,
      max_age = 10,
      linf = 80,
      vbk = 0.2,
      t0 = -0.5,
      m = 0.15,
      weight_a = 0.005,
      weight_b = 3,
      age_mature = 3,
      resolution = resolution,
      seasons = seasons,
      init_explt = 0.1,
      sigma_rec = 0,
      adult_diffusion = 1,
      recruit_diffusion = 5,
      query_fishlife = FALSE
    )
  )

  fleets <- list(
    "fleet_one" = create_fleet(list(
      "critter_a" = Metier$new(
        critter = fauna$critter_a,
        p_explt = 1,
        sel_unit = "p_of_mat",
        sel_start = 0.1,
        sel_delta = 0.1
      ),
      "critter_b" = Metier$new(
        critter = fauna$critter_b,
        p_explt = 1,
        sel_unit = "p_of_mat",
        sel_start = 0.1,
        sel_delta = 0.1
      )
    ), resolution = resolution)
  )

  fleets <- tune_fleets(fauna, fleets, tune_type = "explt", years = 20)

  sim <- simmar(fauna, fleets, years = years)

  # Simulation returns a list with one entry per time step
  expect_type(sim, "list")
  expected_steps <- years * seasons
  expect_length(sim, expected_steps)

  # Each step contains results for both critters
  first_step <- sim[[1]]
  expect_true("critter_a" %in% names(first_step))
  expect_true("critter_b" %in% names(first_step))

  # Population arrays have correct dimensions [patches x ages]
  n_patches <- resolution^2
  n_ages_a <- length(fauna$critter_a$ages)
  expect_equal(nrow(first_step$critter_a$n_p_a), n_patches)
  expect_equal(ncol(first_step$critter_a$n_p_a), n_ages_a)

  # Biomass should be non-negative everywhere
  last_step <- sim[[expected_steps]]
  expect_true(all(last_step$critter_a$b_p_a >= 0))
  expect_true(all(last_step$critter_b$b_p_a >= 0))

  # process_marlin returns fauna and fleets data.tables
  processed <- process_marlin(sim, time_step = time_step)
  expect_true("fauna" %in% names(processed))
  expect_true("fleets" %in% names(processed))
  expect_s3_class(processed$fauna, "data.frame")
  expect_s3_class(processed$fleets, "data.frame")

  # Processed output should contain both critters and the fleet
  expect_true(all(c("critter_a", "critter_b") %in% unique(processed$fauna$critter)))
  expect_true("fleet_one" %in% unique(processed$fleets$fleet))
})

test_that("simulation with MPA reduces catch in closed patches", {
  resolution <- 4
  seasons <- 2
  time_step <- 1 / seasons
  years <- 5

  species_distributions <- sim_habitat(
    critters = "species_a",
    resolution = resolution,
    patch_area = 1,
    kp = 0.5,
    output = "list"
  )

  fauna <- list(
    "species_a" = create_critter(
      common_name = "species_a",
      habitat = species_distributions$critter_distributions$species_a,
      max_age = 8,
      linf = 60,
      vbk = 0.25,
      t0 = -0.5,
      m = 0.2,
      weight_a = 0.01,
      weight_b = 3,
      age_mature = 3,
      resolution = resolution,
      seasons = seasons,
      init_explt = 0.3,
      sigma_rec = 0,
      adult_diffusion = 0.1,
      recruit_diffusion = 1,
      query_fishlife = FALSE
    )
  )

  fleets <- list(
    "fleet_a" = create_fleet(list(
      "species_a" = Metier$new(
        critter = fauna$species_a,
        p_explt = 1,
        sel_unit = "p_of_mat",
        sel_start = 0.1,
        sel_delta = 0.1
      )
    ), resolution = resolution)
  )

  fleets <- tune_fleets(fauna, fleets, tune_type = "explt", years = 20)

  # Create MPA: close half the patches (rows 1-2 closed, rows 3-4 open)
  mpa_locations <- tidyr::expand_grid(x = 1:resolution, y = 1:resolution) |>
    dplyr::mutate(mpa = y <= 2)

  manager <- list(mpas = list(locations = mpa_locations, mpa_year = 1))

  sim_mpa <- simmar(fauna, fleets, years = years, manager = manager)

  processed <- process_marlin(sim_mpa, time_step = time_step)

  # Catch in MPA patches should be zero
  fleet_data <- processed$fleets
  last_year <- fleet_data[fleet_data$step == max(fleet_data$step), ]

  mpa_catch <- sum(last_year$catch[last_year$y <= 2], na.rm = TRUE)
  open_catch <- sum(last_year$catch[last_year$y > 2], na.rm = TRUE)

  expect_equal(mpa_catch, 0)
  expect_true(open_catch > 0)
})

test_that("unfished simulation preserves biomass near B0", {
  resolution <- 4
  seasons <- 2
  time_step <- 1 / seasons
  years <- 5

  species_distributions <- sim_habitat(
    critters = "sole",
    resolution = resolution,
    patch_area = 1,
    kp = 0.5,
    output = "list"
  )

  fauna <- list(
    "sole" = create_critter(
      common_name = "sole",
      habitat = species_distributions$critter_distributions$sole,
      max_age = 10,
      linf = 40,
      vbk = 0.2,
      t0 = -0.5,
      m = 0.2,
      weight_a = 0.01,
      weight_b = 3,
      age_mature = 3,
      resolution = resolution,
      seasons = seasons,
      init_explt = 0,
      sigma_rec = 0,
      adult_diffusion = 1,
      recruit_diffusion = 5,
      query_fishlife = FALSE
    )
  )

  # Fleet with zero exploitation
  fleets <- list(
    "ghost" = create_fleet(list(
      "sole" = Metier$new(
        critter = fauna$sole,
        p_explt = 0,
        sel_unit = "p_of_mat",
        sel_start = 0.1,
        sel_delta = 0.1
      )
    ), resolution = resolution)
  )

  fleets <- tune_fleets(fauna, fleets, tune_type = "explt", years = 20)

  sim <- simmar(fauna, fleets, years = years)

  # Total biomass at start and end should be similar (no fishing, no rec variability)
  b_start <- sum(sim[[1]]$sole$b_p_a)
  b_end <- sum(sim[[length(sim)]]$sole$b_p_a)

  expect_equal(b_end, b_start, tolerance = 0.05)

  # SSB should also be near unfished
  ssb0 <- fauna$sole$ssb0
  ssb_end <- sum(sim[[length(sim)]]$sole$ssb_p_a)
  expect_equal(ssb_end, ssb0, tolerance = 0.05)
})

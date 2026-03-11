test_that("query_fishlife = FALSE with all required params succeeds", {
  resolution <- 4
  seasons <- 2

  species_distributions <- sim_habitat(
    critters = "test_sp",
    resolution = resolution,
    patch_area = 1,
    kp = 0.5,
    output = "list"
  )

  # Should work without any network calls
  critter <- create_critter(
    common_name = "test_sp",
    habitat = species_distributions$critter_distributions$test_sp,
    query_fishlife = FALSE,
    linf = 50,
    vbk = 0.3,
    t0 = -0.5,
    m = 0.2,
    weight_a = 0.01,
    weight_b = 3,
    age_mature = 2,
    max_age = 10,
    resolution = resolution,
    seasons = seasons,
    init_explt = 0.1,
    sigma_rec = 0
  )

  expect_s3_class(critter, "Fish")
  expect_equal(critter$linf, 50)
  expect_equal(critter$vbk, 0.3)
  expect_equal(critter$weight_a, 0.01)
  expect_equal(critter$weight_b, 3)
})

test_that("query_fishlife = FALSE errors when required params are missing", {
  resolution <- 4
  seasons <- 2

  species_distributions <- sim_habitat(
    critters = "test_sp",
    resolution = resolution,
    patch_area = 1,
    kp = 0.5,
    output = "list"
  )

  # Missing linf and vbk for von_bertalanffy
  expect_error(
    create_critter(
      common_name = "test_sp",
      habitat = species_distributions$critter_distributions$test_sp,
      query_fishlife = FALSE,
      m = 0.2,
      weight_a = 0.01,
      weight_b = 3,
      age_mature = 2,
      resolution = resolution,
      seasons = seasons
    ),
    "linf"
  )

  # Missing m
  expect_error(
    create_critter(
      common_name = "test_sp",
      habitat = species_distributions$critter_distributions$test_sp,
      query_fishlife = FALSE,
      linf = 50,
      vbk = 0.3,
      weight_a = 0.01,
      weight_b = 3,
      age_mature = 2,
      resolution = resolution,
      seasons = seasons
    ),
    "m"
  )

  # Missing maturity specification
  expect_error(
    create_critter(
      common_name = "test_sp",
      habitat = species_distributions$critter_distributions$test_sp,
      query_fishlife = FALSE,
      linf = 50,
      vbk = 0.3,
      m = 0.2,
      weight_a = 0.01,
      weight_b = 3,
      resolution = resolution,
      seasons = seasons
    ),
    "age_mature"
  )

  # Missing weight_a and weight_b
  expect_error(
    create_critter(
      common_name = "test_sp",
      habitat = species_distributions$critter_distributions$test_sp,
      query_fishlife = FALSE,
      linf = 50,
      vbk = 0.3,
      m = 0.2,
      age_mature = 2,
      resolution = resolution,
      seasons = seasons
    ),
    "weight_a"
  )
})

test_that("power growth model does not require linf/vbk", {
  resolution <- 4
  seasons <- 2

  species_distributions <- sim_habitat(
    critters = "test_sp",
    resolution = resolution,
    patch_area = 1,
    kp = 0.5,
    output = "list"
  )

  # power growth model doesn't need linf or vbk
  critter <- create_critter(
    common_name = "test_sp",
    habitat = species_distributions$critter_distributions$test_sp,
    query_fishlife = FALSE,
    growth_model = "power",
    length_a = 0.1,
    length_b = 3,
    length_bin_width = 0.1,
    t0 = -0.5,
    m = 0.4,
    max_age = 1.5,
    weight_a = 2,
    weight_b = 1,
    age_mature = 1,
    resolution = resolution,
    seasons = seasons,
    init_explt = 0.1,
    sigma_rec = 0
  )

  expect_s3_class(critter, "Fish")
  expect_equal(critter$growth_model, "power")
})

test_that("m_at_age bypasses the m requirement", {
  resolution <- 4
  seasons <- 2
  max_age <- 5
  time_step <- 1 / seasons

  species_distributions <- sim_habitat(
    critters = "test_sp",
    resolution = resolution,
    patch_area = 1,
    kp = 0.5,
    output = "list"
  )

  ages <- seq(0, max_age, by = time_step)
  m_at_age_vec <- rep(0.2, length(ages))

  # m = NA but m_at_age supplied should work
  critter <- create_critter(
    common_name = "test_sp",
    habitat = species_distributions$critter_distributions$test_sp,
    query_fishlife = FALSE,
    linf = 50,
    vbk = 0.3,
    weight_a = 0.01,
    weight_b = 3,
    age_mature = 2,
    max_age = max_age,
    m_at_age = m_at_age_vec,
    resolution = resolution,
    seasons = seasons,
    init_explt = 0.1,
    sigma_rec = 0
  )

  expect_s3_class(critter, "Fish")
})

test_that("NULL defaults work: no common_name, no spawning_seasons, no recruit_habitat", {
  resolution <- 4
  seasons <- 2

  species_distributions <- sim_habitat(
    critters = "test_sp",
    resolution = resolution,
    patch_area = 1,
    kp = 0.5,
    output = "list"
  )

  # Omit common_name, spawning_seasons, recruit_habitat (all NULL by default)
  critter <- create_critter(
    habitat = species_distributions$critter_distributions$test_sp,
    query_fishlife = FALSE,
    linf = 50,
    vbk = 0.3,
    m = 0.2,
    weight_a = 0.01,
    weight_b = 3,
    age_mature = 2,
    max_age = 10,
    resolution = resolution,
    seasons = seasons,
    init_explt = 0.1,
    sigma_rec = 0
  )

  expect_s3_class(critter, "Fish")
  # common_name should be NULL, not NA
  expect_null(critter$common_name)
  # spawning_seasons should have been filled to 1:seasons
  expect_equal(critter$spawning_seasons, 1:seasons)
})

test_that("error message lists all missing params at once", {
  resolution <- 4
  seasons <- 2

  species_distributions <- sim_habitat(
    critters = "test_sp",
    resolution = resolution,
    patch_area = 1,
    kp = 0.5,
    output = "list"
  )

  # Missing everything: linf, vbk, weight_a, weight_b, m, maturity
  err <- expect_error(
    create_critter(
      common_name = "test_sp",
      habitat = species_distributions$critter_distributions$test_sp,
      query_fishlife = FALSE,
      resolution = resolution,
      seasons = seasons
    )
  )

  # Should mention all missing params in one message
  expect_match(err$message, "linf")
  expect_match(err$message, "vbk")
  expect_match(err$message, "weight_a")
  expect_match(err$message, "weight_b")
  expect_match(err$message, "m")
  expect_match(err$message, "age_mature")
  expect_match(err$message, "query_fishlife")
})

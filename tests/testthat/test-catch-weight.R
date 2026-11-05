test_that("Baranov catch is weighed at the age the fish were caught at", {
  # One swim() step with known starting numbers and F and no movement, compared
  # with the Baranov catch computed by hand. Guards the 1.0.0-1.1.3 bug where
  # age a-1 catch was weighed with weight_at_age[a].
  critter <- create_critter(
    common_name = "test fish",
    query_fishlife = FALSE,
    linf = 100,
    vbk = 0.2,
    t0 = -0.5,
    max_age = 12,
    age_mature = 3,
    weight_a = 1e-5,
    weight_b = 3,
    m = 0.2,
    steepness = 0.8,
    adult_home_range = 1,
    recruit_home_range = 1,
    resolution = c(2, 2),
    patch_area = 1,
    seasons = 2,
    density_dependence = "global_habitat"
  )
  n_ages <- length(critter$ages)
  n_patches <- critter$patches
  set.seed(1)
  start_n <- matrix(runif(n_patches * n_ages, 100, 1000), nrow = n_patches)
  f_at_age <- matrix(runif(n_patches * n_ages, 0, 0.5), nrow = n_patches)
  identity_move <- methods::as(
    methods::as(Matrix::Diagonal(n_patches), "generalMatrix"),
    "CsparseMatrix"
  )
  no_movement <- rep(list(identity_move), length(critter$movement_seasons))

  pop <- critter$swim(
    season = 1,
    adult_movement = no_movement,
    f_p_a = f_at_age,
    last_n_p_a = start_n,
    rec_devs = rep(1, n_patches)
  )

  z_at_age <- sweep(f_at_age, 2, critter$m_at_age, "+")
  catch_numbers <- (f_at_age / z_at_age) * start_n *
    (1 - exp(-critter$time_step * z_at_age))
  expected_catch <- sweep(catch_numbers, 2, critter$weight_at_age, "*")

  expect_equal(unname(pop$c_p_a), unname(expected_catch), tolerance = 1e-12)
})

test_that("catch in numbers is returned and is catch weight divided by weight at age", {
  critter <- create_critter(
    common_name = "test fish",
    query_fishlife = FALSE,
    linf = 100,
    vbk = 0.2,
    t0 = -0.5,
    max_age = 12,
    age_mature = 3,
    weight_a = 1e-5,
    weight_b = 3,
    m = 0.2,
    steepness = 0.8,
    adult_home_range = 1,
    recruit_home_range = 1,
    resolution = c(2, 2),
    patch_area = 1,
    seasons = 2,
    density_dependence = "global_habitat"
  )
  n_ages <- length(critter$ages)
  n_patches <- critter$patches
  set.seed(1)
  start_n <- matrix(runif(n_patches * n_ages, 100, 1000), nrow = n_patches)
  f_at_age <- matrix(runif(n_patches * n_ages, 0, 0.5), nrow = n_patches)
  identity_move <- methods::as(
    methods::as(Matrix::Diagonal(n_patches), "generalMatrix"),
    "CsparseMatrix"
  )
  pop <- critter$swim(
    season = 1,
    adult_movement = rep(list(identity_move), length(critter$movement_seasons)),
    f_p_a = f_at_age,
    last_n_p_a = start_n,
    rec_devs = rep(1, n_patches)
  )

  z_at_age <- sweep(f_at_age, 2, critter$m_at_age, "+")
  expected_numbers <- (f_at_age / z_at_age) * start_n *
    (1 - exp(-critter$time_step * z_at_age))

  expect_equal(unname(pop$cn_p_a), unname(expected_numbers), tolerance = 1e-12)
  expect_equal(
    unname(pop$c_p_a),
    unname(sweep(pop$cn_p_a, 2, critter$weight_at_age, "*")),
    tolerance = 1e-12
  )
})

test_that("simmar and process_marlin carry catch in numbers alongside weight", {
  critter <- create_critter(
    common_name = "test_fish",
    query_fishlife = FALSE,
    linf = 100,
    vbk = 0.2,
    t0 = -0.5,
    max_age = 12,
    age_mature = 3,
    weight_a = 1e-5,
    weight_b = 3,
    m = 0.2,
    resolution = 3,
    seasons = 2,
    init_explt = 0.2,
    sigma_rec = 0,
    adult_diffusion = 1,
    recruit_diffusion = 5
  )
  fauna <- list("test_fish" = critter)
  fleets <- list(
    "fleet_one" = create_fleet(list(
      "test_fish" = Metier$new(
        critter = critter,
        p_explt = 1,
        sel_unit = "p_of_mat",
        sel_start = 0.1,
        sel_delta = 0.1
      )
    ), resolution = 3)
  )
  fleets <- tune_fleets(fauna, fleets, tune_type = "explt", years = 20)
  sim <- simmar(fauna, fleets, years = 4)
  step_out <- sim[[length(sim) - 1]][["test_fish"]]

  expect_equal(
    unname(step_out$c_p_a),
    unname(sweep(step_out$cn_p_a, 2, critter$weight_at_age, "*")),
    tolerance = 1e-12
  )
  expect_equal(sum(step_out$cn_p_a_fl), sum(step_out$cn_p_a), tolerance = 1e-8)

  processed <- process_marlin(sim, time_step = 0.5, keep_age = TRUE)
  expect_true("cn" %in% names(processed$fauna))
  expect_true("catch_numbers" %in% names(processed$fleets))
})

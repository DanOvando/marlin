# Habitat matrices are what-you-see-is-what-you-get: matrix row 1 is the TOP of
# the map. A surface handed to create_critter() and the identical surface handed
# to simmar(habitat = ) must therefore be read the same way.
#
# Every surface in these tests is deliberately y-ASYMMETRIC. A y-flip is
# undetectable on the symmetric/centred habitat used elsewhere in the suite,
# which is how simmar() shipped an as.vector() read that flipped every supplied
# surface about the y axis. Do not "simplify" these fixtures to a centred dnorm.

build_test_critter <- function(hab, resolution) {
  create_critter(
    query_fishlife = FALSE,
    common_name = "test_critter",
    m = 0.3,
    linf = 50,
    vbk = 0.3,
    t0 = -0.5,
    max_age = 10,
    age_mature = 2,
    weight_a = 0.01,
    weight_b = 3,
    habitat = hab,
    # flat recruit habitat isolates adult taxis from the r0s anchor
    recruit_habitat = matrix(1, nrow = resolution[2], ncol = resolution[1]),
    adult_home_range = 5,
    depletion = 0.4,
    density_dependence = "global_habitat",
    seasons = 1,
    init_explt = 0,
    explt_type = "f",
    sigma_rec = 0
  )
}

build_test_fleet <- function(critter, resolution) {
  fleets <- list(
    "test_fleet" = create_fleet(
      metiers = list("test_critter" = Metier$new(
        critter = critter,
        price = 10,
        sel_form = "logistic",
        sel_start = 0.1,
        sel_delta = 0.01,
        catchability = 0,
        p_explt = 1
      )),
      base_effort = prod(resolution),
      resolution = resolution
    )
  )
  tune_fleets(list("test_critter" = critter), fleets)
}

peak_patch <- function(sim) {
  proc <- process_marlin(sim = sim, time_step = 1, keep_age = FALSE)
  out <- proc$fauna[proc$fauna$step == max(proc$fauna$step), ]
  out <- stats::aggregate(b ~ x + y, data = out, FUN = sum)
  hottest <- out[which.max(out$b), ]
  c(x = hottest$x, y = hottest$y)
}

test_that("the two habitat readers agree, and differ only by a row reversal", {
  hab <- matrix(seq_len(50), nrow = 10, ncol = 5)

  fish_class_ordering <- as.data.frame(hab) |>
    dplyr::mutate(y = dplyr::n():1) |>
    tidyr::pivot_longer(-y, names_to = "x", values_to = "value") |>
    dplyr::mutate(x = match(x, unique(x))) |>
    dplyr::arrange(x, y) |>
    dplyr::pull(value)

  # what simmar() must do to match Fish$initialize()
  expect_equal(fish_class_ordering, as.vector(hab[nrow(hab):1, , drop = FALSE]))

  # and a bare as.vector() must NOT match, or this test has stopped biting
  expect_false(isTRUE(all.equal(fish_class_ordering, as.vector(hab))))
})

test_that("habitat in the top-left of the matrix puts fish in the top-left of the map", {
  resolution <- c(5, 10)

  hab <- matrix(0, nrow = resolution[2], ncol = resolution[1])
  hab[1, 1] <- 5 # top-left as literally displayed

  critter <- build_test_critter(hab, resolution)
  fauna <- list("test_critter" = critter)
  fleets <- build_test_fleet(critter, resolution)

  # top-left of the map == x = 1, y = ny
  expect_equal(
    peak_patch(simmar(fauna = fauna, fleets = fleets, years = 15)),
    c(x = 1, y = resolution[2]),
    ignore_attr = TRUE
  )
})

test_that("passing simmar() the habitat create_critter() already has is a no-op", {
  resolution <- c(5, 10)

  # y-asymmetric on purpose: a centred surface cannot detect a y-flip
  hab_df <- tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2])
  hab_df$habitat <- stats::dnorm(
    sqrt((hab_df$x - 2)^2 + (hab_df$y - 3)^2), 0, 3
  ) * 20

  hab <- as.matrix(
    tidyr::pivot_wider(hab_df, names_from = x, values_from = habitat)[, -1]
  )

  critter <- build_test_critter(hab, resolution)
  fauna <- list("test_critter" = critter)
  fleets <- build_test_fleet(critter, resolution)

  without_habitat_arg <- simmar(fauna = fauna, fleets = fleets, years = 10)

  with_habitat_arg <- simmar(
    fauna = fauna,
    fleets = fleets,
    habitat = list("test_critter" = rep(list(hab), 10)),
    years = 10
  )

  expect_equal(
    peak_patch(without_habitat_arg),
    peak_patch(with_habitat_arg)
  )

  expect_equal(
    without_habitat_arg[[length(without_habitat_arg)]]$test_critter$b_p_a,
    with_habitat_arg[[length(with_habitat_arg)]]$test_critter$b_p_a,
    tolerance = 1e-8
  )
})

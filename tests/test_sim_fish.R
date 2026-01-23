library(Rcpp)
library(RcppEigen)
library(here)

Rcpp::sourceCpp(here("src","sim_fish_old.cpp"))
Rcpp::sourceCpp(here("src","sim_fish_new.cpp"))


compare_sim_outputs <- function(old, new, tol = 1e-10) {
  stopifnot(is.list(old), is.list(new))

  # compare a known subset (ignore "tmppop" because it can differ in internal structure)
  keys <- intersect(names(old), names(new))
  keys <- setdiff(keys, "tmppop")

  diffs <- list()

  for (k in keys) {
    a <- old[[k]]
    b <- new[[k]]

    same <- isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE))
    if (!same) {
      diffs[[k]] <- list(
        old_class = class(a),
        new_class = class(b),
        old_dim = if (!is.null(dim(a))) dim(a) else length(a),
        new_dim = if (!is.null(dim(b))) dim(b) else length(b),
        max_abs_diff = if (is.numeric(a) && is.numeric(b)) max(abs(a - b), na.rm = TRUE) else NA_real_
      )
    }
  }

  if (length(diffs) == 0) {
    message("✅ Outputs match (within tolerance = ", tol, ") for all checked keys: ",
            paste(keys, collapse = ", "))
    return(invisible(TRUE))
  } else {
    message("❌ Differences found in: ", paste(names(diffs), collapse = ", "))
    print(diffs)
    return(invisible(diffs))
  }
}

make_test_inputs <- function(patches = 20, ages = 10, time_step = 0.25, seed = 1) {
  set.seed(seed)

  length_at_age <- seq(10, 10 + ages - 1)
  weight_at_age <- seq(1, ages) / ages + 0.5
  fec_at_age <- pmax(0, seq(-2, 2, length.out = ages))
  maturity_at_age <- plogis(seq(-3, 3, length.out = ages))

  semelparous <- FALSE

  f_p_a <- matrix(rexp(patches * ages, rate = 5), patches, ages)
  # keep F modest
  f_p_a <- pmin(f_p_a, 0.5)

  # movement: dense row-stochastic matrix
  M <- matrix(runif(patches * patches), patches, patches)
  M <- M / rowSums(M)

  # movement season mapping: single type used for all seasons
  movement_matrix <- list(M)
  movement_seasons <- list(1:round(1 / time_step))

  # recruit movement matrix: column-stochastic-ish operator for vectors
  Rm <- matrix(runif(patches * patches), patches, patches)
  Rm <- Rm / rowSums(Rm)

  last_n_p_a <- matrix(rexp(patches * ages, rate = 1), patches, ages)

  burn_steps <- 0L
  season <- 1L

  steepness <- 0.75
  r0s <- runif(patches, min = 50, max = 150)

  # test both: ssb0 known vs ssb0 = NA (but you said tuning rare)
  ssb0 <- 1e6
  ssb0_p <- rep(ssb0 / patches, patches)

  m_at_age <- runif(ages, min = 0.05, max = 0.3)

  tune_unfished <- FALSE
  rec_form <- "local_habitat"
  spawning_seasons <- c(1L)  # spawn in season 1
  rec_devs <- rep(1, patches)

  list(
    length_at_age = length_at_age,
    weight_at_age = weight_at_age,
    fec_at_age = fec_at_age,
    maturity_at_age = maturity_at_age,
    semelparous = semelparous,
    f_p_a = f_p_a,
    movement_matrix = movement_matrix,
    movement_seasons = movement_seasons,
    recruit_movement_matrix = Rm,
    last_n_p_a = last_n_p_a,
    patches = patches,
    burn_steps = burn_steps,
    time_step = time_step,
    season = season,
    steepness = steepness,
    r0s = r0s,
    ssb0 = ssb0,
    ssb0_p = ssb0_p,
    m_at_age = m_at_age,
    tune_unfished = tune_unfished,
    rec_form = rec_form,
    spawning_seasons = spawning_seasons,
    rec_devs = rec_devs
  )
}

run_one_test <- function(seed = 1, patches = 20, ages = 10, time_step = 0.25, tol = 1e-10) {
  inp <- make_test_inputs(patches = patches, ages = ages, time_step = time_step, seed = seed)

  old <- do.call(sim_fish_old, inp)
  new <- do.call(sim_fish_new, inp)

  compare_sim_outputs(old, new, tol = tol)
}

# One test
run_one_test(seed = 1)

# A small randomized battery
for (s in 1:20) {
  ok <- run_one_test(seed = s, patches = 30, ages = 12, time_step = 0.2, tol = 1e-10)
  if (!isTRUE(ok)) stop("Mismatch at seed = ", s)
}
message("🎉 All tests passed.")

test_ssb0_tuning <- function(seed = 1, tol = 1e-10) {
  inp <- make_test_inputs(patches = 15, ages = 10, time_step = 0.25, seed = seed)
  inp$ssb0 <- NA_real_
  inp$burn_steps <- 4L
  inp$tune_unfished <- TRUE  # your code expects this internally for burn

  old <- do.call(sim_fish_old, inp)
  new <- do.call(sim_fish_new, inp)

  compare_sim_outputs(old, new, tol = tol)
}

test_ssb0_tuning(seed = 1)


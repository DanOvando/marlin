# Per-unit objectives (rpue/cpue/ppue) are NA in patches where a fleet deployed
# no effort last step -- aggregate_yields() divides by effort and sets NA at
# exactly zero. NA means "no information about this patch", not "this patch is
# bad", and must never leak into other patches' velocities.
#
# Before the fix, mean()/max() without na.rm meant ONE unfished patch turned the
# whole fleet's effort vector into NA.

make_fleet <- function(responsiveness = 0.5, alloc = "rpue") {
  list(fl = list(spatial_allocation = alloc, responsiveness = responsiveness))
}

alloc_with <- function(effort, objective, responsiveness = 0.5) {
  effort <- matrix(effort, ncol = 1)
  objective <- matrix(objective, ncol = 1)
  colnames(effort) <- "fl"
  colnames(objective) <- "fl"

  allocate_effort(
    effort_by_patch = effort,
    buffet = list(rpue_p_fl = objective),
    fleets = make_fleet(responsiveness),
    open_patch = matrix(TRUE, nrow = nrow(effort), ncol = 1)
  )
}

test_that("an NA objective in one patch does not NA the whole fleet", {
  # patch 6 was unfished last step, so its rpue is NA
  res <- alloc_with(
    effort = c(1, 1, 1, 1, 1, 0),
    objective = c(5, 4, 3, 2, 1, NA)
  )

  expect_false(any(is.na(res$effort_new)))
  expect_equal(sum(res$effort_new), 5) # total effort conserved
})

test_that("NA objectives do not distort the velocities of observed patches", {
  # Adding an unobserved patch must not change how the observed patches are
  # ranked against each other -- the NA carries no information either way.
  without_na <- alloc_with(
    effort = c(1, 1, 1, 1),
    objective = c(4, 3, 2, 1)
  )
  with_na <- alloc_with(
    effort = c(1, 1, 1, 1, 0),
    objective = c(4, 3, 2, 1, NA)
  )

  expect_equal(
    order(without_na$effort_new[, 1], decreasing = TRUE),
    order(with_na$effort_new[1:4, 1], decreasing = TRUE)
  )
  expect_false(any(is.na(with_na$velocity_by_patch)))
})

test_that("an unobserved patch receives no velocity", {
  res <- alloc_with(
    effort = c(1, 1, 1, 1, 0),
    objective = c(4, 3, 2, 1, NA)
  )

  # no information => no push
  expect_equal(unname(res$velocity_by_patch[5, 1]), 0)
})

test_that("several NA patches alongside several distinct observed ones is safe", {
  # This is the configuration that used to blow up: >= 2 surviving patches with
  # distinct objectives means the flatness detector does NOT fire, so nothing
  # masked the NA propagation.
  n <- 50
  effort <- c(25, 25, rep(0, n - 2))
  objective <- c(10, 11, rep(NA, n - 2))

  res <- alloc_with(effort = effort, objective = objective)

  expect_false(any(is.na(res$effort_new)))
  expect_equal(sum(res$effort_new), 50)
})

test_that("a fleet that observed nothing anywhere holds position", {
  effort <- c(1, 2, 3, 4)
  res <- alloc_with(
    effort = effort,
    objective = rep(NA_real_, 4)
  )

  expect_false(any(is.na(res$effort_new)))
  expect_equal(as.numeric(res$effort_new), effort)
})

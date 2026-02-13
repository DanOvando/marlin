#' Go Fishing applies a given matrix of effort by patch and fleet to determine the
#' returns from fishing of that effort. This does not actually affect the population, this is
#' essentially an exploratory fishing step used to allocate fishing effort in space.
#'
#' @param e_p_fl a matrix of effort by patch and fleet
#' @param fauna a list of fauna objects
#' @param fleets a list of fleet objects
#' @param groupers the resolution to return the outputs at (used when output_format = "tidy")
#' @param output_format character: "tidy" (default) returns list of tidy data frames
#'   with species summed and effort joined; "matrix" returns list of patches x fleets
#'   matrices (species summed, no effort column). Matrix format is designed for
#'   direct use with \code{\link{allocate_effort}}.
#' @param n_p_a a list of length(fauna) with the numbers by patch and age to go fish
#'
#' @returns a list with six elements: raw catch (\code{c_p_fl}), revenue
#'   (\code{r_p_fl}), profit (\code{prof_p_fl}), and their per-unit-effort
#'   counterparts (\code{cpue_p_fl}, \code{rpue_p_fl}, \code{ppue_p_fl}).
#'   Format depends on \code{output_format}.
#' @export
#'
go_fish <- function(e_p_fl, fauna, n_p_a, fleets, groupers = c("fleet", "patch"),
                    output_format = c("tidy", "matrix")) {

  output_format <- match.arg(output_format)

  fauni <- names(fauna)

  fleet_names <- names(fleets)

  patches <- prod(fauna[[1]]$resolution)

  yields <- vector(mode = "list", length = length(fauna))

  names(yields) <- fauni

  current_season <- 1

  rec_devs <- rep(1, patches)

  for (f in seq_along(fauni)) {
    # setup total f and f by fleet --------------------------------------------

    critter <- fauni[f]

    last_n_p_a <- n_p_a[[critter]]

    ages <- length(fauna[[critter]]$ages)

    movement <- fauna[[critter]]$movement_matrix

    f_p_a <-
      matrix(0, nrow = patches, ncol = ages) # total fishing mortality by patch and age

    f_p_a_fl <-
      array(
        0,
        dim = c(patches, ages, length(fleets)),
        dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
      ) # storage for proportion of fishing mortality by patch, age, and fleet

    p_p_a_fl <-
      array(
        0,
        dim = c(patches, ages, length(fleets)),
        dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
      ) # storage for price by patch, age, and fleet

    for (l in seq_along(fleet_names)) {

      tmp <- fleets[[fleet_names[l]]]$metiers[[critter]]$vul_p_a

      f_p_a <-
        f_p_a + e_p_fl[, fleet_names[l]] * tmp

      f_p_a_fl[, , l] <-
        e_p_fl[, fleet_names[l]] * tmp

      p_p_a_fl[, , l] <- fleets[[l]]$metiers[[critter]]$price
    } # calculate cumulative f at age by patch

    f_p_a_fl <-
      f_p_a_fl / array(
        f_p_a,
        dim = c(patches, ages, length(fleets)),
        dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
      ) # f by patch, age, and fleet

    # apply fishing mortality to population -----------------------------------

    trial_run <- fauna[[critter]]$swim(
      season = current_season,
      adult_movement = movement,
      f_p_a = f_p_a,
      last_n_p_a = last_n_p_a,
      rec_devs = rec_devs
    )

    yields[[critter]] <- allocate_yields(
      f_p_a_fl = f_p_a_fl, e_p_fl = e_p_fl, p_p_a_fl = p_p_a_fl,
      critter = critter, pop = trial_run, fauna = fauna,
      fleets = fleets, patches = patches, ages = ages
    )
  } ## close fauni loop

  aggregate_yields(yields, e_p_fl, output_format = output_format, groupers = groupers)
}

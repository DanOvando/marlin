#' Exploratory Fishing Step (No Population Effect)
#'
#' @description
#' Applies a given matrix of effort by patch and fleet to determine the
#' potential returns from fishing without modifying the underlying population.
#' This is an exploratory step used internally by \code{\link{simmar}} to
#' build the "buffet" that \code{\link{allocate_effort}} uses for spatial
#' effort reallocation.
#'
#' @details
#' \code{go_fish} runs a single, movement-free time step (\code{move_fish = 0})
#' for each critter in \code{fauna} using the supplied population state
#' \code{n_p_a}, then aggregates catch, revenue, and profit across species
#' via \code{\link{aggregate_yields}}.
#'
#' The function is also useful for "what-if" analyses: you can probe the
#' revenue surface under different effort distributions before committing to a
#' simulation run.
#'
#' @param e_p_fl Numeric matrix of effort by patch (rows) and fleet (columns),
#'   with column names matching fleet names. Typically the current step's
#'   effort matrix from \code{\link{simmar}}.
#' @param fauna Named list of fauna objects from \code{\link{create_critter}}.
#' @param n_p_a Named list (one element per species in \code{fauna}) of
#'   numbers-at-age matrices \code{[patch, age]}, representing the current
#'   population state to fish against. Typically extracted from
#'   \code{sim[[step]][[critter]]$n_p_a}.
#' @param fleets Named list of fleet objects from \code{\link{create_fleet}}.
#' @param groupers Character vector. Grouping columns for tidy output format
#'   (default: \code{c("fleet", "patch")}). Ignored when
#'   \code{output_format = "matrix"}.
#' @param output_format Character. Output format:
#'   \describe{
#'     \item{\code{"tidy"}}{(default) Returns a list of tidy data frames
#'       with species summed across and effort joined.}
#'     \item{\code{"matrix"}}{Returns a list of patches x fleets matrices
#'       (species summed). Designed for direct use with
#'       \code{\link{allocate_effort}}. Substantially faster.}
#'   }
#'
#' @return A named list with six elements, each summed across species:
#' \describe{
#'   \item{\code{c_p_fl}}{Catch by patch and fleet}
#'   \item{\code{r_p_fl}}{Revenue by patch and fleet}
#'   \item{\code{prof_p_fl}}{Profit by patch and fleet}
#'   \item{\code{cpue_p_fl}}{Catch per unit effort (\code{NA} where effort = 0)}
#'   \item{\code{rpue_p_fl}}{Revenue per unit effort (\code{NA} where effort = 0)}
#'   \item{\code{ppue_p_fl}}{Profit per unit effort (\code{NA} where effort = 0)}
#' }
#' Format of each element depends on \code{output_format}.
#'
#' @seealso \code{\link{allocate_effort}}, \code{\link{calc_marginal_value}},
#'   \code{\link{simmar}}, \code{\link{aggregate_yields}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Probe the revenue surface under current effort
#' n_p_a_now <- lapply(fauna, function(cr) cr$n_p_a_0)
#' buffet <- go_fish(e_p_fl, fauna, n_p_a_now, fleets, output_format = "matrix")
#'
#' # Use with allocate_effort to redistribute effort spatially
#' result <- allocate_effort(
#'   effort_by_patch = e_p_fl,
#'   buffet          = buffet,
#'   fleets          = fleets,
#'   open_patch      = open_patches
#' )
#' }
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
      rec_devs = rec_devs,
      move_fish = 0
    )

    yields[[critter]] <- allocate_yields(
      f_p_a_fl = f_p_a_fl, e_p_fl = e_p_fl, p_p_a_fl = p_p_a_fl,
      critter = critter, pop = trial_run, fauna = fauna,
      fleets = fleets, patches = patches, ages = ages
    )
  } ## close fauni loop

  aggregate_yields(yields, e_p_fl, output_format = output_format, groupers = groupers)
}

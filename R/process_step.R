#' Parse a Step Name into Year, Season, and Decimal Year
#'
#' Converts a step name of the form \code{"step_YEAR_SEASON"} or
#' \code{"YEAR_SEASON"} (as returned by \code{\link{simmar}}) into a named
#' list containing the clean year-season string, the integer year, the integer
#' season, and the decimal year used in \code{\link{process_marlin}} output.
#'
#' @param step Character. Step name, e.g. \code{"step_5_3"} or \code{"5_3"}.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{year_season}}{Character. Clean \code{"year_season"} string
#'     (any \code{"step_"} prefix removed).}
#'   \item{\code{year}}{Integer. The year component.}
#'   \item{\code{season}}{Integer. The season component.}
#' }
#'
#' @seealso \code{\link{clean_steps}}, \code{\link{simmar}},
#'   \code{\link{process_marlin}}
#'
#' @export
process_step <- function(step) {
  year_season <- marlin::clean_steps(step)

  year <-
    as.integer(stringr::str_remove_all(year_season, "_.*$"))

  season <-
    as.integer(stringr::str_remove_all(year_season, "^.*_"))

  season_length <- 1 / max(season)

  time_step <- year + (season / max(season) - season_length)

  return(list(
    year_season = year_season,
    year = year,
    season = season,
    time_step = time_step
  ))
}

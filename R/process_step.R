#' Process Step
#' Convert step_YEAR_SEASON notation into YEAR_SEASON, YEAR, and SEASON
#'
#' @param step
#'
#' @return a list with year_season, year, and season
#' @export
#'
process_step <- function(step){

  year_season <- marlin::clean_steps(step)

  year <-
    as.integer(stringr::str_remove_all(year_season,"_.*$"))

  season <-
    as.integer(stringr::str_remove_all(year_season,"^.*_"))

  season_length <- 1 / max(season)

  time_step <- year + (season / max(season) - season_length)

  return(list(year_season = year_season,
              year = year,
              season = season,
              time_step = time_step))

}

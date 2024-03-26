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
    as.integer(gsub("_.*$", "", year_season))

  season <-
    as.integer(gsub("^.*_", "", year_season))

  return(list(year_season = year_season,
              year = year,
              season = season))

}

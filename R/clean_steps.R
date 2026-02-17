#' Clean Step Names from simmar Output
#'
#' Strips any \code{"step_"} prefix from step names in the output of
#' \code{\link{simmar}}, returning clean \code{"year_season"} format strings
#' (e.g. \code{"5_3"} for year 5, season 3).
#'
#' @param step Character. Step name to clean, e.g. \code{"step_5_3"} or
#'   \code{"5_3"}.
#'
#' @return Character string in \code{"year_season"} format.
#'
#' @seealso \code{\link{simmar}}, \code{\link{process_marlin}}
#'
#' @export
#'
#' @examples
#' clean_steps("step_1_2")  # Returns "1_2"
#' clean_steps("10_4")      # Returns "10_4" (no-op)
clean_steps <- function(step) {
  step <- stringr::str_remove_all(step, "step_")
}

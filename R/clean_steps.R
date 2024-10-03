#' A small function to trim any characters from step names in the output of `marlin::simmar`
#'
#' @param step The name of the step in question
#'
#' @return the step in units of {year}_{season index}
#' @export
#'
#' @examples
#' clean_steps("step_1_2")
#'
clean_steps <- function(step){

  step <- stringr::str_remove_all(step,"step_")

    # as.numeric(gsub("^\\D*","",step))

}

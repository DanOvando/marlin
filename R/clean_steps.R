#' A small function to trim any characters from step names in the output of `marlin::simmar`
#'
#' @param step The name of the step in question
#'
#' @return a numeric value with only the numeric (plus decimal) parts of the step name
#' @export
#'
#' @examples
#' clean_steps("step_1")
#' 
clean_steps <- function(step){
  
  step <- as.numeric(gsub("^\\D*","",step))
  
}
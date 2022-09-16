#' Find MSY Referemce Points
#' 
#' Still early days
#'
#' @param fauna 
#' @param fleets 
#'
#' @return
#' @export
#'
#' @examples
find_msy <- function(log_emult,fauna, fleets, years){
  
  if (length(fleets > 1) | length(fauna) > 1){
    stop("find_msy only works with one fleet and one critter at the moment")
  }
  
  emult <- exp(log_emult)
  
  sim <- simmar(fauna = hyper_critter,
                      fleets = hyper_fleet,
                      years = years)
  
  
  
  
}
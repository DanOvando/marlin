#' tune fleets is an internal optimizer
#'
#' @param fauna 
#' @param fleets 
#' @param steps 
#'
#' @return
#' @export
#'
#' @examples
tune_fleets <- function(fauna, fleets, years = 50){
  
  
  qs <- nlminb(start = rep(0,length(fauna) * length(fleets)), fleet_tuner, fleets = fleets, fauna = fauna, years = years, lower = rep(0,length(fauna) * length(fleets)))

  cc <- 1
  
  for (f in seq_along(fleets)){
    
    for (ff in seq_along(fauna)){
      
      fleets[[f]][[ff]]$catchability <- qs$par[cc]
      
      cc <- cc + 1
    }
    
  }
  
  return(fleets)
}
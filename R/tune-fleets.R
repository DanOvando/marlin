#' tune fleets is an internal optimizer
#'
#' @param fauna
#' @param fleets
#' @param years
#'
#' @return
#' @export
#'
#' @examples
tune_fleets <- function(fauna,
                        fleets,
                        years = 50,
                        tune_type = "explt") {
    
  storage <- simmar(fauna = fauna,
                      fleets = fleets,
                      years = years)
    
    fauni <- names(fauna)
    
    fleeti <- names(fleets)
    
    for (s in fauni) {
      e_p_fl <- storage[[length(storage)]][[1]]$e_p_fl
      
      b_p = rowSums(storage[[length(storage)]][[s]]$b_p_a)
      
      # now, calculate biomass weighed total effort of each fleet
      
      # b_p <- rep(1, length(b_p))
      #
      # b_p[1:3] <- 0
      
      # e_p_fl <- e_p_fl * (b_p / max(b_p)) # calcualte the total effort weighted by biomass of that species in patches.
      # browser()
      
      weights <- b_p / max(b_p)
      
      e_fl <- colSums((e_p_fl * weights)) / sum(weights)  # calculate the total effort weighted by biomass of that species in patches.
      
        # colSums(e_p_fl * (b_p / max(b_p)))  # calcualte the total effort weighted by biomass of that species in patches.
      
      p_explt <-
        purrr::map_dbl(fleets, c(s, "p_explt"))[names(e_p_fl)]
      
      explt_by_fleet <- fauna[[s]]$init_explt * p_explt
      
      catchability <-  log(1 - explt_by_fleet) / -e_fl
      
      for (f in fleeti) {
        fleets[[f]][[s]]$catchability <- catchability[f]
        
      } # close internal fleet loop
      
    } # close fauna loop
    

  
  if (tune_type == "depletion") {

    fleet_fauna <- length(fauni) * length(fleeti)
    
    qs <- vector(mode = "double", length = (fleet_fauna))
    
    cc <- 0
    for (f in fauni){
      
      for (ff in fleeti){
        
        cc <- cc + 1
        
        qs[cc] <- fleets[[ff]][[f]]$catchability
        
      }
      
    }
        

    
    qs <-
      nlminb(
        start = qs,
        fleet_tuner,
        fleets = fleets,
        fauna = fauna,
        years = years,
        lower = rep(0, length(fauna) * length(fleets)),
        upper = rep(.9, length(fauna) * length(fleets))
      )
    
    cc <- 1
    
    for (f in seq_along(fleets)) {
      for (ff in seq_along(fauna)) {
        fleets[[f]][[ff]]$catchability <- qs$par[cc]
        
        cc <- cc + 1
      } # close internal fauna loop
      
    } # close fleet loop
  } # close depletion if
  
  return(fleets)
}
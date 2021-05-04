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
    
  
  fleet_names <- names(fleets)
  
  fauni <- names(fauna)
  
  # normalize p_explt to make sure it sums to 1
  
  for (s in fauni) {
    p_explts <- purrr::map_dbl(fleets,c("metiers",s,"p_explt"))
    
    p_explts <- p_explts / ifelse(sum(p_explts) > 0, sum(p_explts), 1e-6)
    
    for ( f in fleet_names){
      
      fleets[[f]]$metiers[[s]]$p_explt <- as.numeric(p_explts[f])
      
    } # close fleet loop
    
  } # close fauna loop
  
  storage <- simmar(fauna = fauna,
                      fleets = fleets,
                      years = years)
    
    fauni <- names(fauna)
    
    fleeti <- names(fleets)
    
    for (s in fauni) {
      e_p_fl <- storage[[length(storage)]][[1]]$e_p_fl
      
      b_p = rowSums(storage[[length(storage)]][[s]]$b_p_a)
      
      weights <- b_p / max(b_p)
      
      e_fl <- colSums((e_p_fl * weights)) / sum(weights)  # calculate the total effort weighted by biomass of that species in patches.
      
      p_explt <-
        purrr::map_dbl(fleets, c("metiers",s, "p_explt"))[names(e_p_fl)]
      
      explt_by_fleet <- (fauna[[s]]$init_explt)  * p_explt
      
      # catchability <-  log(1 - explt_by_fleet) / -e_fl
      
      catchability <-  explt_by_fleet / e_fl
      
      for (f in fleeti) {
        fleets[[f]]$metiers[[s]]$catchability <- catchability[f]
        
        if (all(fleets[[f]]$metiers[[s]]$spatial_catchability == 0)) {
          # annoying step: if q = 0 from earlier, then this will be a matrix of zeros and can't get updated
          fleets[[f]]$metiers[[s]]$spatial_catchability <-
            rep(1, length(fleets[[f]]$metiers[[s]]$spatial_catchability))
        }
        
        mean_q <- mean(fleets[[f]]$metiers[[s]]$spatial_catchability)
        
        mean_q <- ifelse(mean_q == 0, 1e-9, mean_q)
        
        fleets[[f]]$metiers[[s]]$spatial_catchability <-    (fleets[[f]]$metiers[[s]]$spatial_catchability  / mean_q) * catchability[f]
        
      } # close internal fleet loop
      
    } # close fauna loop
    

  
  if (tune_type == "depletion") {

    fleet_fauna <- length(fauni) * length(fleeti)
    
    qs <- vector(mode = "double", length = (fleet_fauna))
    
    cc <- 0
    for (f in fauni){
      
      for (ff in fleeti){
        
        cc <- cc + 1
        
        qs[cc] <- fleets[[ff]]$metiers[[f]]$catchability
        
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
        upper = rep(20, length(fauna) * length(fleets))
      )
    
    cc <- 1
    
    for (f in seq_along(fleets)) {
      for (ff in seq_along(fauna)) {
        fleets[[f]]$metiers[[ff]]$catchability <- qs$par[cc]
        
     
        if (all(fleets[[f]]$metiers[[ff]]$spatial_catchability == 0)) {
          # annoying step: if q = 0 from earlier, then this will be a matrix of zeros and can't get updated
          fleets[[f]]$metiers[[ff]]$spatial_catchability <-
            rep(1, length(fleets[[f]]$metiers[[ff]]$spatial_catchability))
        }
        
        mean_q <- mean(fleets[[f]]$metiers[[ff]]$spatial_catchability)
        
        mean_q <- ifelse(mean_q == 0, 1e-9, mean_q)
        
        
        fleets[[f]]$metiers[[ff]]$spatial_catchability <-fleets[[f]]$metiers[[ff]]$spatial_catchability  / mean_q * qs$par[cc]
        cc <- cc + 1
      } # close internal fauna loop
      
    } # close fleet loop
  } # close depletion if
  
  return(fleets)
}
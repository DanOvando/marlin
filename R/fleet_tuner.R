#' fleet tuner
#' 
#' finds catchability (q) by fleet such that
#' target fished depletion is achieved
#'
#' @param qs 
#' @param fauna 
#' @param fleets 
#' @param steps 
#'
#' @return
#' @export
#'
fleet_tuner <- function(qs,fauna, fleets, steps = 50){
  
  cc <- 1
  
  for (f in seq_along(fleets)){
    
    for (ff in seq_along(fauna)){
      
      fleets[[f]][[ff]]$catchability <- qs[cc]
      
      cc <- cc + 1
    }
    
  }
  
  storage <- simmar(fauna = fauna,
                    fleets = fleets,
                    steps = steps)
  
  tmp <- purrr::map_dfr(storage[[length(storage)]], ~as.data.frame(.x$ssb_p_a), .id = "fauna")
  
  b_p <- rowSums(tmp[,2:ncol(tmp)], na.rm = TRUE)
  
  tmp <- data.frame(fauna = tmp$fauna, ssb = b_p) %>% 
    dplyr::group_by(fauna) %>% 
    dplyr::summarise(ssb = sum(ssb)) %>% 
    dplyr::arrange(fauna)
  
  ssb0s <- purrr::map_dbl(fauna, "ssb0")
  
  ssb0s <- ssb0s[sort(names(ssb0s))]
  
  target_depletion <- purrr::map_dbl(fauna, "fished_depletion")
  
  target_depletion <- target_depletion[sort(names(target_depletion))]
  
  tmp$depletion <- tmp$ssb / ssb0s
  
  ss <- sum((log(tmp$depletion) - log(target_depletion))^2)
  
  return(ss)
  
  
}
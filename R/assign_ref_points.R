#' Assign References Points
#'
#' @param fauna a list of critters
#' @param fleets a list of cleets
#'
#' @return a fauna object but with MSY based reference points included for each critter
#' @export
#'
assign_ref_points <- function(fauna, fleets){
  
  for (f in seq_along(fauna)){
    
    msy_mult <- optim(1e-3,find_msy, lower = 0, upper = 10, fauna = fauna, fleets = fleets, opt = TRUE,target_critter = names(fauna)[f])
    
    msy_state <- find_msy(msy_mult$par, fauna = fauna, fleets = fleets, opt = FALSE, target_critter =  names(fauna)[f])
    
    ref_points <- msy_state$fauna %>% 
      dplyr::filter(critter == names(fauna)[f]) %>% 
      dplyr::summarise(ssb_msy = sum(ssb), 
                       b_msy  = sum(b),
                       n_msy = sum(n),
                       msy = sum(c),
                       u_msy = sum(c) / sum(b))
    
    base_e_msy <- mean(purrr::map_dbl(fleets,"base_effort")) * msy_mult$par
    
    ref_points$base_e_msy_mult <- msy_mult$par
    
    ref_points$base_e_msy <- base_e_msy
    
    fauna[[f]]$ref_points <- ref_points
    
  }
  
  return(fauna)
}

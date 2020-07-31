#' Create Fleet
#' 
#' Creates a fleet object, mostly by adding in 
#' selectivity at age for each fleet and species
#'
#' @param fleets 
#' @param fauna 
#' @param base_effort 
#'
#' @return
#' @export
#'
create_fleet <-
  function(fleets,
           fauna,
           base_effort = 0) {
    # idea: each fleet has a list of fauna inside of it specifying the price, selectivity, q for that species
    
    fleet_names <- names(fleets)
    
    fauni <- names(fauna)
    
    for (f in seq_along(fleet_names)) {
      tmp_fleet <- fleets[[f]]
      
      for (s in seq_along(fauni)) {
        tmp_critter <- tmp_fleet[[fauni[s]]]
        
        length_bins <-
          as.numeric(colnames(fauna[[fauni[s]]]$length_at_age_key))
        
        l_50_sel <-
          fauna[[fauni[s]]]$length_50_mature * tmp_critter$sel_start
        
        l_95_sel <-
          fauna[[fauni[s]]]$length_50_mature * (tmp_critter$sel_start + tmp_critter$sel_delta)
        
        if (tmp_critter$sel_form == "logistic") {
        
          sel_at_bin <- ((1 / (1 + exp(
            -log(19) * ((length_bins - l_50_sel) / (l_95_sel - l_50_sel))
          ))))
          
          p_sel_at_age <-
            as.matrix(fauna[[fauni[s]]]$length_at_age_key) %*% sel_at_bin
          
          sel_at_age <- p_sel_at_age
          
          fleets[[f]][[fauni[s]]]$sel_at_age <- as.numeric(sel_at_age)
          
        } else if (tmp_critter$sel_form == "dome"){ # close logistic form if
        
          sel_at_bin <- dnorm(length_bins, l_50_sel, sd = (l_95_sel - l_50_sel))
          
          p_sel_at_age <-
            (as.matrix(fauna[[fauni[s]]]$length_at_age_key) %*% sel_at_bin)
          
          sel_at_age <- p_sel_at_age / max(p_sel_at_age)
          
          fleets[[f]][[fauni[s]]]$sel_at_age <- as.numeric(sel_at_age)
          
        } # close dome shaped
        
        
      } # close fauni loop
      
    
      fleets[[f]]$base_effort <-  base_effort
      
      
    } # close fleet loop
    
    
    return(fleets)
    
  } # close function
#' R6 Class Representing a fishing fleet
#'
#' @description
#' A fleet object has all the required characteristics of a fishing fleet
#'
#' @details
#' creates fleet object with spaces for selectivity, mpa response, etc. 

Metier <- R6::R6Class("metier",
                     lock_objects = FALSE,
                public = list(
                  initialize = function(critter = NA, # this might be redundant
                                        price = 10,
                                        sel_form = "logistic",
                                        sel_start = 1, 
                                        sel_delta = .1,
                                        catchability = 0.01,
                                        spatial_catchability = NA,
                                        p_explt = 1) {
                    
                      self$price <- price
                      
                      self$sel_form <-  sel_form
                      
                      self$sel_start <-  sel_start
                      
                      self$sel_delta <-  sel_delta
                      
                      self$catchability <- catchability
                    
                      self$p_explt <- p_explt
                      
                      length_bins <-
                        as.numeric(colnames(critter$length_at_age_key))
                      
                      l_50_sel <-
                        critter$length_50_mature * sel_start
                      
                      l_95_sel <-
                        critter$length_50_mature * (sel_start + sel_delta)
                      
                      if (sel_form == "logistic") {
                        
                        sel_at_bin <- ((1 / (1 + exp(
                          -log(19) * ((length_bins - l_50_sel) / (l_95_sel - l_50_sel))
                        ))))
                        
                        p_sel_at_age <-
                          as.matrix(critter$length_at_age_key) %*% sel_at_bin
                        
                        sel_at_age <- p_sel_at_age
                        
                        self$sel_at_age <- as.numeric(sel_at_age)
                        
                      } else if (sel_form == "dome"){ # close logistic form if
                        
                        sel_at_bin <- dnorm(length_bins, l_50_sel, sd = (l_95_sel - l_50_sel))
                        
                        p_sel_at_age <-
                          (as.matrix(critter$length_at_age_key) %*% sel_at_bin)
                        
                        sel_at_age <- p_sel_at_age / max(p_sel_at_age)
                        
                        self$sel_at_age <- as.numeric(sel_at_age)
                        
                      } # close dome shaped
                      
                      if (all(is.na(spatial_catchability))){
                        
                        self$spatial_catchability <- rep(1, critter$patches)
                        
                        
                      } else {
                        if (unique(dim(spatial_catchability)) != sqrt(critter$patches)) {
                          stop(
                            glue::glue(
                              "spatial_catchability must either be NA or a {sqrt(critter$patches}) by {sqrt(critter$patches)} matrix"
                            )
                          )
                        } #close dim check
                        
                        tmp <- spatial_catchability %>%
                          as.data.frame() %>%
                          dplyr::mutate(x = 1:nrow(.)) %>%
                          tidyr::pivot_longer(
                            -x,
                            names_to = "y",
                            values_to = "catchability",
                            names_prefix = "V"
                          ) %>%
                          dplyr::mutate(catchability = catchability / mean(catchability))
                        
                        self$spatial_catchability <- tmp$catchability * catchability
                        
                      } # close deal with spatial q
                  } # close initialize
                ) # close public
) # close class
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
                  #' @description
                  #' Parameters for specific metier of a given fleet
                  #'
                  #' @param critter the name of the critter in the fauna object this metier applies to
                  #' @param price the price per unit weight of the critter in question
                  #' @param sel_form the selectivity form, one of "logistic", "dome", or "manual"
                  #' @param sel_unit the unit of selectivity, one of "p_of_mat" which means selectivity is in proportion of age at maturity, or "length" where selectivity is in units of length
                  #' @param sel_start the value of sel_unit at which selectivity "starts"
                  #' @param sel_delta the delta parameter in the selectivity function
                  #' @param catchability the catchability per uni effort pararmeter, generally overwritten by tune_fleet
                  #' @param spatial_catchability a matrix of spatial q
                  #' @param sel_at_age a manual vector of gear (contact) selectivity at age, where values are between 0 and 1
                  #' @param p_explt the proportion of total exploitation for a given critter coming from this metier. This value is relaive to all other p_explt values for the critter in question. Set to 0 to have metier not catch critter at all
                  initialize = function(critter = NA, # this might be redundant
                                        price = 10,
                                        sel_form = "logistic",
                                        sel_unit = "p_of_mat",
                                        sel_start = 1,
                                        sel_delta = .1,
                                        catchability = 0.2,
                                        spatial_catchability = NA,
                                        p_explt = 1,
                                        sel_at_age = NULL) {


                      catchability <- pmax(1e-9,catchability)

                      self$price <- price

                      self$sel_form <-  sel_form

                      self$sel_start <-  sel_start

                      self$sel_delta <-  sel_delta

                      self$catchability <- catchability

                      self$p_explt <- p_explt

                      self$port_distance <- NA

                      length_bins <-
                        as.numeric(colnames(critter$length_at_age_key))

                      if (sel_unit == "p_of_mat"){
                      l_50_sel <-
                        critter$length_50_mature * sel_start

                      l_95_sel <-
                        critter$length_50_mature * (sel_start + sel_delta)

                      } else if (sel_unit == "length"){

                        l_50_sel <- sel_start

                        l_95_sel <- sel_start + sel_delta

                      }

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

                      } else if (sel_form == "manual") {

                        if (is.null(sel_at_age)){
                          stop("sel_form = 'manual' but no manual set_at_age provided")
                        }

                        self$sel_at_age <- sel_at_age

                      } # close sel_form things



                      if (all(is.na(spatial_catchability))){

                        self$spatial_catchability <- rep(catchability, critter$patches)


                      } else {
                        if (unique(dim(spatial_catchability)) != sqrt(critter$patches)) {
                          stop(
                            glue::glue(
                              "spatial_catchability must either be NA or a {sqrt(critter$patches}) by {sqrt(critter$patches)} matrix"
                            )
                          )
                        } #close dim check


                        if (any(spatial_catchability < 0)){
                          spatial_catchability <- spatial_catchability - min(spatial_catchability)

                        }

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

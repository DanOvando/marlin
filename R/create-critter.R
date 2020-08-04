#' Create Critter
#' 
#' Creates a critter object. If only a scientific name is provided
#' create_critter will try and look up relevant life history from
#' FishLife
#' 
#' Critical inputs are adult_movement, adult_movement_sigma, and resolution
#' 
#'

create_critter <- function(common_name = 'white seabass',
                           scientific_name = NA,
                           critter_type = "fish",
                           seasonal_habitat = list(),
                           habitat_seasons = list(seasons),
                           rec_habitat = NA,
                           seasons = 1,
                           rec_form = 1,
                           adult_movement = 0, 
                           adult_movement_sigma = 2,
                           fished_depletion = 0.4,
                           ...) {

  
  if (critter_type == "fish"){
    
    critter <-
      marlin::Fish$new(
        common_name = common_name,
        scientific_name = scientific_name,
        seasonal_habitat = seasonal_habitat,
        habitat_seasons = habitat_seasons,
        rec_habitat = rec_habitat,
        seasons = seasons,
        rec_form = rec_form,
        fished_depletion = fished_depletion,
        adult_movement = adult_movement,
        adult_movement_sigma = adult_movement_sigma
      )
    
  }

  return(critter)
}

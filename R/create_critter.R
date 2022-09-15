#' Create Critter
#' 
#' Creates a critter object. If only a scientific name is provided
#' create_critter will try and look up relevant life history from
#' FishLife
#' 
#' Critical inputs are adult_movement, adult_movement_sigma, and resolution
#' 
#'
#' @param common_name 
#' @param scientific_name 
#' @param critter_type 
#' @param base_habitat 
#' @param season_blocks 
#' @param rec_habitat 
#' @param seasons 
#' @param rec_form 
#' @param adult_movement 
#' @param adult_movement_sigma 
#' @param fished_depletion 
#' @param init_explt initial annual exploitation rate
#' @param explt_type f or fmsy
#' @param ... 

create_critter <- function(common_name = NA,
                           scientific_name = NA,
                           get_common_name = FALSE,
                           critter_type = "fish",
                           base_habitat = list(),
                           season_blocks = list(),
                           recruit_habitat = NA,
                           seasons = 1,
                           rec_form = 1,
                           fec_form = "weight",
                           adult_diffusion = 2,
                           recruit_diffusion = 10,
                           fished_depletion = 0.4,
                           init_explt = .1,
                           explt_type = "f",
                           burn_years = 50,
                           weight_a = NA,
                           fec_expo = 1,
                           resolution = 20,
                           spawning_seasons = NA,
                           ...) {

  
  if (!is.list(base_habitat)){
    base_habitat <-  list(base_habitat)
  }

  init_explt <- init_explt / seasons # convert to seasonal exploitation rate
  if (critter_type == "fish"){
    critter <-
      marlin::Fish$new(
        common_name = common_name,
        scientific_name = scientific_name,
        base_habitat = base_habitat,
        season_blocks = season_blocks,
        recruit_habitat = recruit_habitat,
        seasons = seasons,
        rec_form = rec_form,
        fec_form = fec_form,
        fec_expo = fec_expo,
        weight_a = weight_a,
        fished_depletion = fished_depletion,
        adult_diffusion = adult_diffusion,
        recruit_diffusion = recruit_diffusion,
        init_explt = init_explt,
        explt_type = explt_type,
        burn_years = burn_years,
        get_common_name = get_common_name,
        resolution = resolution,
        spawning_seasons = spawning_seasons,
        ...
      )
    
  }

  return(critter)
}
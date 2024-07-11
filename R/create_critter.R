#' Create Critter
#'
#' Creates a critter object. If only a scientific name is provided
#' create_critter will try and look up relevant life history from
#' FishLife
#'
#' Critical inputs are adult_movement, adult_movement_sigma, and resolution
#'
#' See ?Fish for documentation of all possible parameters
#'
#'
#' @param common_name the common name of the species
#' @param scientific_name the scientific name of the species, preferable to common name
#' @param habitat a list with adult habitat per season
#' @param season_blocks a list with seasons per block
#' @param seasons number of seasons per year (integer)
#' @param fished_depletion depletion (biomass / unfished biomass) at start of simulation
#' @param init_explt initial annual exploitation rate (fraction of exploitable population killed)
#' @param explt_type f or fmsy
#' @param recruit_habitat habitat for recruitment
#' @param fec_form one of "weight" or "pups"
#' @param adult_diffusion adult diffusion rate
#' @param recruit_diffusion recruit diffusion rate
#' @param burn_years number of years to burn in the simulation prior to starting things
#' @param weight_a alpha parameter in the allometric weight function alpha x length ^ beta
#' @param fec_expo exponent for fecundity relationship. >1 means hyperallometric fecundity
#' @param resolution the resolution of the system, either an integer or a vector integers of length two specifying the dimensions of the system in width and height (e.g. `c(10,100)`)
#' @param patch_area the area of each patch
#' @param spawning_seasons which seasons spawning occurs
#' @param density_dependence one of 'global_habitat','local_habitat','pre_dispersal','post_dispersal','global_ssb'
#' @param get_common_name TRUE or FALSE to lookup common name from scientific name. Requires internet connection
#' @param critter_type placeholder for someday if non-Fish objects are implemented
#' @param ... additional parameters passed to `Fish` class, see `?Fish`

create_critter <- function(common_name = NA,
                           scientific_name = NA,
                           get_common_name = FALSE,
                           critter_type = "fish",
                           habitat = list(),
                           season_blocks = list(),
                           recruit_habitat = NA,
                           seasons = 1,
                           fec_form = "weight",
                           adult_diffusion = 2,
                           recruit_diffusion = 10,
                           fished_depletion = 0.4,
                           init_explt = .1,
                           explt_type = "f",
                           burn_years = 50,
                           weight_a = NA,
                           fec_expo = 1,
                           resolution = c(10,10),
                           patch_area = 1,
                           spawning_seasons = NA,
                           density_dependence = "global_habitat",
                           ...) {


  if (!is.list(habitat)){
    habitat <-  list(habitat)
  }

  if (critter_type == "fish"){
    critter <-
      marlin::Fish$new(
        common_name = common_name,
        scientific_name = scientific_name,
        habitat = habitat,
        season_blocks = season_blocks,
        recruit_habitat = recruit_habitat,
        seasons = seasons,
        density_dependence = density_dependence,
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

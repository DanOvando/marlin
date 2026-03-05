#' Create a Critter (Species) Object
#'
#' @description
#' Constructs a critter (species population) object for use in
#' \code{\link{simmar}}. Life history parameters (growth, maturity, natural
#' mortality, weight-at-age) are looked up from FishLife when
#' \code{scientific_name} is supplied, or can be passed manually via \code{...}.
#' The returned object is an R6 \code{Fish} instance initialised to unfished
#' equilibrium and burned in for \code{burn_years} years.
#'
#' @details
#' ## Spatial setup
#' \code{resolution} determines the 2-D grid: \code{c(nx, ny)} produces
#' \code{nx * ny} patches indexed by \code{expand_grid(x, y)} (x varies
#' slowest). \code{habitat} should be a list (one matrix per season block) of
#' \code{[ny, nx]} habitat-quality matrices; these drive adult movement taxis
#' and, optionally, recruit distribution.
#'
#' ## Movement
#' Adult and recruit movement is modelled via diffusion-taxis. If
#' \code{adult_home_range} / \code{recruit_home_range} are provided (km), they
#' are used to derive diffusion coefficients automatically and take precedence
#' over \code{adult_diffusion} / \code{recruit_diffusion}. If a home range is
#' \code{NULL}, the corresponding diffusion value is used directly.
#'
#' ## Density dependence options
#' \describe{
#'   \item{\code{"global_habitat"}}{Beverton-Holt spawner-recruit using global
#'     SSB; recruits distributed proportionally to habitat.}
#'   \item{\code{"local_habitat"}}{Recruitment density dependence acts locally
#'     within each patch.}
#'   \item{\code{"pre_dispersal"}}{Competition occurs before larval dispersal.}
#'   \item{\code{"post_dispersal"}}{Competition occurs after larval dispersal.}
#'   \item{\code{"global_ssb"}}{Beverton-Holt using global SSB with uniform
#'     spatial recruit distribution.}
#' }
#'
#' @param common_name Character. Common name of the species (e.g.
#'   \code{"yellowfin tuna"}). Used to look up life history from FishLife when
#'   \code{scientific_name} is not supplied.
#' @param scientific_name Character. Scientific name (e.g.
#'   \code{"Thunnus albacares"}). Preferred over \code{common_name} for
#'   FishLife lookup; case-insensitive.
#' @param get_common_name Logical. If \code{TRUE}, resolves the common name
#'   from \code{scientific_name} via an internet lookup. Default \code{FALSE}.
#' @param critter_type Character. Currently only \code{"fish"} is supported.
#'   Placeholder for future non-fish implementations.
#' @param habitat List of habitat-quality matrices (one per season block).
#'   Each matrix must be \code{[ny, nx]}. Values represent relative habitat
#'   quality on a log scale; they drive movement taxis and, depending on
#'   \code{density_dependence}, recruit distribution. Generate with
#'   \code{\link{sim_habitat}}.
#' @param season_blocks List of integer vectors specifying which seasons belong
#'   to each habitat block (e.g. \code{list(1:2, 3:4)} for a two-block,
#'   four-season year). Must align with \code{habitat}.
#' @param seasons Integer. Number of seasons per year. Sets
#'   \code{time_step = 1 / seasons}. All critters in a \code{fauna} list must
#'   share the same value.
#' @param depletion Numeric in (0, 1]. Target spawning stock biomass relative
#'   to unfished SSB (B/B0) at the start of the simulation.
#' @param fished_depletion Numeric. Alias for \code{depletion}; takes
#'   precedence when both are supplied.
#' @param init_explt Numeric. Initial annual exploitation rate (fraction of the
#'   exploitable population removed per year).
#' @param f Numeric. Alias for \code{init_explt}.
#' @param explt_type Character. One of \code{"f"} (instantaneous fishing
#'   mortality; default) or \code{"fmsy"} (as a multiple of Fmsy).
#' @param recruit_habitat Matrix or \code{NA}. Habitat matrix used to
#'   distribute new recruits. Defaults to adult habitat of the first season
#'   block when \code{NA}.
#' @param fec_form Character. Fecundity form: \code{"weight"} (fecundity
#'   proportional to body weight; default) or \code{"pups"} (fixed litter
#'   size, e.g. for elasmobranchs).
#' @param lorenzen_c Numeric. Exponent in the Lorenzen (1996) size-dependent
#'   natural mortality function. Negative values (default \code{-1}) give
#'   higher mortality for smaller fish.
#' @param adult_home_range Numeric or \code{NULL}. Characteristic adult home range
#'   (km). When not \code{NULL}, this is used to derive the adult diffusion
#'   coefficient and takes precedence over \code{adult_diffusion}. Set to
#'   \code{NULL} to use \code{adult_diffusion} directly.
#' @param recruit_home_range Numeric or \code{NULL}. Characteristic recruit home
#'   range (km). When not \code{NULL}, this is used to derive the recruit
#'   diffusion coefficient and takes precedence over \code{recruit_diffusion}.
#'   Set to \code{NULL} to use \code{recruit_diffusion} directly.
#' @param adult_diffusion Numeric or \code{NULL}. Adult diffusion rate
#'   (km^2 per time step). Used only when \code{adult_home_range} is \code{NULL}.
#' @param recruit_diffusion Numeric or \code{NULL}. Recruit diffusion rate.
#'   Used only when \code{recruit_home_range} is \code{NULL}.
#' @param burn_years Integer. Number of years to burn in the unfished
#'   population to equilibrium before starting the simulation. Default 50.
#' @param weight_a Numeric or \code{NA}. Intercept alpha in the allometric
#'   weight-at-length relationship \eqn{W = \alpha L^\beta}. Looked up from
#'   FishLife when \code{NA}.
#' @param fec_expo Numeric. Exponent for the fecundity-weight relationship.
#'   Values > 1 produce hyperallometric fecundity. Default \code{1}.
#' @param resolution Integer scalar or length-2 integer vector \code{c(nx, ny)}
#'   giving grid dimensions. A scalar is replicated to a square grid.
#' @param patch_area Numeric. Area of each patch (km^2). Used to scale
#'   diffusion and compute density-dependent quantities.
#' @param spawning_seasons Integer vector. Which seasons spawning occurs in.
#'   Defaults to all seasons when \code{NA}.
#' @param density_dependence Character. Density dependence form; see Details.
#'   One of \code{"global_habitat"}, \code{"local_habitat"},
#'   \code{"pre_dispersal"}, \code{"post_dispersal"}, \code{"global_ssb"}.
#' @param steepness Numeric in (0.2, 1). Beverton-Holt steepness parameter
#'   \eqn{h} governing the shape of the spawner-recruit relationship.
#' @param growth_model Character. Growth model for length-at-age:
#'   \code{"von_bertalanffy"} (default), \code{"power"}, or
#'   \code{"growth_cessation"}.
#' @param ... Additional parameters forwarded to the \code{\link{Fish}} R6
#'   class constructor. Common options include \code{k}, \code{linf},
#'   \code{t0}, \code{m}, \code{age_mature}, \code{weight_b}, \code{ssb0},
#'   \code{sigma_rec} (recruitment standard deviation), and \code{ac_rec}
#'   (recruitment autocorrelation). See \code{?Fish} for the full list.
#'
#' @return An R6 \code{Fish} object containing the age-structured population
#'   state at unfished equilibrium, movement matrices, life-history schedules,
#'   and reference points. Pass inside a named list as \code{fauna} to
#'   \code{\link{simmar}}.
#'
#' @seealso \code{\link{simmar}}, \code{\link{create_fleet}},
#'   \code{\link{tune_fleets}}, \code{\link{sim_habitat}}, \code{\link{Fish}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Minimal: look up life history from FishLife
#' tuna <- create_critter(
#'   scientific_name = "Thunnus albacares",
#'   resolution      = c(10, 10),
#'   seasons         = 4
#' )
#'
#' # With explicit spatial structure
#' snapper <- create_critter(
#'   scientific_name    = "Lutjanus campechanus",
#'   adult_home_range   = 10,
#'   recruit_home_range = 20,
#'   resolution         = c(10, 10),
#'   seasons            = 1,
#'   density_dependence = "global_habitat",
#'   fished_depletion   = 0.4
#' )
#'
#' fauna <- list(snapper = snapper)
#' }
create_critter <- function(common_name = NA,
                           scientific_name = NA,
                           get_common_name = FALSE,
                           critter_type = "fish",
                           habitat = list(),
                           season_blocks = list(),
                           recruit_habitat = NA,
                           seasons = 1,
                           lorenzen_c = -1,
                           fec_form = "weight",
                           adult_home_range = 1,
                           recruit_home_range = 2,
                           adult_diffusion = NULL,
                           recruit_diffusion = NULL,
                           depletion = 0.4,
                           fished_depletion = NULL,
                           init_explt = .1,
                           f = NULL,
                           explt_type = "f",
                           burn_years = 50,
                           weight_a = NA,
                           fec_expo = 1,
                           resolution = c(10, 10),
                           patch_area = 1,
                           spawning_seasons = NA,
                           density_dependence = "global_habitat",
                           steepness = 0.8,
                           growth_model = "von_bertalanffy",
                           ...) {
  if (!is.list(habitat)) {
    habitat <- list(habitat)
  }
  
  # Create checks for land (NAs) in habitat layers 
  ## Find NAs in adult habitat layers - should be the same for each 
  ## item in the list so we can just use the first item
  habitat_df <- habitat[[1]] |>
    as.data.frame() |>
    dplyr::mutate(y = dplyr::n():1) |>
    tidyr::pivot_longer(-y, values_to = "value") |>
    dplyr::group_by(y) |>
    dplyr::mutate(x = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::arrange(x, y)
  
  habitat_NAs <- which(is.na(habitat_df$value))
  
  ## Find NAs in recruit habitat layer
  recruit_habitat_df <- recruit_habitat |>
    as.data.frame() |>
    dplyr::mutate(y = dplyr::n():1) |>
    tidyr::pivot_longer(-y, values_to = "value") |>
    dplyr::group_by(y) |>
    dplyr::mutate(x = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::arrange(x, y)
  
  recruit_habitat_NAs <- which(is.na(recruit_habitat_df$value))
  
  ## If one has NAs and one doesn't (due to forgetfullness) - 
  ## add the NAs and add a warning
  if(length(recruit_habitat_NAs) == 0 & length(habitat_NAs) > 0) { 
    # Add NAs to the data.frame
    recruit_habitat_df$value[habitat_NAs] <- NA
    
    # Re-arrange for proper matricing
    recruit_habitat_df <- recruit_habitat_df %>% 
      dplyr::arrange(x, desc(y))
    
    # Overwrite original habitat using the new one
    recruit_habitat <- matrix(recruit_habitat_df$value, 
                              nrow = length(unique(recruit_habitat_df$y)), 
                              ncol = length(unique(recruit_habitat_df$x)))
    
    warning("Land areas (NAs) are present in the adult habitat layer, but not the recruit habitat layer. Adding land areas to recruit habitat layer...")
    recruit_habitat_NAs <- habitat_NAs
  }
  
  if(length(recruit_habitat_NAs) > 0 & length(habitat_NAs) == 0) { 
    # This one is a little trickier since it can be a list with many layers
    habitat <- purrr::map(.x = 1:length(habitat), 
                          .f = ~{
                            temp_df <- habitat |>
                              as.data.frame() |>
                              dplyr::mutate(y = dplyr::n():1) |>
                              tidyr::pivot_longer(-y, values_to = "value") |>
                              dplyr::group_by(y) |>
                              dplyr::mutate(x = dplyr::row_number()) |>
                              dplyr::ungroup() |>
                              dplyr::arrange(x, y)
                            
                            # Add NAs to the data.frame
                            temp_df$value[recruit_habitat_NAs] <- NA
                            
                            # Re-arrange for proper matricing
                            temp_df <- temp_df %>% 
                              dplyr::arrange(x, desc(y))
                            
                            # Overwrite original habitat using the new one
                           matrix(temp_df$value,
                                  nrow = length(unique(temp_df$y)), 
                                  ncol = length(unique(temp_df$x)))
                            })
    warning("Land areas (NAs) are present in the recruit habitat layer, but not the adult habitat layer. Adding land areas to adult habitat layer...")
    habitat_NAs <- recruit_habitat_NAs
  }
  
  ## If the locations of NAs are different, throw an error
  if(length(habitat_NAs) > 0 & length(recruit_habitat_NAs) > 0 & (!identical(sort(habitat_NAs), sort(recruit_habitat_NAs)))) { 
    stop("Land areas (NAs) must be identical in the supplied adult and recruit habitat layers")
  }

  
  if (!is.null(f)){
    init_explt = f
  }
  if (!is.null(fished_depletion)){
    depletion = fished_depletion
  }

  # Home range takes precedence over explicit diffusion rates
  if (!is.null(adult_home_range)) {
    adult_diffusion <- NULL
  }
  if (!is.null(recruit_home_range)) {
    recruit_diffusion <- NULL
  }

  if (critter_type == "fish") {
    critter <-
      marlin::Fish$new(
        common_name = common_name,
        scientific_name = scientific_name,
        habitat = habitat,
        season_blocks = season_blocks,
        recruit_habitat = recruit_habitat,
        seasons = seasons,
        density_dependence = density_dependence,
        lorenzen_c = lorenzen_c,
        fec_form = fec_form,
        fec_expo = fec_expo,
        weight_a = weight_a,
        depletion = depletion,
        adult_home_range = adult_home_range,
        recruit_home_range = recruit_home_range,
        adult_diffusion = adult_diffusion,
        recruit_diffusion = recruit_diffusion,
        init_explt = init_explt,
        explt_type = explt_type,
        burn_years = burn_years,
        get_common_name = get_common_name,
        resolution = resolution,
        spawning_seasons = spawning_seasons,
        growth_model = growth_model,
        patch_area = patch_area,
        steepness = steepness,
        ...
      )
  }

  return(critter)
}

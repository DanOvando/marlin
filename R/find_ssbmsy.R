#' Find Baseline SSBmsy
#' 
#' And assigns it to each critter in fauna object
#'
#' @param mult multiplier to effort
#' @param sel_start selectivity as a multiplier of length at maturity
#' @param fauna a fauna object
#' @param years the number of years to run things out for
#' @param use one of "graphs" or "optim"
#'
#' @return results of simulation
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' years = 50
#' 
#' seasons = 4
#' 
#' fauna <- 
#'list(
#'  "bigeye" = create_critter(
#'    scientific_name =  "thunnus obesus",
#'    adult_diffusion = 10,
#'    density_dependence = "post_dispersal",
#'    seasons = seasons,
#'    resolution = resolution,
#'    age_mature = 1,
#'    steepness = 0.9,
#'    ssb0 = 1000
#'  )
#')
#' 
#' baseline_ssbmsy <- find_ssbmsy(fauna = fauna)
#' 
#' }
#' 

assign_ssbmsy <- function(fauna, sel_start = 0.01, years = 50){
  
  inner_find_ssbmsy <- function(mult = 1,sel_start = .01,fauni, years = 50, use = "graphs"){
    
    
    tmp <- list(fauni)
    
    names(tmp) <- 'a'
    
    tmp_metier <- list(Metier$new(
      critter = tmp[[1]],
      sel_form = "logistic",
      sel_start = sel_start,
      sel_delta = .01,
      p_explt = 1,
      catchability = 1
    ))
    
    names(tmp_metier) <- names(tmp)
    
    fleets <- list(
      "longline" = create_fleet(
        tmp_metier,
        base_effort = mult * prod(resolution),
        resolution = resolution
      )
    )
    
    
    sim <- simmar(fauna = tmp, fleets = fleets, years = years)
    
    eqish <- sim[[length(sim)]]
    if (use == "graphs"){
      out <- data.frame(ssb = sum(eqish[[1]]$ssb_p_a),
                        b = sum(eqish[[1]]$b_p_a),
                        yield = sum(eqish[[1]]$c_p_a))
      
    } else if (use == "optim") {
      out <- -sum(eqish[[1]]$c_p_a)
    }
    
  }
  
  
  inner_inner_find_ssbmsy <- function(fauni,  sel_start, years){
    emsy <- optim(
      0.1,
      inner_find_ssbmsy,
      fauni = fauni,
      lower = 0,
      upper = 1,
      use = "optim",
      sel_start = sel_start,
      years = years,
      method = "L-BFGS-B"
    )
    
    baseline_ssbmsy <- inner_find_ssbmsy(emsy$par, fauni = fauni,sel_start = sel_start, years = years)$ssb
    
  }
  
  ssbmsys <- purrr::map(fauna, inner_inner_find_ssbmsy, sel_start = sel_start, years = years)
  
  for (f in names(fauna)){
    
    fauna[[f]]$baseline_ssbmsy <- ssbmsys[[f]]
    
  }
  return(fauna)
}



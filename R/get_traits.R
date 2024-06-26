#' Get Traits
#'
#' Retrieves life history traits from FishLife
#'
#' This function returns the mean un-logged life history traits for the closest match to the
#' supplied taxonomic information.
#'
#' @param Class Character input for taxonomic class
#' @param Order Character input for taxonomic class
#' @param Family Character input for taxonomic class
#' @param Genus Character input for taxonomic class
#' @param Species Character input for taxonomic class
#' @param verbose logical where TRUE prints closest match, FALSE does not
#'
#' @return a dataframe of mean trait values
#' @export
#'
#' @examples
#' \dontrun{
#' life_traits <- Get_traits(Genus = "Lutjanus", Species = "campechanus")
#' }
get_traits <-
  function(Class = "predictive",
           Order = "predictive",
           Family = "predictive",
           Genus = "predictive",
           Species = "predictive") {
    closest_match <-
      marlin::search_species(
        Class = Class,
        Order = Order,
        Family = Family,
        Genus = Genus,
        Species = Species
      )
    
    closest_taxa_match <- closest_match$closest_match
    
    trait_table <-
      as.data.frame(t(marlin::FishBase_and_RAM$ParHat$beta_gj[closest_match$GroupNum[[1]], ]))
    
    trait_table[colnames(trait_table) != 'Temperature'] <-
      exp(trait_table[colnames(trait_table) != 'Temperature'])
    
    trait_table$closest_taxa_match <- closest_taxa_match
    
    return(trait_table)
  }

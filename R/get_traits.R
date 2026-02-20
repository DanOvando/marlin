#' Get Life History Traits from FishLife
#'
#' @description
#' Retrieves mean (un-logged) life history traits for the closest taxonomic
#' match in the FishLife database. Traits are looked up by matching the
#' supplied taxonomic hierarchy against the \code{FishBase_and_RAM}
#' database bundled with marlin.
#'
#' @details
#' Traits are returned on the natural scale (exponentiated from the log scale
#' used internally by FishLife), except \code{Temperature} which is returned
#' as-is. Useful for inspecting the life history that \code{\link{create_critter}}
#' will use when \code{scientific_name} is supplied, or for populating custom
#' parameter sets for less well-known species.
#'
#' @param Class Character. Taxonomic class (e.g. \code{"Actinopterygii"}).
#'   Default \code{"predictive"} (let FishLife infer from lower ranks).
#' @param Order Character. Taxonomic order. Default \code{"predictive"}.
#' @param Family Character. Taxonomic family. Default \code{"predictive"}.
#' @param Genus Character. Genus name (e.g. \code{"Lutjanus"}). Default
#'   \code{"predictive"}.
#' @param Species Character. Species epithet (e.g. \code{"campechanus"}).
#'   Default \code{"predictive"}. Supply at least \code{Genus} and
#'   \code{Species} for species-level lookup.
#'
#' @return A data frame with one row containing mean trait values for the
#'   closest taxonomic match, including columns for growth parameters
#'   (\code{K}, \code{Linf}, \code{t0}), natural mortality (\code{M}),
#'   maturity (\code{tm}), weight-length (\code{a}, \code{b}),
#'   \code{Temperature}, and \code{closest_taxa_match} (the matched taxon
#'   string from FishLife).
#'
#' @seealso \code{\link{create_critter}}, \code{\link{search_species}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Look up red snapper traits
#' traits <- get_traits(Genus = "Lutjanus", Species = "campechanus")
#' traits$K       # von Bertalanffy growth rate
#' traits$Linf    # asymptotic length (cm)
#' traits$M       # natural mortality rate
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

    trait_table[colnames(trait_table) != "Temperature"] <-
      exp(trait_table[colnames(trait_table) != "Temperature"])

    trait_table$closest_taxa_match <- closest_taxa_match

    return(trait_table)
  }

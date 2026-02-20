#' Search FishLife Taxonomy for the Closest Species Match
#'
#' @description
#' Matches taxonomic inputs against the \code{ParentChild_gz} taxonomy tree
#' in the FishLife database, returning the row number(s) of the closest
#' ancestor match. Ported and adapted from the FishLife package by James
#' Thorson.
#'
#' @details
#' The function queries WORMS (World Register of Marine Species) via
#' \code{taxize} to resolve the full taxonomy for the supplied genus/species,
#' then sweeps from Order down to Species looking for the closest match in
#' FishLife's \code{ParentChild_gz}. Unspecified levels are padded with
#' \code{"predictive"} and the function falls back to the nearest ancestor
#' if an exact species match is not found.
#'
#' This function requires an internet connection. It is called internally
#' by \code{\link{get_traits}} and \code{\link{create_critter}}.
#'
#' @param Class Character. Taxonomic class. Default \code{"predictive"}.
#' @param Order Character. Taxonomic order. Default \code{"predictive"}.
#' @param Family Character. Taxonomic family. Default \code{"predictive"}.
#' @param Genus Character. Genus name. Default \code{"predictive"}.
#' @param Species Character. Species epithet. Default \code{"predictive"}.
#' @param add_ancestors Logical. If \code{TRUE} (default), returns row
#'   indices for all ancestor taxa as well as the matched taxon. Used by
#'   FishLife's hierarchical prediction.
#' @param Database List. The FishLife database to search; defaults to
#'   \code{marlin::FishBase_and_RAM}.
#' @param ParentChild_gz Data frame. The parent-child taxonomy table from
#'   \code{Database}; extracted automatically.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{GroupNum}}{Integer vector of row indices in
#'     \code{ParentChild_gz} for the matched taxon and its ancestors.}
#'   \item{\code{match_taxonomy}}{Character vector of matched taxon strings
#'     (in \code{"Genus_Species_predictive_..."} format).}
#'   \item{\code{closest_match}}{Character. The best-matching taxon string.}
#' }
#'
#' @references
#' Thorson, J.T. (FishLife R package).
#' \url{https://github.com/James-Thorson-NOAA/FishLife}
#'
#' @seealso \code{\link{get_traits}}, \code{\link{create_critter}}
#'
#' @export
search_species <- function(Class = "predictive",
                           Order = "predictive",
                           Family = "predictive",
                           Genus = "predictive",
                           Species = "predictive",
                           add_ancestors = TRUE,
                           Database = marlin::FishBase_and_RAM,
                           ParentChild_gz = Database$ParentChild_gz) {
  quiet_query <- purrr::quietly(taxize::get_wormsid)
  worms_id <-
    (quiet_query(paste0(
      Genus, ifelse(Species == "predictive", "", paste0(" ", Species))
    ), accepted = TRUE, ask = FALSE, rows = 1)$result)

  taxonomy <- taxize::classification(worms_id, db = "worms")[[1]]

  taxonomy$rank <- tolower(taxonomy$rank)

  # add missing taxonomic levels from taxize
  full_taxonomy <- c(Class, Order, Family, Genus, Species)
  if (!all(c(Species) == "predictive")) {
    full_taxonomy[5] <- taxonomy$name[taxonomy$rank == "species"]
  }
  if (!all(c(Species, Genus) == "predictive")) {
    full_taxonomy[4] <- taxonomy$name[taxonomy$rank == "genus"]
  }
  if (!all(c(Species, Genus, Family) == "predictive")) {
    full_taxonomy[3] <- taxonomy$name[taxonomy$rank == "family"]
  }
  if (!all(c(Species, Genus, Family, Order) == "predictive")) {
    full_taxonomy[2] <- taxonomy$name[taxonomy$rank == "order"]
  }
  if (!all(c(Species, Genus, Family, Order, Class) == "predictive")) {
    # superclass seems to more frequently match the values in FishLife so defaulting to this for now
    class <- taxonomy$name[taxonomy$rank == "superclass"]

    if (length(class) == 0) {
      class <- taxonomy$name[taxonomy$rank == "class"]
    }

    full_taxonomy[1] <- class
  }


  # check if species has been returned as Genus + species
  if (length(strsplit(full_taxonomy[[5]], "\\s+")[[1]]) > 1) {
    full_taxonomy[[5]] <-
      strsplit(full_taxonomy[[5]], "\\s+")[[1]][2] # If "species" entry contains Genus + species, convert to the second element, assuming this is species, splitting by whitespace
  }

  match_taxonomy <- full_taxonomy

  # find closest matches in FishLife database, based on entries with the most number of taxonomic matches.
  # Doing this instead of exact matches since taxonomy in FishLife no longer matches taxonomy in Fishbase / taxize, e.g.
  # FishLife says Thunnus opesus is order Perciformes, whereas resources now classify them as in order Scombriformes

  fishlife_match_counter <-
    purrr::map_dbl(ParentChild_gz$ChildName,
      \(x, match_taxonomy) sum(match_taxonomy %in% strsplit(as.character(x), "_")[[1]]),
      match_taxonomy = match_taxonomy
    )

  Group <- which.max(fishlife_match_counter)

  closest_match <- as.character(ParentChild_gz[Group, "ChildName"])

  # Pick out ancestors
  if (add_ancestors == TRUE) {
    Group <- marlin::find_ancestors(child_num = Group, ParentChild_gz = ParentChild_gz)
  }

  # Function to add predictive to taxon name
  Add_predictive <- function(char_vec) {
    return_vec <- char_vec
    for (i in 1:length(return_vec)) {
      vec <- strsplit(as.character(return_vec[i]), "_")[[1]]
      return_vec[i] <- paste(c(vec, rep("predictive", 5 - length(vec))),
        collapse =
          "_"
      )
    }
    return(return_vec)
  }
  match_taxonomy <- unique(as.character(Add_predictive(ParentChild_gz[Group, "ChildName"])))
  # Find new matches
  GroupNum <- match(match_taxonomy, ParentChild_gz[, "ChildName"])

  # Return match
  Return <- list("GroupNum" = GroupNum, "match_taxonomy" = match_taxonomy, closest_match = closest_match)
  return(Return)
}

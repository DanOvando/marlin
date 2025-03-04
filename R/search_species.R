#' Search species: Ported from https://github.com/James-Thorson-NOAA/FishLife/blob/main/R/Search_species.R
#'
#' Match taxonomic inputs to a given row of \code{ParentChild_gz} or its closest ancestor
#'
#' This function attempts to do a smart match to elements of \code{ParentChild_gz}.  It sweeps from Order to Species
#' and ignores any taxonomic input listed as \code{"predictive"} until it finds something else.  It then appends
#' \code{"predictive"} to any lower taxonomic level that is missing, and checks whether this specification yields a single,
#' unique taxon.  If it does, it then returns the row number and potentially any ancestors (higher taxonomic levels)
#'
#' @param Class Character input for taxonomic class
#' @param Order Character input for taxonomic class
#' @param Family Character input for taxonomic class
#' @param Genus Character input for taxonomic class
#' @param Species Character input for taxonomic class
#' @param add_ancestors Boolean whether to add ancestors for matching species or not
#' @param ParentChild_gz vector providing index of parent-taxon for every child-taxa
#'
#' @return integer of row numbers of \code{ParentChild_gz} matching \code{genus_species}

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

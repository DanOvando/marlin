#' Find ancestors Ported fromhttps://github.com/James-Thorson-NOAA/FishLife/blob/main/R/Find_ancestors.R
#'
#' Find higher taxonomic levels for a given taxon (e.g., Class and Order for a given Family)
#'
#' @param child_num row number of \code{ParentChild_gz} for which to find ancestors

#' @return vector of row numbers of \code{ParentChild_gz} for ancestors (including \code{child_num})

#' @export
find_ancestors = function( child_num, Database=marlin::FishBase_and_RAM, ParentChild_gz=Database$ParentChild_gz ){
  
  # Search for all ancestors in the taxonomic tree
  family_nums = child_num
  while(TRUE){
    if( is.na(ParentChild_gz[rev(family_nums)[1],'ParentRowNumber'])==TRUE ) break()
    family_nums = c(family_nums, ParentChild_gz[rev(family_nums)[1],'ParentRowNumber'])
  }
  
  # Return vector of ancestors
  return( family_nums )
}
# Search FishLife Taxonomy for the Closest Species Match

Matches taxonomic inputs against the `ParentChild_gz` taxonomy tree in
the FishLife database, returning the row number(s) of the closest
ancestor match. Ported and adapted from the FishLife package by James
Thorson.

## Usage

``` r
search_species(
  Class = "predictive",
  Order = "predictive",
  Family = "predictive",
  Genus = "predictive",
  Species = "predictive",
  add_ancestors = TRUE,
  Database = marlin::FishBase_and_RAM,
  ParentChild_gz = Database$ParentChild_gz
)
```

## Arguments

- Class:

  Character. Taxonomic class. Default `"predictive"`.

- Order:

  Character. Taxonomic order. Default `"predictive"`.

- Family:

  Character. Taxonomic family. Default `"predictive"`.

- Genus:

  Character. Genus name. Default `"predictive"`.

- Species:

  Character. Species epithet. Default `"predictive"`.

- add_ancestors:

  Logical. If `TRUE` (default), returns row indices for all ancestor
  taxa as well as the matched taxon. Used by FishLife's hierarchical
  prediction.

- Database:

  List. The FishLife database to search; defaults to
  [`marlin::FishBase_and_RAM`](https://danovando.github.io/marlin/reference/FishBase_and_RAM.md).

- ParentChild_gz:

  Data frame. The parent-child taxonomy table from `Database`; extracted
  automatically.

## Value

A named list with three elements:

- `GroupNum`:

  Integer vector of row indices in `ParentChild_gz` for the matched
  taxon and its ancestors.

- `match_taxonomy`:

  Character vector of matched taxon strings (in
  `"Genus_Species_predictive_..."` format).

- `closest_match`:

  Character. The best-matching taxon string.

## Details

The function queries WORMS (World Register of Marine Species) via
`taxize` to resolve the full taxonomy for the supplied genus/species,
then sweeps from Order down to Species looking for the closest match in
FishLife's `ParentChild_gz`. Unspecified levels are padded with
`"predictive"` and the function falls back to the nearest ancestor if an
exact species match is not found.

This function requires an internet connection. It is called internally
by
[`get_traits`](https://danovando.github.io/marlin/reference/get_traits.md)
and
[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md).

## References

Thorson, J.T. (FishLife R package).
<https://github.com/James-Thorson-NOAA/FishLife>

## See also

[`get_traits`](https://danovando.github.io/marlin/reference/get_traits.md),
[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md)

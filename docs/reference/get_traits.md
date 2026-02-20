# Get Life History Traits from FishLife

Retrieves mean (un-logged) life history traits for the closest taxonomic
match in the FishLife database. Traits are looked up by matching the
supplied taxonomic hierarchy against the `FishBase_and_RAM` database
bundled with marlin.

## Usage

``` r
get_traits(
  Class = "predictive",
  Order = "predictive",
  Family = "predictive",
  Genus = "predictive",
  Species = "predictive"
)
```

## Arguments

- Class:

  Character. Taxonomic class (e.g. `"Actinopterygii"`). Default
  `"predictive"` (let FishLife infer from lower ranks).

- Order:

  Character. Taxonomic order. Default `"predictive"`.

- Family:

  Character. Taxonomic family. Default `"predictive"`.

- Genus:

  Character. Genus name (e.g. `"Lutjanus"`). Default `"predictive"`.

- Species:

  Character. Species epithet (e.g. `"campechanus"`). Default
  `"predictive"`. Supply at least `Genus` and `Species` for
  species-level lookup.

## Value

A data frame with one row containing mean trait values for the closest
taxonomic match, including columns for growth parameters (`K`, `Linf`,
`t0`), natural mortality (`M`), maturity (`tm`), weight-length (`a`,
`b`), `Temperature`, and `closest_taxa_match` (the matched taxon string
from FishLife).

## Details

Traits are returned on the natural scale (exponentiated from the log
scale used internally by FishLife), except `Temperature` which is
returned as-is. Useful for inspecting the life history that
[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md)
will use when `scientific_name` is supplied, or for populating custom
parameter sets for less well-known species.

## See also

[`create_critter`](https://danovando.github.io/marlin/reference/create_critter.md),
[`search_species`](https://danovando.github.io/marlin/reference/search_species.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Look up red snapper traits
traits <- get_traits(Genus = "Lutjanus", Species = "campechanus")
traits$K       # von Bertalanffy growth rate
traits$Linf    # asymptotic length (cm)
traits$M       # natural mortality rate
} # }
```

# Find ancestors Ported fromhttps://github.com/James-Thorson-NOAA/FishLife/blob/main/R/Find_ancestors.R

Find higher taxonomic levels for a given taxon (e.g., Class and Order
for a given Family)

## Usage

``` r
find_ancestors(
  child_num,
  Database = marlin::FishBase_and_RAM,
  ParentChild_gz = Database$ParentChild_gz
)
```

## Arguments

- child_num:

  row number of `ParentChild_gz` for which to find ancestors

## Value

vector of row numbers of `ParentChild_gz` for ancestors (including
`child_num`)

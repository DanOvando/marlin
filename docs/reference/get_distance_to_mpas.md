# Measures the distance of each cell to the nearest MPA edge and to all MPA cells

Measures the distance of each cell to the nearest MPA edge and to all
MPA cells

## Usage

``` r
get_distance_to_mpas(mpa_locations, resolution, patch_area = 10)
```

## Arguments

- mpa_locations:

  a dataframe with at least column coordinate columns x,y,and mpa (TRUE
  or FALSE)

- resolution:

  the resolution of the simulation system. Can supply a vector `c(x,y)`
  to denote an X by Y system, or one number to denote an X by X system

- patch_area:

  the area of each patch

## Value

a dataframe with two columns added. `distance_to_mpa_edge` measures the
distance from each patch to the nearest MPA edge in units of
sqrt(patch_area), with negative values indicating areas inside an MPA
`total_mpa_distance` measure to total distance to all MPA cells from
each patch, in units of sqrt(patch_area)

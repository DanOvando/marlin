# marlin Colour Palettes

Returns a `colorRampPalette` function for one of seven marlin-themed
colour palettes, all derived from the colours of a blue marlin (*Makaira
nigricans*).

## Usage

``` r
marlin_pal(palette = "fish_scales", reverse = FALSE, ...)
```

## Arguments

- palette:

  Character. Name of the palette; see Details. Default `"fish_scales"`.

- reverse:

  Logical. If `TRUE`, reverses the palette colour order. Default
  `FALSE`.

- ...:

  Additional arguments passed to `colorRampPalette`.

## Value

A `colorRampPalette` function that accepts an integer `n` and returns
`n` hex colour strings.

## Details

Available palettes:

- `"fish_scales"`:

  11-colour sequential palette spanning the full dorsal-to-ventral
  colour range of the fish.

- `"diverging_fish"`:

  7-colour diverging palette suitable for difference maps or comparisons
  between two treatments.

- `"lateral_lines"`:

  8-colour palette inspired by the lateral line markings.

- `"dark_blues"`:

  4-colour palette of deep blues from the dorsal and fin regions.

- `"sea_blues"`:

  5-colour palette of mid-range blues.

- `"sky_blues"`:

  4-colour palette of lighter blues and grey-blues from the ventral
  region.

- `"sands"`:

  4-colour palette of warm sandy neutrals.

## See also

[`theme_marlin`](https://danovando.github.io/marlin/reference/theme_marlin.md),
[`plot_marlin`](https://danovando.github.io/marlin/reference/plot_marlin.md)

## Examples

``` r
# 5-colour diverging palette
cols <- marlin_pal("diverging_fish")(5)
scales::show_col(cols)


# Use in a ggplot
library(ggplot2)
ggplot(mtcars, aes(mpg, wt, colour = factor(cyl))) +
  geom_point() +
  scale_colour_manual(values = marlin_pal("sea_blues")(3))
```

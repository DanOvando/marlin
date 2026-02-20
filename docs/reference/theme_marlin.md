# ggplot2 Theme for marlin

A clean `ggplot2` theme based on `theme_classic` with marlin styling:
bordered panel, light dashed gridlines, bold-italic titles, dark blue
facet strip backgrounds, and italicised legend text.

## Usage

``` r
theme_marlin(base_size = 14, ...)
```

## Arguments

- base_size:

  Numeric. Base font size in points. Default `14`.

- ...:

  Additional arguments passed to
  [`ggplot2::theme_classic`](https://ggplot2.tidyverse.org/reference/ggtheme.html).

## Value

A `ggplot2` `theme` object. Add to a ggplot with `+`.

## See also

[`plot_marlin`](https://danovando.github.io/marlin/reference/plot_marlin.md),
[`marlin_pal`](https://danovando.github.io/marlin/reference/marlin_pal.md)

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(mpg, wt)) +
  geom_point() +
  theme_marlin()
```

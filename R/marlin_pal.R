#' marlin Colour Palettes
#'
#' @description
#' Returns a \code{colorRampPalette} function for one of seven marlin-themed
#' colour palettes, all derived from the colours of a blue marlin
#' (\emph{Makaira nigricans}).
#'
#' @details
#' Available palettes:
#' \describe{
#'   \item{\code{"fish_scales"}}{11-colour sequential palette spanning the
#'     full dorsal-to-ventral colour range of the fish.}
#'   \item{\code{"diverging_fish"}}{7-colour diverging palette suitable for
#'     difference maps or comparisons between two treatments.}
#'   \item{\code{"lateral_lines"}}{8-colour palette inspired by the lateral
#'     line markings.}
#'   \item{\code{"dark_blues"}}{4-colour palette of deep blues from the dorsal
#'     and fin regions.}
#'   \item{\code{"sea_blues"}}{5-colour palette of mid-range blues.}
#'   \item{\code{"sky_blues"}}{4-colour palette of lighter blues and
#'     grey-blues from the ventral region.}
#'   \item{\code{"sands"}}{4-colour palette of warm sandy neutrals.}
#' }
#'
#' @param palette Character. Name of the palette; see Details.
#'   Default \code{"fish_scales"}.
#' @param reverse Logical. If \code{TRUE}, reverses the palette colour order.
#'   Default \code{FALSE}.
#' @param ... Additional arguments passed to \code{colorRampPalette}.
#'
#' @return A \code{colorRampPalette} function that accepts an integer \code{n}
#'   and returns \code{n} hex colour strings.
#'
#' @seealso \code{\link{theme_marlin}}, \code{\link{plot_marlin}}
#'
#' @export
#'
#' @examples
#' # 5-colour diverging palette
#' cols <- marlin_pal("diverging_fish")(5)
#' scales::show_col(cols)
#'
#' # Use in a ggplot
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt, colour = factor(cyl))) +
#'   geom_point() +
#'   scale_colour_manual(values = marlin_pal("sea_blues")(3))
marlin_pal <- function(palette = "fish_scales", reverse = FALSE, ...) {
  # specify colors
  marlin_colors <- c(
    "fin_tip" = "#13191E",
    "fin_mid" = "#2F385C",
    "fin_base" = "#515680",
    "dorsal_tip" = "#224E5F",
    "dorsal_mid" = "#286073",
    "dorsal_light" = "#5DACC7",
    "ventral1" = "#3E749E",
    "ventral2" = "#8A9CB8",
    "ventral_light" = "#C7D0DD",
    "lateral_line_dark" = "#6A6344",
    "lateral_line" = "#C8C3AA",
    "lateral_line_light" = "#D7D4C2",
    "lateral_line_highlight1" = "#839F6D",
    "lateral_line_lowlight1" = "#4C7462",
    "lateral_line_highlight2" = "#EDDEA4",
    "lateral_line_lowlight2" = "#E0C764",
    "lateral_line_highlight3" = "#E7E7E7",
    "lateral_line_lowlight3" = "#A2A2A2"
  )

  marlin_cols <- function(...) {
    cols <- c(...)

    if (is.null(cols)) {
      return(marlin_colors)
    }

    marlin_colors[cols]
  }

  # specify palettes
  marlin_palettes <- list(
    "fish_scales" = marlin_cols("fin_tip", "fin_mid", "fin_base", "dorsal_tip", "dorsal_mid", "ventral1", "ventral2", "lateral_line_highlight1", "lateral_line", "lateral_line_highlight2", "lateral_line_highlight3"),
    "diverging_fish" = marlin_cols("fin_tip", "fin_base", "dorsal_light", "ventral1", "lateral_line_highlight1", "lateral_line_lowlight2", "lateral_line_lowlight3"),
    "lateral_lines" = marlin_cols("lateral_line_lowlight1", "lateral_line_highlight1", "lateral_line_lowlight2", "lateral_line_highlight2", "lateral_line", "lateral_line_dark", "lateral_line_lowlight3", "lateral_line_highlight3"),
    "dark_blues" = marlin_cols("fin_tip", "fin_base", "ventral2", "lateral_line_highlight3"),
    "sea_blues" = marlin_cols("dorsal_tip", "dorsal_mid", "dorsal_light", "lateral_line_lowlight1", "lateral_line_highlight1"),
    "sky_blues" = marlin_cols("ventral1", "ventral2", "ventral_light", "lateral_line_highlight3"),
    "sands" = marlin_cols("lateral_line_dark", "lateral_line", "lateral_line_light", "lateral_line_highlight3")
  )

  pal <- marlin_palettes[[palette]]

  if (reverse) pal <- rev(pal)

  return(colorRampPalette(pal, ...))
}

#'  marlin pal palette
#'
#' @param palette
#' @param reverse
#' @param ...
#'
#' @return color palletes
#' @export
#'
#' @examples
#' image(1:11, 9:10, as.matrix(1:11),
#'   col = marlin_pal("fish_scales")(11),
#'   xlim = c(-4.5, 12), ylim = c(0, 10),
#'   xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n"
#' )
#' text(-1.8, 9.5, 'marlin_pal("fish_scales")', cex = 0.7)
#' image(1:7, 7.5:8.5, as.matrix(1:7),
#'   col = marlin_pal("diverging_fish")(7),
#'   xlim = c(-4.5, 12), ylim = c(0, 10),
#'   xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", add = TRUE
#' )
#' text(-1.8, 8, 'marlin_pal("diverging_fish")', cex = 0.7)
#' image(1:11, 6:7, as.matrix(1:11),
#'   col = marlin_pal("lateral_lines")(11),
#'   xlim = c(-4.5, 12), ylim = c(0, 10),
#'   xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", add = TRUE
#' )
#' text(-1.8, 6.5, 'marlin_pal("lateral_lines")', cex = 0.7)
#' image(1:11, 4.5:5.5, as.matrix(1:11),
#'   col = marlin_pal("dark_blues")(11),
#'   xlim = c(-4.5, 9), ylim = c(0, 10),
#'   xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", add = TRUE
#' )
#' text(-1.8, 5, 'marlin_pal("dark_blues")', cex = 0.7)
#' image(1:11, 3:4, as.matrix(1:11),
#'   col = marlin_pal("sea_blues")(11),
#'   xlim = c(-4.5, 9), ylim = c(0, 10),
#'   xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", add = TRUE
#' )
#' text(-1.8, 3.5, 'marlin_pal("sea_blues")', cex = 0.7)
#' image(1:11, 1.5:2.5, as.matrix(1:11),
#'   col = marlin_pal("sky_blues")(11),
#'   xlim = c(-4.5, 9), ylim = c(0, 10),
#'   xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", add = TRUE
#' )
#' text(-1.8, 2, 'marlin_pal("sky_blues")', cex = 0.7)
#' image(1:11, 0:1, as.matrix(1:11),
#'   col = marlin_pal("sands")(11),
#'   xlim = c(-4.5, 9), ylim = c(0, 10),
#'   xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", add = TRUE
#' )
#' text(-1.8, 0.5, 'marlin_pal("sands")', cex = 0.7)
#'
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

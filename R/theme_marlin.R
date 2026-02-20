#' ggplot2 Theme for marlin
#'
#' @description
#' A clean \code{ggplot2} theme based on \code{theme_classic} with marlin
#' styling: bordered panel, light dashed gridlines, bold-italic titles, dark
#' blue facet strip backgrounds, and italicised legend text.
#'
#' @param base_size Numeric. Base font size in points. Default \code{14}.
#' @param ... Additional arguments passed to \code{ggplot2::theme_classic}.
#'
#' @return A \code{ggplot2} \code{theme} object. Add to a ggplot with \code{+}.
#'
#' @seealso \code{\link{plot_marlin}}, \code{\link{marlin_pal}}
#'
#' @export
#' @importFrom ggplot2 %+replace%
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   theme_marlin()
theme_marlin <- function(base_size = 14, ...) {
  ggplot2::theme_classic(...) %+replace%
    ggplot2::theme(
      # background
      panel.background = ggplot2::element_blank(),
      # borders and axis lines
      panel.border = ggplot2::element_rect(linewidth = 1, fill = "NA"),
      axis.line = ggplot2::element_blank(),
      # gridlines
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.1,
        linetype = 2,
        colour = "lightgray"
      ),
      panel.grid.minor = ggplot2::element_line(
        linewidth = 0.1,
        linetype = 2,
        colour = "lightgray"
      ),
      # title
      plot.title = ggplot2::element_text(
        face = "bold.italic",
        size = 14,
        margin = ggplot2::margin(t = 5, b = 5),
        vjust = 0.5,
        hjust = 0
      ),
      plot.title.position = "panel",
      # strips for faceted plots
      strip.background = ggplot2::element_rect(fill = "#2F385C"),
      strip.text = ggplot2::element_text(
        colour = "white",
        face = "italic",
        margin = ggplot2::margin(t = 2, b = 2),
        vjust = 0.5
      ),
      # legends
      legend.title = ggplot2::element_text(face = "italic", hjust = 0.5),
      # text
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = 0.75 * base_size)
    )
}

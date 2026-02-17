#' Aggregate Per-Species Yield Outputs Into Fleet-Level Totals
#'
#' @description
#' Sums catch, revenue, and profit across species from the outputs of
#' \code{\link{allocate_yields}} to produce a fleet-level "buffet" of spatial
#' yield metrics. Also computes per-unit-effort (PUE) counterparts.
#' This is the final step that builds the buffet consumed by
#' \code{\link{allocate_effort}}.
#'
#' @details
#' Used in two contexts:
#' \enumerate{
#'   \item Inside \code{\link{go_fish}}, to aggregate exploratory fishing
#'     results across species before returning the buffet.
#'   \item Inside the \code{\link{simmar}} step loop, after
#'     \code{\link{allocate_yields}} has been called for each species.
#' }
#'
#' @param yields Named list of per-species yield outputs (one element per
#'   species). Each element must be the output of \code{\link{allocate_yields}}
#'   and contain patches x fleets matrices: \code{r_p_fl}, \code{c_p_fl},
#'   and \code{prof_p_fl}. Names should match species names in \code{fauna}.
#' @param e_p_fl Numeric matrix of effort by patch (rows) and fleet (columns),
#'   with fleet names as column names. Used to compute per-unit-effort metrics.
#' @param output_format Character. Output format:
#'   \describe{
#'     \item{\code{"matrix"}}{(default) Returns list of patches x fleets
#'       matrices with fleet names as column names. Fast; designed for
#'       internal use in \code{\link{simmar}}.}
#'     \item{\code{"tidy"}}{Returns list of tidy data frames with species
#'       summed and effort joined. Slower; used inside \code{\link{go_fish}}
#'       for user-facing output.}
#'   }
#' @param groupers Character vector. Grouping columns for tidy format
#'   (default: \code{c("fleet", "patch")}). Ignored when
#'   \code{output_format = "matrix"}.
#'
#' @return A named list with six elements (summed across all species):
#' \describe{
#'   \item{\code{r_p_fl}}{Revenue by patch and fleet}
#'   \item{\code{c_p_fl}}{Catch by patch and fleet}
#'   \item{\code{prof_p_fl}}{Profit by patch and fleet}
#'   \item{\code{rpue_p_fl}}{Revenue per unit effort (\code{NA} where effort = 0)}
#'   \item{\code{cpue_p_fl}}{Catch per unit effort (\code{NA} where effort = 0)}
#'   \item{\code{ppue_p_fl}}{Profit per unit effort (\code{NA} where effort = 0)}
#' }
#'
#' @seealso \code{\link{allocate_yields}}, \code{\link{go_fish}},
#'   \code{\link{simmar}}, \code{\link{allocate_effort}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Inside simmar after the species loop:
#' buffet <- aggregate_yields(yields_this_step, updated_e_p_f)
#'
#' # Inside go_fish with tidy output:
#' buffet <- aggregate_yields(yields, e_p_fl, output_format = "tidy")
#' }
aggregate_yields <- function(
    yields,
    e_p_fl,
    output_format = c("matrix", "tidy"),
    groupers = c("fleet", "patch")
) {

  output_format <- match.arg(output_format)

  fauni <- names(yields)
  fleet_names <- colnames(e_p_fl)
  n_fleets <- ncol(e_p_fl)
  n_patches <- nrow(e_p_fl)

  nms <- c("r_p_fl", "c_p_fl", "prof_p_fl")
  pue_nms <- c("rpue_p_fl", "cpue_p_fl", "ppue_p_fl")

  if (output_format == "matrix") {

    # Sum each matrix across species
    raw <- setNames(
      lapply(nms, function(mat_name) {
        mat_out <- matrix(0, nrow = n_patches, ncol = n_fleets,
                          dimnames = list(NULL, fleet_names))
        for (sp in fauni) {
          mat_out <- mat_out + yields[[sp]][[mat_name]]
        }
        mat_out
      }),
      nms
    )

    # Per-unit-effort: NA where effort is zero
    pue <- setNames(
      lapply(nms, function(mat_name) {
        out <- raw[[mat_name]] / e_p_fl
        out[e_p_fl == 0] <- NA
        out
      }),
      pue_nms
    )

    c(raw, pue)

  } else {

    # Tidy format
    e_p_fl_df <-
      tibble::as_tibble(e_p_fl) |>
      dplyr::mutate(patch = dplyr::row_number(), .before = 1) |>
      tidyr::pivot_longer(
        cols = -patch,
        names_to = "fleet",
        values_to = "effort"
      )

    raw <- setNames(
      purrr::map(nms, function(mat_name) {

        purrr::imap_dfr(yields, function(species_element, species_name) {

          mat <- species_element[[mat_name]]

          tibble::as_tibble(mat) |>
            dplyr::mutate(patch = dplyr::row_number(), .before = 1) |>
            tidyr::pivot_longer(
              cols = -patch,
              names_to = "fleet",
              values_to = "value"
            ) |>
            dplyr::mutate(species = species_name, .before = 1)

        }) |>
          dplyr::group_by(dplyr::across(dplyr::all_of(groupers))) |>
          dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop") |>
          dplyr::left_join(e_p_fl_df, by = c("patch", "fleet"))

      }),
      nms
    )

    # Per-unit-effort: NA where effort is zero
    pue <- setNames(
      lapply(nms, function(mat_name) {
        raw[[mat_name]] |>
          dplyr::mutate(value = dplyr::if_else(effort == 0, NA_real_, value / effort))
      }),
      pue_nms
    )

    c(raw, pue)
  }
}

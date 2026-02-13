#' Aggregate Yields Across Species
#'
#' @description
#' Takes per-species yield outputs (from \code{\link{allocate_yields}}) and sums
#' catch, revenue, and profit across species to produce fleet-level totals.
#' Also computes per-unit-effort versions of each. This is the step that
#' builds the \code{buffet} that \code{\link{allocate_effort}} selects from.
#'
#' Can be used in two contexts:
#' \itemize{
#'   \item Inside \code{\link{go_fish}} to aggregate exploratory fishing results
#'   \item Inside the \code{\link{simmar}} step loop after \code{allocate_yields}
#'     has been called for each species
#' }
#'
#' @param yields Named list of per-species yield outputs. Each element should be
#'   the output of \code{\link{allocate_yields}} and must contain matrices
#'   \code{r_p_fl}, \code{c_p_fl}, and \code{prof_p_fl} (patches x fleets).
#'   Names should be species names.
#' @param e_p_fl Numeric matrix of effort by patch (rows) and fleet (columns).
#'   Used to compute per-unit-effort metrics. Column names should be fleet names.
#' @param output_format Character: \code{"matrix"} (default) returns list of
#'   patches x fleets matrices; \code{"tidy"} returns list of tidy data frames
#'   with species summed and effort joined.
#' @param groupers Character vector of grouping columns for tidy output
#'   (default: \code{c("fleet", "patch")}). Ignored when
#'   \code{output_format = "matrix"}.
#'
#' @return A named list with six elements:
#' \describe{
#'   \item{r_p_fl}{Revenue by patch and fleet (summed across species)}
#'   \item{c_p_fl}{Catch by patch and fleet (summed across species)}
#'   \item{prof_p_fl}{Profit by patch and fleet (summed across species)}
#'   \item{rpue_p_fl}{Revenue per unit effort by patch and fleet (NA where effort = 0)}
#'   \item{cpue_p_fl}{Catch per unit effort by patch and fleet (NA where effort = 0)}
#'   \item{ppue_p_fl}{Profit per unit effort by patch and fleet (NA where effort = 0)}
#' }
#'
#' Format of each element depends on \code{output_format}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Inside simmar, after the fauni loop:
#' buffet <- aggregate_yields(yields, e_p_fl)
#'
#' # Inside go_fish:
#' buffet <- aggregate_yields(yields, e_p_fl, output_format = "tidy")
#' }
#'
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

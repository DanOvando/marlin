#' Calculate Marginal Revenue and Marginal Profit by Patch and Fleet
#'
#' @description
#' Uses forward finite differences on \code{\link{go_fish}} to approximate
#' \eqn{\partial R / \partial E_{p,f}} (marginal revenue) and
#' \eqn{\partial \Pi / \partial E_{p,f}} (marginal profit) for every
#' patch-fleet combination. The output matrices are designed to be merged
#' directly into the \code{buffet} list consumed by
#' \code{\link{allocate_effort}} when a fleet's \code{spatial_allocation}
#' is \code{"marginal_revenue"} or \code{"marginal_profit"}, or when
#' \code{fleet_model = "sole_owner"}.
#'
#' @details
#' Two approximation strategies are available:
#'
#' \describe{
#'   \item{\code{"separable"}}{Bumps effort in **all patches** of one fleet
#'     simultaneously by \code{epsilon}, then computes derivatives in a single
#'     extra \code{go_fish} call per fleet. Requires \eqn{n_{fleets}} extra
#'     evaluations. Fast, but assumes the marginal return in patch \eqn{i}
#'     does not depend on the perturbation in patch \eqn{j} (valid when
#'     per-patch depletion effects are negligible, which is typical for large
#'     grids).}
#'   \item{\code{"patch_loop"}}{Bumps effort in **one patch at a time**,
#'     requiring \eqn{n_{patches} \times n_{fleets}} extra evaluations. Fully
#'     accounts for cross-patch depletion effects; appropriate for small,
#'     dense grids where within-step movements are significant.}
#' }
#'
#' Called automatically by \code{\link{simmar}} when any fleet uses marginal
#' allocation or the sole-owner model. The results are appended to the buffet
#' as \code{mr_p_fl} and \code{mp_p_fl}.
#'
#' @param e_p_fl Numeric matrix of effort by patch (rows) and fleet (columns).
#' @param fauna Named list of fauna objects (as in \code{\link{go_fish}}).
#' @param n_p_a Named list of numbers-at-age matrices by patch (one element per
#'   species), representing the current population state.
#' @param fleets Named list of fleet objects (as in \code{\link{go_fish}}).
#' @param baseline Optional. List returned by a prior call to
#'   \code{\link{go_fish}}\code{(output_format = "matrix")} with the same
#'   inputs. Supplying this avoids one redundant baseline evaluation and
#'   can improve performance when marginals are computed alongside an
#'   existing buffet.
#' @param method Character. Approximation strategy; see Details.
#'   One of \code{"separable"} (default) or \code{"patch_loop"}.
#' @param epsilon Numeric. Additive effort perturbation for the finite
#'   difference. Should be small relative to typical effort values but large
#'   enough to avoid floating-point noise. Default \code{1e-3}.
#'
#' @return A named list with two elements, each a patches x fleets matrix
#'   with the same dimension names as \code{e_p_fl}:
#' \describe{
#'   \item{\code{mr_p_fl}}{Marginal revenue
#'     \eqn{\partial R / \partial E_{p,f}} by patch and fleet.}
#'   \item{\code{mp_p_fl}}{Marginal profit
#'     \eqn{\partial \Pi / \partial E_{p,f}} by patch and fleet.}
#' }
#' Append to the buffet from \code{\link{go_fish}} as
#' \code{buffet$mr_p_fl} and \code{buffet$mp_p_fl} for use by
#' \code{\link{allocate_effort}}.
#'
#' @seealso \code{\link{go_fish}}, \code{\link{allocate_effort}},
#'   \code{\link{simmar}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Compute baseline buffet
#' buffet <- go_fish(e_p_fl, fauna, n_p_a, fleets, output_format = "matrix")
#'
#' # Append marginal values (reuse baseline to avoid redundant evaluation)
#' marginals <- calc_marginal_value(
#'   e_p_fl   = e_p_fl,
#'   fauna    = fauna,
#'   n_p_a    = n_p_a,
#'   fleets   = fleets,
#'   baseline = buffet,
#'   method   = "separable"
#' )
#'
#' buffet$mr_p_fl <- marginals$mr_p_fl
#' buffet$mp_p_fl <- marginals$mp_p_fl
#'
#' # Now allocate effort using marginal profit
#' result <- allocate_effort(
#'   effort_by_patch = e_p_fl,
#'   buffet          = buffet,
#'   fleets          = fleets,
#'   open_patch      = open_patches
#' )
#' }
calc_marginal_value <- function(e_p_fl,
                                fauna,
                                n_p_a,
                                fleets,
                                baseline = NULL,
                                method = c("separable", "patch_loop"),
                                epsilon = 1e-3) {

  method <- match.arg(method)

  fleet_names <- colnames(e_p_fl)
  patches     <- nrow(e_p_fl)
  n_fleets    <- ncol(e_p_fl)

  # --- baseline evaluation (reuse if caller already has it) ------------------
  if (is.null(baseline)) {
    baseline <- go_fish(
      e_p_fl        = e_p_fl,
      fauna         = fauna,
      n_p_a         = n_p_a,
      fleets        = fleets,
      output_format = "matrix"
    )
  }

  base_r    <- baseline$r_p_fl      # patches x fleets
  base_prof <- baseline$prof_p_fl

  # --- output matrices -------------------------------------------------------
  mr_p_fl <- matrix(0, nrow = patches, ncol = n_fleets,
                    dimnames = dimnames(e_p_fl))
  mp_p_fl <- matrix(0, nrow = patches, ncol = n_fleets,
                    dimnames = dimnames(e_p_fl))

  if (method == "separable") {
    # -- Separable: bump ALL patches for one fleet at a time ------------------
    for (fl in seq_len(n_fleets)) {

      e_bumped       <- e_p_fl
      e_bumped[, fl] <- e_bumped[, fl] + epsilon

      bumped <- go_fish(
        e_p_fl        = e_bumped,
        fauna         = fauna,
        n_p_a         = n_p_a,
        fleets        = fleets,
        output_format = "matrix"
      )

      mr_p_fl[, fl] <- (bumped$r_p_fl[, fl]    - base_r[, fl])    / epsilon
      mp_p_fl[, fl] <- (bumped$prof_p_fl[, fl]  - base_prof[, fl]) / epsilon
    }

  } else {
    # -- Patch loop: bump ONE patch-fleet cell at a time ----------------------
    for (fl in seq_len(n_fleets)) {
      for (p in seq_len(patches)) {

        e_bumped        <- e_p_fl
        e_bumped[p, fl] <- e_bumped[p, fl] + epsilon

        bumped <- go_fish(
          e_p_fl        = e_bumped,
          fauna         = fauna,
          n_p_a         = n_p_a,
          fleets        = fleets,
          output_format = "matrix"
        )

        mr_p_fl[p, fl] <- (bumped$r_p_fl[p, fl]    - base_r[p, fl])    / epsilon
        mp_p_fl[p, fl] <- (bumped$prof_p_fl[p, fl]  - base_prof[p, fl]) / epsilon
      }
    }
  }

  list(
    mr_p_fl = mr_p_fl,
    mp_p_fl = mp_p_fl
  )
}

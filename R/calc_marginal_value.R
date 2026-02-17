#' Calculate Marginal Revenue or Profit by Patch and Fleet
#'
#' Uses forward finite differences on \code{\link{go_fish}} to approximate
#' \eqn{\partial R / \partial E_{p,f}} (or the profit analogue) for every
#' patch–fleet combination.  The output is a list of \code{patches × fleets}
#' matrices designed to be merged directly into the \code{buffet} list consumed
#' by \code{\link{allocate_effort}}.
#'
#' Two strategies are available via \code{method}:
#' \describe{
#'   \item{\code{"separable"}}{Bumps effort in **all** patches simultaneously
#'     by \code{epsilon} and computes derivatives in a single
#'     \code{\link{go_fish}} call per fleet.
#'     Fast (\emph{n_fleets} extra evaluations) but assumes the marginal
#'     return in patch \emph{i} is independent of the effort perturbation in
#'     patch \emph{j}, conditional on covariates.}
#'   \item{\code{"patch_loop"}}{Bumps effort in **one patch at a time**,
#'     requiring \emph{patches × fleets} extra evaluations.  Slower but fully
#'     accounts for cross-patch depletion effects.}
#' }
#'
#' @param e_p_fl Numeric matrix of effort by patch (rows) and fleet (columns).
#' @param fauna A list of fauna objects (as in \code{go_fish}).
#' @param n_p_a A list of numbers-at-age by patch (as in \code{go_fish}).
#' @param fleets A list of fleet objects (as in \code{go_fish}).
#' @param baseline Optional. The list returned by a previous call to
#'   \code{go_fish(output_format = "matrix")} with the same inputs.  Passing
#'   this avoids one redundant evaluation.
#' @param method Character, one of \code{"separable"} or \code{"patch_loop"}.
#'   See \strong{Details}.
#' @param epsilon Numeric scalar; the additive effort perturbation used for the
#'   finite difference.
#'
#' @return A named list with two elements, each a \code{patches × fleets}
#'   matrix with the same dimension names as \code{e_p_fl}:
#'   \describe{
#'     \item{\code{mr_p_fl}}{Marginal revenue
#'       (\eqn{\partial \text{Revenue} / \partial E_{p,f}}).}
#'     \item{\code{mp_p_fl}}{Marginal profit
#'       (\eqn{\partial \text{Profit} / \partial E_{p,f}}).}
#'   }
#'   These can be inserted directly into the \code{buffet} list and picked up
#'   by \code{\link{allocate_effort}} when a fleet's
#'   \code{spatial_allocation} is \code{"marginal_revenue"} or
#'   \code{"marginal_profit"}.
#'
#' @export
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

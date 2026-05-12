#' Create a Fleet Object
#'
#' @description
#' Builds a fleet object specifying the fishing behaviour, spatial dynamics,
#' and economic parameters for one group of vessels targeting one or more
#' species. Always run \code{\link{tune_fleets}} after creating fleets to
#' calibrate catchability and costs to the desired initial conditions.
#'
#' @details
#' ## Fleet models
#' \describe{
#'   \item{\code{"constant_effort"}}{Total effort is fixed at \code{base_effort}
#'     each time step. Use for scenarios where fishing pressure is prescribed
#'     externally.}
#'   \item{\code{"open_access"}}{Total effort adjusts each step based on a
#'     normalised average-profitability signal. Equilibrates where total profits
#'     = 0. Entry/exit speed is controlled by \code{oa_max_growth_per_year},
#'     \code{oa_max_decline_per_year}, and \code{oa_signal_half}.}
#'   \item{\code{"sole_owner"}}{Identical dynamics to \code{"open_access"} but
#'     uses the marginal (not average) profit signal. Equilibrates at MEY where
#'     marginal profit = 0. Requires \code{calc_marginal_value()} to be
#'     computed each step; \code{\link{simmar}} handles this automatically.}
#'   \item{\code{"manual"}}{Total effort each step is taken directly from a
#'     user-supplied vector; see the \code{manager} argument of
#'     \code{\link{simmar}}.}
#' }
#'
#' ## Spatial allocation
#' The \code{spatial_allocation} argument determines how total fleet effort is
#' distributed among patches each step. Options:
#' \describe{
#'   \item{\code{"rpue"}}{Revenue per unit effort (default). Effort concentrates
#'     in high-revenue patches.}
#'   \item{\code{"revenue"}}{Total revenue. Similar to \code{"rpue"} but favours
#'     larger patches.}
#'   \item{\code{"ppue"}}{Profit per unit effort (cost-aware).}
#'   \item{\code{"profit"}}{Total profit (cost-aware).}
#'   \item{\code{"cpue"} / \code{"catch"}}{Catch-based variants.}
#'   \item{\code{"marginal_revenue"} / \code{"marginal_profit"}}{Uses finite-
#'     difference marginal returns from \code{\link{calc_marginal_value}}.
#'     Required for \code{"sole_owner"} fleets. Requires
#'     \code{fleet_model = "sole_owner"} or explicit pre-computation.}
#'   \item{\code{"manual"}}{Effort distributed proportionally to continuous
#'     weights in \code{fishing_grounds$fishing_ground} (0--1 valued).}
#'   \item{\code{"uniform"}}{Effort spread equally across all open patches.}
#' }
#'
#' ## Costs
#' Total cost per fleet is:
#' \deqn{C = c_0 \, E^{ref} \sum_l \left[\left(\frac{E_l}{E^{ref}}\right)^\gamma + \theta \, \tilde{d}_l \frac{E_l}{E^{ref}}\right]}
#' where \eqn{c_0} is \code{cost_per_unit_effort}, \eqn{E^{ref}} is the
#' reference effort per patch, \eqn{\gamma} is \code{effort_cost_exponent},
#' \eqn{\theta} is \code{travel_weight} (derived from \code{travel_fraction}),
#' and \eqn{\tilde{d}_l} is the normalised distance from patch \eqn{l} to the
#' nearest port.
#'
#' @param metiers Named list of \code{\link{Metier}} R6 objects, one per
#'   species in \code{fauna}. Each metier specifies price, selectivity, and
#'   relative catchability for that fleet-species combination. Names must match
#'   species names in \code{fauna}.
#' @param fleet_model Character. Effort dynamics model; see Details.
#'   One of \code{"constant_effort"} (default), \code{"open_access"},
#'   \code{"sole_owner"}, or \code{"manual"}.
#' @param base_effort Numeric. Total effort units available to the fleet.
#'   Defaults to \code{prod(resolution)} (one unit per patch). Catchability is
#'   calibrated relative to this value by \code{\link{tune_fleets}}.
#' @param spatial_allocation Character. Spatial effort allocation strategy;
#'   see Details. Default \code{"rpue"}.
#' @param cr_ratio Numeric. Target cost-to-revenue ratio at equilibrium.
#'   \code{1} implies zero profits (open-access equilibrium). Used by
#'   \code{\link{tune_fleets}} to calibrate \code{cost_per_unit_effort}.
#' @param effort_cost_exponent Numeric. Exponent \eqn{\gamma} in the effort
#'   cost function, controlling congestion / convexity. Values > 1 make
#'   additional units of effort progressively more expensive. Default \code{1.2}.
#' @param travel_fraction Numeric in [0, 1). Fraction of total costs
#'   attributable to travel at equilibrium. \code{0} means no spatial cost
#'   heterogeneity (all patches equally costly). Controls how strongly port
#'   proximity shapes the spatial cost surface.
#' @param ports Data frame with columns \code{x} and \code{y} giving port
#'   patch coordinates. Minimum distances from each patch to the nearest port
#'   are used to compute the travel-cost component when
#'   \code{travel_fraction > 0}.
#' @param mpa_response Character. Vessel response to MPA closures: \code{"stay"}
#'   (vessels stay in the closed area) or \code{"leave"} (vessels redistribute
#'   to open patches, concentrating effort).
#' @param resolution Integer scalar or length-2 integer vector \code{c(nx, ny)}.
#'   Must match the resolution of the \code{fauna} objects.
#' @param patch_area Numeric. Area of each patch (km^2). Used to compute
#'   port distances.
#' @param fishing_grounds Data frame with columns \code{x}, \code{y}, and
#'   \code{fishing_ground} (logical or numeric). Restricts where effort can be
#'   deployed; \code{TRUE}/\code{1} = open, \code{FALSE}/\code{0} = closed.
#'   When \code{spatial_allocation = "manual"}, numeric values in
#'   \code{fishing_ground} are used as effort weights. Defaults to all patches
#'   open.
#' @param oa_max_growth_per_year Numeric. Maximum fractional increase in total
#'   effort per year under highly profitable conditions for open-access and
#'   sole-owner fleets (e.g. \code{0.5} = +50% per year). Converted to a
#'   per-step multiplier internally. Must be > 0.
#' @param oa_max_decline_per_year Numeric in (0, 1). Maximum fractional
#'   decrease in total effort per year under highly unprofitable conditions
#'   (e.g. \code{0.5} = -50% per year). Must be in (0, 1).
#' @param oa_signal_half Numeric in (0, 1). Profitability signal value at which
#'   effort adjustment speed reaches half its maximum. Lower values make fleets
#'   more sensitive. Typical range \code{0.15}--\code{0.4}. Default \code{0.3}.
#' @param cost_per_unit_effort Numeric. Base cost per unit of effort. Overridden
#'   by \code{\link{tune_fleets}} when \code{tune_costs = TRUE}; rarely needs
#'   manual adjustment.
#' @param cost_per_distance Numeric. Deprecated; use \code{travel_fraction}
#'   instead.
#' @param responsiveness Numeric. Per-step responsiveness of patch effort to
#'   the objective signal in \code{\link{allocate_effort}} (the \eqn{\eta}
#'   parameter of the multiplicative update). Larger values move effort more
#'   aggressively toward high-objective patches each step; if too large, can
#'   drive period-2 sawtooth oscillation. Default \code{0.025}.
#' @param memory_halflife Non-negative numeric. Half-life **in
#'   years** of the exponential smoothing applied to this fleet's spatial objective surface
#'   inside \code{\link{simmar}}. \code{0} (default) disables smoothing — the
#'   fleet sees only the previous step's objective, matching legacy behavior.
#'   The parameter is season-agnostic: simmar converts it internally to time
#'   steps (\code{halflife_steps = halflife * steps_per_year}) so a given
#'   value produces the same calendar-time smoothing regardless of how many
#'   seasons per year the model uses. Larger values dampen high-frequency
#'   feedback oscillations by blending in past objective surfaces; the
#'   per-step weight on the current surface is
#'   \eqn{\alpha = 1 - 0.5^{1/halflife_{steps}}}. Smoothing updates only
#'   patches that are currently open; closed patches retain their last-seen
#'   smoothed value ("freeze and resume"). Early warm-up steps use a
#'   Welford-style ramp (effective
#'   \eqn{\alpha_\mathrm{eff} = \max(\alpha, 1/n)}) so the smoothed surface
#'   is not anchored to the first observed buffet column.
#'
#'   Practical guidance: values around \code{0.5}–\code{1.5} years are the
#'   typical sweet spot. Larger halflives introduce phase lag of roughly
#'   \eqn{1.44 \times \mathrm{halflife}} years between a true change in patch
#'   marginal value and the fleet's perceived value, which can produce
#'   low-frequency overshoot/undershoot — a different pathology from the
#'   high-frequency sawtooth that motivates the parameter. If catch
#'   trajectories under a given halflife show slow swings that aren't present
#'   at \code{halflife = 0}, reduce it.
#'
#' @return A named list (fleet object) with all parameters needed by
#'   \code{\link{simmar}} and \code{\link{tune_fleets}}, including
#'   computed travel weights, normalised cost-per-patch, and (for open-access /
#'   sole-owner fleets) the derived annual effort-adjustment parameters.
#'
#' @seealso \code{\link{tune_fleets}}, \code{\link{simmar}},
#'   \code{\link{create_critter}}, \code{\link{Metier}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a metier for a single species
#' met <- Metier$new(
#'   critter     = fauna[["tuna"]],
#'   price       = 10,
#'   sel_form    = "logistic",
#'   sel_start   = 0.3,
#'   sel_delta   = 0.1,
#'   catchability = 0.01,
#'   p_explt     = 1
#' )
#'
#' # Constant-effort fleet
#' fleet <- create_fleet(
#'   metiers    = list(tuna = met),
#'   resolution = c(10, 10)
#' )
#'
#' # Open-access fleet with port-based travel costs
#' ports <- data.frame(x = 1, y = 1)
#' oa_fleet <- create_fleet(
#'   metiers              = list(tuna = met),
#'   fleet_model          = "open_access",
#'   spatial_allocation   = "ppue",
#'   travel_fraction      = 0.3,
#'   ports                = ports,
#'   cr_ratio             = 0.9,
#'   resolution           = c(10, 10)
#' )
#'
#' fleets <- list(fleet = fleet)
#' fleets <- tune_fleets(fauna, fleets, tune_type = "depletion")
#' }
create_fleet <-
  function(metiers,
           mpa_response = "stay",
           fleet_model = "constant_effort",
           oa_max_growth_per_year = 0.5,
           oa_max_decline_per_year = 0.5,
           oa_signal_half = 0.3,
           cost_per_unit_effort = 1,
           spatial_allocation = "rpue",
           effort_cost_exponent = 1.2,
           travel_fraction = 0,
           ports = NULL,
           cost_per_distance = 1,
           cr_ratio = 1,
           resolution,
           patch_area = 1,
           base_effort = NULL,
           fishing_grounds = NULL,
           responsiveness = 0.025,
           memory_halflife = 0) {

    fleet_model <- stringr::str_replace_all(fleet_model, " ", "_") # in case someone used spaces accidentally (like dumbass old dan)

    if (!is.numeric(memory_halflife) ||
        length(memory_halflife) != 1L ||
        is.na(memory_halflife) ||
        memory_halflife < 0) {
      stop("`memory_halflife` must be a single non-negative number (0 disables smoothing).")
    }

    if (length(resolution) == 1) {
      resolution <- rep(resolution, 2)
    }
    if (is.null(base_effort)) {
      base_effort <- prod(resolution)
    }

    if (is.null(fishing_grounds)) {
      fishing_grounds <- tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) |>
        dplyr::mutate(fishing_ground = TRUE) |>
        dplyr::arrange(x, y)
    }

    fishing_grounds <- fishing_grounds |>
      dplyr::arrange(x, y)

    if (nrow(fishing_grounds) != prod(resolution)) {
      stop("supplied fishing_grounds do not match the spatial dimensions of the simualted domain. make sure that number of rows and columns match supplied resolution.")
    }

    # Count number of open patches (fishing_ground > 0)
    n_open_patches <- sum(fishing_grounds$fishing_ground > 0)

    # Set default effort_reference as mean per-patch effort across open patches
    effort_reference <- base_effort / n_open_patches

    # Initialize cost_per_patch based on ports or as uniform
    if (is.null(ports)) {

      cost_per_patch <- rep(0, prod(resolution))

    } else {

      patches <- tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) |>
        dplyr::mutate(patch = 1:dplyr::n()) # extra step to make sure patch ordering is consistent

      if (any(ports$x > resolution[1]) | any(ports$y > resolution[2])) {
        stop("one or more port locations is outside of spatial grid")
      }

      ports <- ports |>
        dplyr::left_join(patches, by = c("x", "y"))

      # calculate the distance between each of the patches
      port_distance <- tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) |>
        dist(diag = TRUE) |>
        as.matrix() * sqrt(mean(patch_area))

      # minimum distance from each patch to any port patch
      port_distances <- apply(
        matrix(port_distance[ports$patch, ], nrow = length(ports$patch)),
        2,
        min
      )

      cost_per_patch <- port_distances * cost_per_distance
    }

    # Normalize cost_per_patch to mean 1 across open patches
    # This makes it a pure spatial shape, interpretable with travel_fraction
    open_patches <- fishing_grounds$fishing_ground > 0

    # Handle edge case: if cost_per_patch is all zeros, convert to travel_fraction = 0
    if (all(cost_per_patch == 0)) {

      travel_fraction <- 0
      cost_per_patch_normalized <- rep(1, length(cost_per_patch))

    } else {

      mean_cost_open <- mean(cost_per_patch[open_patches])
      if (mean_cost_open > 0) {
        cost_per_patch_normalized <- cost_per_patch / mean_cost_open
      } else {
        cost_per_patch_normalized <- rep(1, length(cost_per_patch))
      }
    }

    # Compute travel_weight (theta) from travel_fraction
    # theta = travel_fraction / (1 - travel_fraction)
    # This parameterization ensures travel costs are travel_fraction of total costs
    if (travel_fraction < 0 || travel_fraction >= 1) {
      stop("travel_fraction must be in [0, 1)")
    }

    if (travel_fraction == 0) {
      travel_weight <- 0
    } else {
      travel_weight <- travel_fraction / (1 - travel_fraction)
    }

    # ---------------------------------------------------------------------
    # open access parameters (annual, intuitive knobs)
    #
    # simmar() will convert annual quantities to per-step values using time_step:
    #   rho_step = oa_rho_year * time_step
    #   m_max_step = oa_m_max_year ^ time_step
    #   m_min_step = oa_m_min_year ^ time_step
    #
    # and will use:
    #   multiplier = exp(rho_step * tanh(oa_signal / oa_k))
    #
    # where oa_signal is the normalized profitability signal in [-1, 1]
    # ---------------------------------------------------------------------

    oa_m_max_year <- NULL
    oa_m_min_year <- NULL
    oa_rho_year <- NULL
    oa_k <- NULL

    if (fleet_model %in% c("open_access", "sole_owner")) {

      if (oa_max_growth_per_year <= 0) {
        stop("oa_max_growth_per_year must be > 0 (e.g. 0.5 for +50% per year)")
      }
      if (oa_max_decline_per_year <= 0 || oa_max_decline_per_year >= 1) {
        stop("oa_max_decline_per_year must be in (0, 1) (e.g. 0.5 for -50% per year)")
      }
      if (oa_signal_half <= 0 || oa_signal_half >= 1) {
        stop("oa_signal_half must be in (0, 1) (typical values are 0.2-0.4)")
      }

      # annual caps expressed as multipliers
      oa_m_max_year <- 1 + oa_max_growth_per_year
      oa_m_min_year <- 1 - oa_max_decline_per_year

      # choose oa_rho_year so that when tanh() saturates (-> 1), multiplier_year = exp(oa_rho_year) = oa_m_max_year
      oa_rho_year <- log(oa_m_max_year)

      # choose oa_k so that tanh(oa_signal_half / oa_k) = 0.5 (half of maximum response speed)
      oa_k <- oa_signal_half / atanh(0.5)
    }

    fleet <- list(
      metiers = metiers,
      base_effort = base_effort,
      mpa_response = mpa_response,
      cr_ratio = cr_ratio,
      cost_per_unit_effort = cost_per_unit_effort,
      spatial_allocation = spatial_allocation,
      fleet_model = fleet_model,
      effort_cost_exponent = effort_cost_exponent,
      travel_fraction = travel_fraction,
      travel_weight = travel_weight,
      cost_per_patch = cost_per_patch,  # Keep original for backwards reference
      cost_per_patch_normalized = cost_per_patch_normalized,
      effort_reference = effort_reference,
      fishing_grounds = fishing_grounds,

      # open access knobs (annual, intuitive)
      oa_max_growth_per_year = oa_max_growth_per_year,
      oa_max_decline_per_year = oa_max_decline_per_year,
      oa_signal_half = oa_signal_half,

      # open access derived params used by simmar()
      oa_m_max_year = oa_m_max_year,
      oa_m_min_year = oa_m_min_year,
      oa_rho_year = oa_rho_year,
      oa_k = oa_k,

      responsiveness = responsiveness,

      memory_halflife = memory_halflife
    )

    return(fleet)
  }

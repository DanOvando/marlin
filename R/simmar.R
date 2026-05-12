#' Run the marlin Spatially Explicit Fisheries Simulation
#'
#' @description
#' \code{simmar} is the main simulation engine for the marlin package.
#' Given initialised \code{fauna} and \code{fleets} objects, it advances the
#' age-structured, spatially explicit population forward in time for a
#' specified number of years, applying fishing mortality, natural mortality,
#' movement, and recruitment each time step.
#'
#' @details
#' ## Simulation loop (per time step)
#' \enumerate{
#'   \item **Effort dynamics**: For each fleet, total effort for this step is
#'     determined by the fleet model (\code{"constant_effort"},
#'     \code{"open_access"}, \code{"sole_owner"}, or \code{"manual"}).
#'     Open-access and sole-owner fleets use a tanh-based normalised
#'     profitability signal derived from the previous step's revenues and costs.
#'   \item **Spatial allocation**: \code{\link{allocate_effort}} redistributes
#'     each fleet's total effort across patches using a multiplicative
#'     velocity-field update driven by the chosen \code{spatial_allocation}
#'     metric (\code{"rpue"}, \code{"profit"}, etc.) from the previous step's
#'     buffet.
#'   \item **Fishing & population dynamics**: For each critter,
#'     age-structured fishing mortality \eqn{F_{p,a}} is applied, followed by
#'     natural mortality, movement, and Beverton-Holt recruitment.
#'   \item **Yield accounting**: \code{\link{allocate_yields}} and
#'     \code{\link{aggregate_yields}} build the per-fleet buffet of catch,
#'     revenue, profit, and CPUE by patch.
#'   \item **Marginal values (if needed)**: When any fleet uses
#'     \code{spatial_allocation = "marginal_profit"} or
#'     \code{"marginal_revenue"}, or \code{fleet_model = "sole_owner"},
#'     \code{\link{calc_marginal_value}} appends patch-level marginal returns
#'     to the buffet.
#'   \item **Quota enforcement**: If \code{manager$quotas} is set for a
#'     species, a scalar effort multiplier is solved via \code{optim} to keep
#'     total catch at or below the quota.
#' }
#'
#' ## Step naming
#' Steps are named \code{"year_season"} (e.g. \code{"5_3"} = year 5, season
#' 3). Use \code{\link{clean_steps}} to strip any \code{"step_"} prefix.
#'
#' ## Initial conditions
#' By default the simulation starts from the unfished equilibrium embedded in
#' each critter object. Pass \code{initial_conditions = sim[[length(sim)]]}
#' from a prior run to chain simulations (e.g. a baseline followed by an MPA
#' scenario).
#'
#' @param fauna Named list of critter objects from \code{\link{create_critter}}.
#'   All critters must share the same \code{time_step} (i.e. the same number
#'   of \code{seasons}).
#' @param fleets Named list of fleet objects from \code{\link{create_fleet}},
#'   tuned with \code{\link{tune_fleets}}.
#' @param manager Named list of management actions. Supported elements:
#'   \describe{
#'     \item{\code{mpas}}{List with \code{$locations} (a data frame with
#'       columns \code{x}, \code{y}, \code{mpa}; \code{mpa = TRUE} = closed)
#'       and \code{$mpa_year} (integer year when closures are activated).}
#'     \item{\code{quotas}}{Named numeric vector of annual catch caps per
#'       species (names match \code{fauna}). Effort is scaled down by an
#'       optimised multiplier when the cap would be exceeded.}
#'     \item{\code{effort_cap}}{Named numeric vector capping total effort per
#'       fleet in open-access / sole-owner models.}
#'     \item{\code{closed_seasons}}{Named list of integer vectors specifying
#'       which seasons are closed to fishing for each critter.}
#'   }
#' @param habitat Named list (one entry per critter) of time-varying habitat.
#'   Each entry is a list of \code{[ny, nx]} matrices (one per year or per
#'   time step). When provided, the movement matrix is updated each step based
#'   on the difference in habitat quality between adjacent patches.
#' @param years Numeric. Number of years to simulate. Ignored when \code{steps}
#'   is provided.
#' @param steps Integer. Alternative to \code{years}: number of time steps to
#'   simulate. Takes precedence over \code{years} when not \code{NA}.
#' @param starting_season Integer. Starting season within a year. Rarely
#'   needed; used when chaining partial-year runs.
#' @param initial_conditions Initial population state, typically
#'   \code{sim[[length(sim)]]} from a previous \code{simmar} call. Defaults to
#'   the unfished equilibrium embedded in each critter object.
#' @param starting_step Character. Step name of the form \code{"year_season"}
#'   to start from. Controls step labelling when chaining runs.
#' @param keep_starting_step Logical. If \code{TRUE} (default), the starting
#'   (initial-conditions) step is included in the output.
#' @param log_rec_devs Numeric matrix of pre-generated log recruitment
#'   deviates with rows = steps and columns = critters (names must match
#'   \code{fauna}). When \code{NULL} (default), deviates are generated
#'   internally using each critter's \code{sigma_rec} and \code{ac_rec}.
#' @param cor_rec Numeric \eqn{n \times n} correlation matrix for recruitment
#'   deviates across species. Default \code{diag(length(fauna))} (independent
#'   recruitment).
#'
#' @return A named list of length \code{years / time_step} (or \code{steps}),
#'   where each element is named \code{"year_season"} and contains a named
#'   sub-list (one per critter) with the population state at that step.
#'   Each critter's state includes:
#'   \describe{
#'     \item{\code{n_p_a}}{Numbers by patch and age (matrix).}
#'     \item{\code{b_p_a}}{Biomass by patch and age (matrix).}
#'     \item{\code{ssb_p_a}}{Spawning stock biomass by patch and age (matrix).}
#'     \item{\code{c_p_a}}{Catch in numbers by patch and age (matrix).}
#'     \item{\code{c_p_fl}}{Catch by patch and fleet (matrix).}
#'     \item{\code{r_p_fl}}{Revenue by patch and fleet (matrix).}
#'     \item{\code{prof_p_fl}}{Profit by patch and fleet (matrix).}
#'     \item{\code{e_p_fl}}{Effort by patch and fleet (data frame).}
#'   }
#'   Pass output to \code{\link{process_marlin}} for tidying.
#'
#' @seealso \code{\link{process_marlin}}, \code{\link{plot_marlin}},
#'   \code{\link{create_critter}}, \code{\link{create_fleet}},
#'   \code{\link{tune_fleets}}, \code{\link{go_fish}},
#'   \code{\link{allocate_effort}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic run
#' sim <- simmar(fauna = fauna, fleets = fleets, years = 50)
#'
#' # With MPA starting in year 10
#' mpa_locs <- expand_grid(x = 1:10, y = 1:10) |>
#'   mutate(mpa = x <= 5)
#' sim_mpa <- simmar(
#'   fauna   = fauna,
#'   fleets  = fleets,
#'   years   = 50,
#'   manager = list(
#'     mpas = list(locations = mpa_locs, mpa_year = 10)
#'   )
#' )
#'
#' # Chain: baseline + policy scenario
#' sim_base    <- simmar(fauna = fauna, fleets = fleets, years = 50)
#' sim_policy  <- simmar(fauna = fauna, fleets = fleets, years = 20,
#'                       initial_conditions = sim_base[[length(sim_base)]],
#'                       manager = list(quotas = list(tuna = 500)))
#' }
simmar <- function(fauna = list(),
                   fleets = list(),
                   manager = list(),
                   habitat = list(),
                   years = 100,
                   steps = NA,
                   starting_season = NA,
                   initial_conditions = NA,
                   starting_step = NA,
                   keep_starting_step = TRUE,
                   log_rec_devs = NULL,
                   cor_rec = diag(length(fauna))) {

  init_cond_provided <-
    !all(is.na(initial_conditions)) # marker in case initial conditions were provided

  fauni <- names(fauna)

  # if initial_conditions is passed, enforce names/order too (cheap and catches a ton)
  if (!all(is.na(initial_conditions))) {
    initial_conditions <- rlang::set_names(initial_conditions, fauni)
    initial_conditions <- initial_conditions[fauni]
  }

  if (is.null(fauni) || anyNA(fauni) || any(fauni == "") || any(duplicated(fauni))) {
    stop("fauna must be a *named* list with unique, non-empty species names.")
  }

  n_critters <- length(fauna)

  fleet_names <- names(fleets)

  time_step <- unique(purrr::map_dbl(fauna, "time_step"))
  if (length(time_step) != 1) {
    stop("all fauna must share the same time_step")
  }
  steps_per_year <- 1 / time_step


  steps_per_year <- 1 / time_step

  patch_area_vec <- fauna[[1]]$grid$patch_area

  sigma_recs <- purrr::map_dbl(fauna, "sigma_rec") # gather recruitment standard deviations

  ac_recs <- purrr::map_dbl(fauna, "ac_rec") # gather autocorrelation in recruitment standard deviations

  if (any(abs(ac_recs) >= 1)) {
    stop("all sigma_ac values must be between -1 and 1")
  }

  covariance_rec <- cor_rec * (sigma_recs %o% sigma_recs)

  if (length(time_step) > 1) {
    stop(
      paste(
        "All critters in fauna must have the same time step: current time steps are",
        paste(time_step, collapse = " ")
      )
    )
  }

  if (is.na(steps)) {
    steps <-
      (years) / time_step + 1 # tack on extra step for accounting
  } else {
    steps <- steps + 1 # to store initial conditions
  } # if steps are specified instead of years

  steps <- pmax(steps, 3) # need to be at least 1 initial + 2 running steps

  # generate recruitment deviates
  if (!is.na(starting_step)) {
    year_season <- marlin::clean_steps(starting_step)

    starting_year <- as.integer(sub("_.*$", "", year_season)) - 1
    starting_season <- as.integer(sub("^.*_", "", year_season)) - 1

    offset <- (starting_year * steps_per_year) + starting_season
  } else {
    offset <- 0
  }

  # create step names of form "year_season"
  step_names <- paste(
    rep(1:(steps + offset), each = steps_per_year),
    1:steps_per_year,
    sep = "_"
  )
  step_names <- step_names[(1:steps) + offset]

  # --- SPEED: parse year/season ONCE (avoid repeated stringr in the step loop) ---
  step_year   <- as.integer(sub("_.*$", "", step_names))
  step_season <- as.integer(sub("^.*_", "", step_names))
  # --- end step parsing ---

  # recruitment deviates
  if (is.null(log_rec_devs)) {
    log_rec_devs <- matrix(
      NA_real_,
      nrow = length(step_names),
      ncol = n_critters,
      dimnames = list(step_names, fauni)
    )

    # SPEED: draw all MVN innovations once, then do AR(1) recursion
    eps <- mvtnorm::rmvnorm(length(step_names), mean = rep(0, n_critters), sigma = covariance_rec)
    log_rec_devs[1, ] <- eps[1, ]
    if (length(step_names) > 1) {
      ar_scale <- sqrt(1 - ac_recs^2)
      for (i in 2:length(step_names)) {
        log_rec_devs[i, ] <- ac_recs * log_rec_devs[i - 1, ] + ar_scale * eps[i, ]
      }
    }
  } else {

    if (is.null(colnames(log_rec_devs)) ||
        !setequal(colnames(log_rec_devs), names(fauna))) {
      stop("column names of log_rec_devs must match those of fauna")
    }
    log_rec_devs <- log_rec_devs[, names(fauna), drop = FALSE]


  }

  # step_names <- step_names[1:steps] # chop back to the actual number of steps used; works this way in case there are multiple steps per year but an incomplete number of years
  patches <- unique(purrr::map_dbl(fauna, "patches"))

  # not building checks for same resolution here since almost redundant to patches unless product of resolution is identical
  resolution <-
    (purrr::map(fauna[1], ~ data.frame(t(.x$resolution)))) |>
    purrr::list_rbind() |>
    unlist()

  if (all(is.na(initial_conditions))) {
    initial_conditions <-
      purrr::map(fauna, c("unfished")) # pull out unfished conditions created by create_critter

    initial_e_p_fl <- purrr::imap(fleets, ~ data.frame(x = rep(.x$base_effort / patches, patches)) |>
                                    setNames(.y)) |>
      purrr::list_cbind() |>
      as.matrix()


    initial_n_p_a <- purrr::map(fauna, "n_p_a_0")

    initial_exploration <- go_fish(
      e_p_fl = initial_e_p_fl,
      fauna = fauna,
      fleets = fleets,
      n_p_a = initial_n_p_a,
      output_format = "matrix"
    )

    open_patches <- matrix(nrow = patches, ncol = length(fleets))

    for (i in seq_along(fleets)){
      open_patches[,i] <- as.numeric(fleets[[i]]$fishing_grounds$fishing_ground)
    }

    needs_marginals <- any(
      purrr::map_chr(fleets, "spatial_allocation") %in%
        c("marginal_revenue", "marginal_profit")
    ) || any(
      purrr::map_chr(fleets, "fleet_model") == "sole_owner"
    )

    if (needs_marginals) {
      # Gather current n_p_a from storage for go_fish
      marginal_n_p_a <- initial_n_p_a

      marginals <- calc_marginal_value(
        e_p_fl   = initial_e_p_fl,
        fauna    = fauna,
        n_p_a    = marginal_n_p_a,
        fleets   = fleets,
        baseline = NULL,          # could pass go_fish result if available
        method   = "separable",   # fast default; "patch_loop" for small grids
        epsilon  = 1e-3
      )

      initial_exploration$mr_p_fl <- marginals$mr_p_fl
      initial_exploration$mp_p_fl <- marginals$mp_p_fl
    }

    initial_effort <- allocate_effort(
      effort_by_patch = initial_e_p_fl,
      fleets = fleets,
      buffet = initial_exploration,
      open_patch = open_patches,
      flatness_tol = 1e-2
    )

  }

  if (length(patches) > 1) {
    stop(
      "fauna have different habitat resolutions: set resolution to same number for all species!"
    )
  }

  #  determine open fishing seasons

  season_foo <- function(fauni, name, manager) {
    open_seasons <- 1:fauni$seasons

    if (length(manager$closed_seasons[name]) > 0) {
      open_seasons <-
        open_seasons[!(open_seasons %in% manager$closed_seasons[[name]])]
    }

    return(open_seasons)
  }

  fishing_seasons <-
    purrr::imap(fauna, season_foo, manager = manager) # figure out which seasons are open for each critter

  if (!is.null(names(fishing_seasons)) && all(fauni %in% names(fishing_seasons))) {
    fishing_seasons <- fishing_seasons[fauni]
  }

  # --- Habitat handling (same logic, faster inner use) ---
  # We keep your expansion logic exactly, but ALSO precompute the per-step habitat vector
  # in the SAME patch order you were enforcing with:
  #   pivot_longer(as.data.frame(mat), everything(), names_to="x", names_transform=list(x=as.integer)) |> arrange(x)
  #
  # NOTE (spirit-of-comment): because pivot_longer over a data.frame stacks *columns* in order,
  # and arrange(x) doesn’t change that column order, the resulting $value is identical to as.vector(mat).
  habitat_vecs <- list()

  if (length(habitat) > 0) {
    if (!any(names(habitat) %in% fauni)) {
      stop(
        "names of critters in habitat must must match name of at least one critter in fauni list"
      )
    }

    for (f in seq_along(habitat)) {
      if (length(habitat[[f]]) != years &
          length(habitat[[f]]) != (steps - 1)) {
        stop(
          "supplied habitat vector must be same length as either number of years or number of years times number of seasons"
        )
      } else if (length(habitat[[f]]) != (steps - 1)) {
        new_habitat <- vector(mode = "list", length = steps - 1)

        for (i in 1:(steps - 1)) {
          supplied_years <- length(habitat[[f]])
          if (!is.na(starting_step)) {
            years_in <-
              as.integer(sub("_.*$", "", step_names[i])) - as.integer(sub("_.*$", "", starting_step)) + 1
            # put year in index form not named form
          } else {
            years_in <-
              as.integer(sub("_.*$", "", step_names[i]))
          }
          years_in <- min(years_in, supplied_years)

          # year <-  floor(step_names[i] - starting_step) # put year in index form not named form
          new_habitat[[i]] <- habitat[[f]][[years_in]]
        }

        habitat[[f]] <- new_habitat
      } # close expansion of habitat if needed
    } # close fauni loop

    # Precompute vectors once (speeds up the step loop)
    for (sp in intersect(names(habitat), fauni)) {
      habitat_vecs[[sp]] <- lapply(habitat[[sp]], function(hmat) {
        as.vector(hmat) # equivalent ordering to your pivot_longer(... ) |> arrange(x)
      })
    }
  } # close check on supplied habitat
  # --- end habitat handling ---

  storage <- replicate(
    steps,
    setNames(vector("list", length(fauna)), names(fauna)),
    simplify = FALSE
  )

  storage[[1]] <-
    initial_conditions # start populations at initial conditions


  # fleets <- purrr::map(
  #   fleets,
  #   ~ purrr::list_modify(.x, e_p_s = matrix((.x$base_effort / patches), nrow = patches, ncol = steps))
  # )


  last_e_p_f <- matrix(NA,nrow = patches, ncol = length(fleets))

  if (init_cond_provided) {

    last_e_p_f <- initial_conditions[[1]]$e_p_fl

  } else {
    last_e_p_f <- initial_effort$effort_new
  }

  # pre_quota_e_p_f stores the fleet's "intended" effort before any quota

  # reduction. For constant_effort fleets, this is used as the starting point

  # each step so that quota reductions are transient rather than permanent.
  # For open_access / sole_owner fleets, effort rebounds naturally via
  # profitability dynamics so they use last_e_p_f (realized effort) instead.
  pre_quota_e_p_f <- last_e_p_f

  r_p_f <- matrix(0, patches, length(fauni))
  e_p_f <- matrix(0, patches, length(fauni))
  f_q <- rep(0, length(fauni))

  c_p_fl <- matrix(nrow = patches, ncol = length(fleets), dimnames = list(NULL, fleet_names))
  r_p_fl <- matrix(nrow = patches, ncol = length(fleets), dimnames = list(NULL, fleet_names))
  prof_p_fl <- matrix(nrow = patches, ncol = length(fleets), dimnames = list(NULL, fleet_names))

  fishable <- rep(1, patches)

  # Per-fleet exponentially-smoothed history of the spatial allocation objective. NULL until the
  # first main-loop allocation populates it. Used only when a fleet has
  # memory_halflife > 0; see create_fleet() for the mechanism.
  # `_n` is the per-fleet count of post-seed update calls, used to ramp
  # the effective alpha during warm-up (Welford-style) so the smoothed surface
  # isn't anchored to the first observed objective.
  objective_history   <- vector("list", length(fleets))
  objective_history_n <- integer(length(fleets))
  names(objective_history)   <- fleet_names
  names(objective_history_n) <- fleet_names

  # Map of spatial_allocation -> buffet column name (mirrors allocate_effort()).
  alloc_to_mat <- c(
    revenue          = "r_p_fl",
    catch            = "c_p_fl",
    profit           = "prof_p_fl",
    rpue             = "rpue_p_fl",
    cpue             = "cpue_p_fl",
    ppue             = "ppue_p_fl",
    marginal_revenue = "mr_p_fl",
    marginal_profit  = "mp_p_fl"
  )

  # loop over steps
  for (s in 2:steps) {

    updated_e_p_f <- last_e_p_f

    year <- step_year[s]
    current_season <- step_season[s]

    # NOTE: starting_season is in signature but not used in your original script body here;
    # leaving behavior unchanged.

    if (length(manager$mpas) > 0) {
      manager$mpas$locations <- manager$mpas$locations |>
        dplyr::arrange(x, y)

      if (year == manager$mpas$mpa_year) {
        fishable <- manager$mpas$locations$mpa == 0
      }
    }

    fleet_fishable <- vector(mode = "list", length = length(fleet_names))
    fleet_concentrator <- vector(mode = "list", length = length(fleet_names))

    for (l in seq_along(fleet_names)) {
      tmp_fleet_fishable <- fishable

      if (!is.null(fleets[[l]]$fishing_grounds)) {
        tmp_fleet_fishable[fleets[[l]]$fishing_grounds$fishing_ground == FALSE] <- 0
      }

      fleet_fishable[[l]] <- tmp_fleet_fishable

      concentrator <- rep(1, patches)
      if (length(manager$mpas) > 0) {
        if (year >= manager$mpas$mpa_year && fleets[[l]]$mpa_response == "leave") {
          concentrator <- as.numeric(fishable)
        }
      }
      fleet_concentrator[[l]] <- concentrator
    }

    for (l in seq_along(fleet_names)) {
      # xx add ability to incorporate past revenues here. Idea. have a marker that lets you know if initial conditions were passed in. If s<= 2 or there were no initial conditions, pull this. If s<= 2 but there were initial conditions, pull the initial conditions for those steps

      if (fleets[[l]]$fleet_model == "manual") {
        total_effort <- fleets[[l]]$effort[s - 1]
      } else if (fleets[[l]]$fleet_model == "constant_effort") {
        # Use pre-quota effort so that quota reductions are transient:
        # the fleet "tries" its full intended effort each step, and only
        # gets reduced if the quota binds again this step.
        total_effort <- sum(pre_quota_e_p_f[,l] * fleet_concentrator[[l]])
      } else {
        # open_access / sole_owner: use realized effort; their dynamics
        # naturally allow effort to rebound via profitability signals
        total_effort <- sum(last_e_p_f[,l] * fleet_concentrator[[l]])
      }

      # -------------------------------------------------------------------------
      # OPEN ACCESS FLEET EFFORT DYNAMICS
      #
      # Concept:
      # Total fleet effort changes multiplicatively each step based on a
      # profitability signal derived from last step's realized revenue and cost.
      #
      # Key design goals:
      # - Works when revenue can be negative (e.g., species penalty / negative price)
      # - Scale-free (invariant to currency or cost scaling)
      # - Stable under stochastic recruitment and spatial heterogeneity
      # - Time-step invariant (same annual behavior regardless of step size)
      #
      # Structure:
      # 1) Compute normalized profitability signal in [-1, 1]
      # 2) Map signal -> behavioral response via tanh() (bounded adjustment speed)
      # 3) Convert annual responsiveness + caps to per-step equivalents using time_step
      # 4) Apply multiplicative effort update
      # -------------------------------------------------------------------------

      if (fleets[[l]]$fleet_model == "open_access" & s > 2) {

        # --- Effort cap (manager override, if provided) ---
        effort_cap_val <- Inf
        if (length(manager$effort_cap[[l]]) > 0) {
          effort_cap_val <- manager$effort_cap[[l]]
        }

        # --- Last step fleet-wide economic outcomes ---
        # NOTE: buffet contains outcomes from previous timestep
        last_revenue <- sum(buffet$r_p_fl[, names(fleets)[l]], na.rm = TRUE)

        last_cost <- last_revenue - sum(buffet$prof_p_fl[, names(fleets)[l]], na.rm = TRUE)

        # --- Sign-safe, scale-free profitability signal ---
        #
        # We use:
        #   signal = (R - C) / (|R| + C + eps)
        #
        # Why:
        # - bounded in [-1, 1]
        # - handles negative revenue safely
        # - invariant to rescaling economic units
        # - denominator represents total economic "throughput magnitude"
        #
        eps <- 1e-12 * (abs(last_revenue) + last_cost + 1)

        profit <- last_revenue - last_cost
        denom <- abs(last_revenue) + last_cost + eps

        oa_signal <- profit / denom

        # --- Convert annual parameters to per-step behavior ---
        #
        # time_step is fraction of a year per step
        #
        # Continuous-time interpretation:
        #   multiplier_step = exp( rho_year * time_step * response(signal) )
        #
        # So:
        #   rho_step = rho_year * time_step
        #
        # Caps also scale multiplicatively:
        #   m_step = m_year ^ time_step
        #

        rho_step <- fleets[[l]]$oa_rho_year * time_step
        k <- fleets[[l]]$oa_k

        m_max_step <- fleets[[l]]$oa_m_max_year^time_step
        m_min_step <- fleets[[l]]$oa_m_min_year^time_step

        # --- Behavioral response ---
        #
        # tanh() gives:
        # - linear response near break-even profitability
        # - saturated response under extreme profit or loss
        #
        # Economic interpretation:
        # - entry / exit speed is bounded by real-world frictions
        # - avoids unrealistic boom/bust effort dynamics
        #
        multiplier <- exp(rho_step * tanh(oa_signal / k))

        # --- Apply symmetric caps (per step) ---
        multiplier <- pmax(m_min_step, pmin(m_max_step, multiplier))

        # --- Apply effort update ---
        total_effort <- pmin(effort_cap_val, total_effort * multiplier)
      }

      # -------------------------------------------------------------------------
      # SOLE OWNER FLEET EFFORT DYNAMICS
      #
      # Identical to open access except the profitability signal is based on
      # the fleet-level marginal profit (effort-weighted mean of patch-level
      # marginal profits from calc_marginal_value) rather than average profit.
      #
      # Equilibrium: total effort where marginal profit = 0, i.e., MEY
      # (maximum economic yield), compared to open access which equilibrates
      # where average profit = 0.
      # -------------------------------------------------------------------------

      if (fleets[[l]]$fleet_model == "sole_owner" & s > 2) {

        # --- Effort cap (manager override, if provided) ---
        effort_cap_val <- Inf
        if (length(manager$effort_cap[[l]]) > 0) {
          effort_cap_val <- manager$effort_cap[[l]]
        }

        fl_name <- names(fleets)[l]

        # --- Fleet-level marginal profit signal ---
        # Effort-weighted mean of patch-level marginal profits:
        #   dΠ/dE = Σ_p (∂Π_p/∂e_p) * (e_p / E)
        # This represents the expected return of one additional unit of total
        # effort, distributed proportionally to the current spatial allocation.

        e_weights <- last_e_p_f[, l]
        e_total <- sum(e_weights)

        if (e_total > 0 && !is.null(buffet$mp_p_fl)) {

          mp_fleet <- sum(buffet$mp_p_fl[, fl_name] * e_weights, na.rm = TRUE) / e_total
          mr_fleet <- sum(buffet$mr_p_fl[, fl_name] * e_weights, na.rm = TRUE) / e_total

          # Marginal cost per unit effort (MR - MP)
          mc_fleet <- mr_fleet - mp_fleet

          # --- Sign-safe, scale-free marginal profitability signal ---
          # Mirrors the OA signal structure:
          #   OA:         signal = (R - C)  / (|R|  + C  + eps)   → 0 at avg profit = 0
          #   Sole owner: signal = MP / (|MR| + |MC| + eps)       → 0 at marginal profit = 0
          eps <- 1e-12 * (abs(mr_fleet) + abs(mc_fleet) + 1)
          oa_signal <- mp_fleet / (abs(mr_fleet) + abs(mc_fleet) + eps)

        } else {
          oa_signal <- 0  # no effort or marginals not yet available
        }

        # --- Convert annual parameters to per-step behavior ---
        rho_step <- fleets[[l]]$oa_rho_year * time_step
        k <- fleets[[l]]$oa_k

        m_max_step <- fleets[[l]]$oa_m_max_year^time_step
        m_min_step <- fleets[[l]]$oa_m_min_year^time_step

        # --- Behavioral response (same tanh machinery as open access) ---
        multiplier <- exp(rho_step * tanh(oa_signal / k))
        multiplier <- pmax(m_min_step, pmin(m_max_step, multiplier))

        # --- Apply effort update ---
        total_effort <- pmin(effort_cap_val, total_effort * multiplier)
      }


      ### allocate fleet in space ###
      if (s > 2){

        e_p <- last_e_p_f[,l]

        # Build a smoothed objective override when this fleet has memory enabled.
        halflife <- fleets[[l]]$memory_halflife
        if (is.null(halflife)) halflife <- 0
        alloc_type_l <- fleets[[l]]$spatial_allocation
        objective_override_l <- NULL
        if (halflife > 0 && alloc_type_l %in% names(alloc_to_mat)) {
          mat_name_l <- alloc_to_mat[[alloc_type_l]]
          obj_t      <- buffet[[mat_name_l]][, fleet_names[l]]
          open_l     <- as.logical(fleet_fishable[[l]])
          if (is.null(objective_history[[l]])) {
            # Seed: initialize with the first observed surface.
            smoothed <- obj_t
            objective_history_n[[l]] <- 1L
          } else {
            # Welford-style warm-up ramp: alpha_eff = max(alpha_target, 1/n) so
            # early steps behave like a running mean (no anchoring to the seed
            # surface) and only later relax to the asymptotic halflife. Once
            # 1/n drops below alpha_target, the ramp is done.
            objective_history_n[[l]] <- objective_history_n[[l]] + 1L
            n_obs          <- objective_history_n[[l]]
            # halflife is specified in years; convert to time steps so the
            # parameter is invariant to `seasons`.
            halflife_steps <- halflife * steps_per_year
            alpha_target   <- 1 - 0.5^(1 / halflife_steps)
            alpha_eff      <- max(alpha_target, 1 / n_obs)
            prev         <- objective_history[[l]]
            smoothed     <- prev
            # Update only on currently open patches; closed patches keep the
            # last-seen smoothed value ("freeze and resume").
            smoothed[open_l] <- alpha_eff * obj_t[open_l] + (1 - alpha_eff) * prev[open_l]
          }
          objective_history[[l]] <- smoothed
          objective_override_l <- matrix(
            smoothed, ncol = 1L,
            dimnames = list(NULL, fleet_names[l])
          )
        }

        current_effort <- allocate_effort(
          effort_by_patch = e_p,
          total_effort_by_fleet = total_effort,
          fleets = fleets[l],
          buffet = buffet,
          open_patch = fleet_fishable[[l]],
          flatness_tol = 1e-3,
          objective_override = objective_override_l
        )

        # warning("this is the bananas messy part. storage is indexed s-1 but fleet effort is indexed s. So, the effort sotred in s here is actually in storage in s-1, the effort that produced the outcomes in that time step of storage")

        updated_e_p_f[,l] <- current_effort$effort_new[,1]
      }

    } # close loop over fleets

    # Snapshot the intended effort before any quota reductions.
    # This is what fleets "wanted" to deploy this step.
    pre_quota_e_p_f_step <- updated_e_p_f

    for (f in seq_along(fauni)) {

      critter <- fauni[f]

      ages <- length(fauna[[critter]]$length_at_age)

      last_n_p_a <- storage[[s - 1]][[critter]]$n_p_a

      f_p_a <-
        matrix(0, nrow = patches, ncol = ages) # total fishing mortality by patch and age

      f_p_a_fl <-
        array(
          0,
          dim = c(patches, ages, length(fleets)),
          dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
        ) # storage for proportion of fishing mortality by patch, age, and fleet

      p_p_a_fl <-
        array(
          0,
          dim = c(patches, ages, length(fleets)),
          dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
        ) # storage for price by patch, age, and fleet

      for (l in seq_along(fleet_names)) {

        # tmp <-  outer(fleets[[l]]$metiers[[critter]]$spatial_catchability,  fleets[[l]]$metiers[[critter]]$sel_at_age, `*`)

        tmp <- fleets[[fleet_names[l]]]$metiers[[critter]]$vul_p_a

        f_p_a <-
          f_p_a + updated_e_p_f[,l] * tmp

        f_p_a_fl[, , l] <-
          updated_e_p_f[,l] * tmp

        p_p_a_fl[, , l] <- fleets[[l]]$metiers[[critter]]$price
      } # calculate cumulative f at age by patch

      f_p_a_fl <-
        f_p_a_fl / array(
          f_p_a,
          dim = c(patches, ages, length(fleets)),
          dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
        ) # f by patch, age, and fleet

      movement <- fauna[[critter]]$movement_matrix

      # if there is updated habitat for the critter in question in current time step, update habitat
      if ((length(habitat) > 0) &
          (critter %in% names(habitat))) {
        season_block <-
          which(sapply(fauna[[critter]]$movement_seasons, function(x, y) {
            any(y %in% x)
          }, x = current_season)) # figure out which season block you are in

        # update habitat in this time step
        current_habitat <- habitat[[critter]][[s - 1]]

        # need to use pivot_longer to match patch order from expand_grid
        # SPEED (spirit-of-comment): we precomputed the equivalent patch-ordered vector with as.vector()
        hab_vals <- habitat_vecs[[critter]][[s - 1]]

        neighbors <- find_neighbors(resolution = resolution, water = is.finite(hab_vals))

        edges <- Matrix::summary(neighbors)

        to <- edges$i

        from   <- edges$j

        delta_h <- hab_vals[to] - hab_vals[from]

        mean_area <- (patch_area_vec[to] + patch_area_vec[from]) / 2 # arithmetic mean of mean area

        mult <- exp((time_step * delta_h) / sqrt(mean_area))

        mult <- pmin(mult, fauna[[critter]]$max_hab_mult)

        P <- length(hab_vals)

        current_habitat <-  Matrix::sparseMatrix(
          i = to,
          j = from,
          x = mult,
          dims = c(P, P),
          dimnames = dimnames(neighbors)
        )

        diffusion_and_taxis <-
          fauna[[critter]]$diffusion_foundation[[season_block]] * current_habitat

        inst_movement_matrix <-
          prep_movement(diffusion_and_taxis, resolution = resolution)

        movement[[season_block]] <- sparsify_transition(as.matrix(expm::expm(as.matrix(inst_movement_matrix))))

        if (any(!is.finite(movement[[season_block]]))) {
          stop("scale of supplied habitat differences are too extreme, try rescaling so that the exponent of the differences are less extreme in magnitude")
        }
      }

      if (!(current_season %in% fishing_seasons[[critter]])) {
        f_p_a <- f_p_a * 0
      }

      fauna_rec_devs <- rep(exp(log_rec_devs[s, critter] - fauna[[critter]]$sigma_rec^2 / 2), patches)

      pop <- fauna[[critter]]$swim(
        season = current_season,
        adult_movement = movement,
        f_p_a = f_p_a,
        last_n_p_a = last_n_p_a,
        rec_devs = fauna_rec_devs
      )

      # tmp_e_p_fl <-
      #   purrr::list_cbind(unname(purrr::map(
      #     fleets, ~ data.frame(x = as.numeric(.x$e_p_s[, s]))
      #   )), name_repair = "unique_quiet")
      # colnames(tmp_e_p_fl) <- names(fleets)

      buffet <- allocate_yields(f_p_a_fl = f_p_a_fl, e_p_fl = updated_e_p_f,p_p_a_fl = p_p_a_fl, critter = critter,pop = pop, fauna = fauna, fleets = fleets, patches = patches, ages = ages )

      c_p_a_fl <- buffet$c_p_a_fl

      fmult <- 1
      if (length(manager$quotas[critter]) > 0) {
        if (manager$quotas[[critter]] < sum(c_p_a_fl, na.rm = TRUE)) {
          quota <- manager$quotas[[critter]]

          fmulter <- optim(
            par = 0.9,
            fn = marlin::quota_finder,
            quota = quota,
            fauna = fauna,
            current_season = current_season,
            movement = movement,
            f_p_a = f_p_a,
            last_n_p_a = last_n_p_a,
            f_p_a_fl = f_p_a_fl,
            critter = critter,
            patches = patches,
            ages = ages,
            fleets = fleets,
            rec_devs = fauna_rec_devs,
            lower = 0,
            upper = 1,
            method = "L-BFGS-B"
          )

          fmult <- fmulter$par
          f_p_a <- f_p_a * fmult

          pop <- fauna[[critter]]$swim(
            season = current_season,
            adult_movement = movement,
            f_p_a = f_p_a,
            last_n_p_a = last_n_p_a,
            rec_devs = fauna_rec_devs
          )
          # warning("quotas are broken right no until you update e and f by fleet post quota adaptation")
          # i think it's as simple as just downscaling by the same amount, since everything is a uniform scalar up and down?

          f_p_a_fl <- f_p_a_fl * fmult

          updated_e_p_f <- updated_e_p_f * fmult

          buffet <- allocate_yields(f_p_a_fl = f_p_a_fl, e_p_fl = updated_e_p_f,p_p_a_fl = p_p_a_fl, critter = critter,pop = pop, fauna = fauna, fleets = fleets, patches = patches, ages = ages )

        }
      }

      storage[[s - 1]][[critter]]$c_p_fl <-
        buffet$c_p_fl # store catch by patch  by fleet

      storage[[s - 1]][[critter]]$r_p_fl <-
        buffet$r_p_fl # store revenue by patch  by fleet

      storage[[s - 1]][[critter]]$prof_p_fl <-
        buffet$prof_p_fl # store profits by patch by fleet

      storage[[s - 1]][[critter]]$c_p_a_fl <-
        buffet$c_p_a_fl # catch stored in each model is the catch that came from the last time step, so put in the right place here

      storage[[s - 1]][[critter]]$r_p_a_fl <-
        buffet$r_p_a_fl # revenue stored in each model is the revenue that came from the last time step, so put in the right place here

      storage[[s - 1]][[critter]]$c_p_a <-
        pop$c_p_a # catch stored in each model is the catch that came from the last time step, so put in the right place here

      storage[[s - 1]][[critter]]$e_p_fl <-
        as.data.frame(updated_e_p_f) # store effort by patch by fleet (note that this is the same across species)

      storage[[s - 1]][[critter]]$f_p_a_fl <-
        f_p_a_fl # store effort by patch by fleet (note that this is the same across species)

      if (any(updated_e_p_f < 0)) {
        stop("something hase gone very wrong, effort is negative")
      }

      storage[[s]][[critter]] <- pop

      storage[[s]][[critter]]$b0 <- fauna[[critter]]$b0
    } # close fauni, much faster this way than dopar, who knew

    # gather per-critter yields from storage into a list
    yields_this_step <- setNames(
      lapply(fauni, function(cr) storage[[s - 1]][[cr]]),
      fauni
    )

    buffet <- aggregate_yields(yields_this_step, updated_e_p_f)

    last_e_p_f <- updated_e_p_f

    # Update pre_quota_e_p_f: for constant_effort fleets, preserve the
    # pre-quota intended effort so next step starts from full effort.
    # For dynamic fleets (OA/SO), use realized effort (they rebound naturally).
    for (l in seq_along(fleet_names)) {
      if (fleets[[l]]$fleet_model == "constant_effort") {
        pre_quota_e_p_f[, l] <- pre_quota_e_p_f_step[, l]
      } else {
        pre_quota_e_p_f[, l] <- updated_e_p_f[, l]
      }
    }

    # --- Marginal value signals for effort allocation -------------------------
    # Check if any fleet uses marginal-value allocation before paying the cost
    needs_marginals <- any(
      purrr::map_chr(fleets, "spatial_allocation") %in%
        c("marginal_revenue", "marginal_profit")
    ) || any(
      purrr::map_chr(fleets, "fleet_model") == "sole_owner"
    )

    if (needs_marginals) {
      # Gather current n_p_a from storage for go_fish
      marginal_n_p_a <- setNames(
        lapply(fauni, function(cr) storage[[s]][[cr]]$n_p_a),
        fauni
      )

      marginals <- calc_marginal_value(
        e_p_fl   = updated_e_p_f,
        fauna    = fauna,
        n_p_a    = marginal_n_p_a,
        fleets   = fleets,
        baseline = NULL,          # could pass go_fish result if available
        method   = "separable",   # fast default; "patch_loop" for small grids
        epsilon  = 1e-3
      )

      buffet$mr_p_fl <- marginals$mr_p_fl
      buffet$mp_p_fl <- marginals$mp_p_fl
    }

  } # close steps

  trimmed_names <- step_names[ifelse(keep_starting_step, 1, 2):steps]

  trimmed_names <-
    trimmed_names[1:(length(trimmed_names) - 1)] # a seriously annoying edge case in case you are only running one step

  storage <-
    storage[ifelse(keep_starting_step, 1, 2):pmax(2, steps - 1)] # since catch is retrospective, chop off last time step to ensure that every step has a catch history, and drop starting step is specified

  storage <-
    rlang::set_names(storage, nm = paste0(trimmed_names))

  storage <- purrr::map(storage, ~ rlang::set_names(.x, fauni))

  return(storage)
} # close function

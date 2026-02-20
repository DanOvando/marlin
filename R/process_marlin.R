#' Tidy Simulation Output from simmar
#'
#' @description
#' Converts the nested list output of \code{\link{simmar}} into two tidy
#' tibbles: one for population state (\code{$fauna}) and one for fleet
#' outcomes (\code{$fleets}). This is the primary post-processing step before
#' plotting or analysis.
#'
#' @details
#' ## Output columns
#'
#' **\code{$fauna}** — one row per critter × patch × age × step:
#' \itemize{
#'   \item \code{critter}: species name
#'   \item \code{patch}: patch index
#'   \item \code{x}, \code{y}: spatial coordinates
#'   \item \code{age}: age class (numeric)
#'   \item \code{mean_length}: mean length at age (from life history)
#'   \item \code{step}: decimal year (e.g. 5.75 = year 5, season 4 of 4)
#'   \item \code{n}: numbers at age
#'   \item \code{b}: biomass at age
#'   \item \code{ssb}: spawning stock biomass at age
#'   \item \code{c}: catch in numbers at age
#' }
#'
#' **\code{$fleets}** — one row per critter × patch × age × fleet × step:
#' \itemize{
#'   \item \code{critter}, \code{patch}, \code{x}, \code{y}, \code{age},
#'     \code{fleet}, \code{step}, \code{mean_length}: as above
#'   \item \code{catch}: catch in numbers
#'   \item \code{revenue}: revenue
#'   \item \code{effort}: effort units
#'   \item \code{cpue}: catch per unit effort (\code{NA} where effort = 0)
#' }
#'
#' When \code{keep_age = FALSE}, age classes are summed and a single
#' \code{age = "all"} row is returned per spatial unit per step.
#'
#' @param sim Named list returned by \code{\link{simmar}}. Step names should
#'   follow the \code{"year_season"} convention (with optional \code{"step_"}
#'   prefix, which is automatically stripped).
#' @param steps_to_keep Character or integer vector of step names / indices to
#'   include in the output. Useful for reducing memory use when only the final
#'   years are needed. Default \code{NULL} keeps all steps.
#' @param time_step Numeric. Fraction of a year per step (e.g. \code{0.25} for
#'   quarterly seasons). When \code{NULL}, inferred from the step names in
#'   \code{sim} (requires at least two steps).
#' @param keep_age Logical. If \code{TRUE} (default), return age-structured
#'   output. If \code{FALSE}, aggregate across all ages before returning.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{fauna}}{Tibble of population state; see Details.}
#'   \item{\code{fleets}}{Tibble of fleet outcomes; see Details.}
#' }
#'
#' @seealso \code{\link{simmar}}, \code{\link{plot_marlin}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- simmar(fauna = fauna, fleets = fleets, years = 50)
#'
#' # Full output (all steps, all ages)
#' proc <- process_marlin(sim, time_step = 1)
#'
#' # Last 10 years only, aggregated across ages
#' last_steps <- tail(names(sim), 10)
#' proc_last  <- process_marlin(sim,
#'                              steps_to_keep = last_steps,
#'                              time_step     = 1,
#'                              keep_age      = FALSE)
#'
#' # Plot SSB over time
#' plot_marlin(proc, plot_var = "ssb", plot_type = "time")
#' }
process_marlin <- function(sim,
                                steps_to_keep = NULL,
                                time_step = NULL,
                                keep_age = TRUE) {
  # requireNamespace("data.table", quietly = TRUE)
  # requireNamespace("stringi", quietly = TRUE)
  # requireNamespace("tibble", quietly = TRUE)

  dt <- data.table::as.data.table

  # ---- Step name cleanup ----
  nm0 <- names(sim)
  nm  <- stringi::stri_replace_all_fixed(nm0, "step_", "", vectorize_all = FALSE)
  names(sim) <- nm

  # ---- steps_to_keep default ----
  if (is.null(steps_to_keep)) steps_to_keep <- names(sim)
  sim <- sim[as.character(steps_to_keep)]

  # ---- infer seasons + time_step if needed ----
  seasons <- as.integer(stringi::stri_extract_last_regex(names(sim), "[0-9]+$"))
  seasons <- unique(seasons)

  if (is.null(time_step)) {
    if (length(sim) > 1) {
      time_step <- 1 / length(unique(seasons))
    } else {
      time_step <- 1
      warning(
        "unclear what time step length is; assuming it is 1. Consider setting manually with time_step parameter"
      )
    }
  }

  # ---- grid/coords (computed once) ----
  resolution <- sim[[1]][[1]]$resolution
  coords <- data.table::CJ(x = seq_len(resolution[1]),
                           y = seq_len(resolution[2]),
                           unique = TRUE)
  data.table::setorder(coords, x, y)
  coords[, patch := .I]
  data.table::setcolorder(coords, c("patch", "x", "y"))

  # ---- length lookup from sim[[1]] ----
  length_lookup <- data.table::rbindlist(
    lapply(names(sim[[1]]), function(critter) {
      obj <- sim[[1]][[critter]]
      dt(list(
        critter = critter,
        age = obj$ages,
        mean_length = obj$length_at_age
      ))
    }),
    use.names = TRUE
  )

  # ---- helper: tidy one critter's pop metrics into long, then wide ----
  tidy_pop_one_critter <- function(obj, critter_name) {
    ages <- obj$ages

    mats <- list(
      n_p_a   = obj$n_p_a,
      b_p_a   = obj$b_p_a,
      ssb_p_a = obj$ssb_p_a,
      c_p_a   = obj$c_p_a
    )

    out_long <- data.table::rbindlist(
      lapply(names(mats), function(metric_nm) {
        m <- mats[[metric_nm]]

        if (!keep_age) {
          v <- rowSums(m)
          dt_metric <- dt(list(patch = seq_along(v), age = "all", value = v))
        } else {
          dt_metric <- dt(m)
          dt_metric[, patch := .I]
          dt_metric <- data.table::melt(
            dt_metric,
            id.vars = "patch",
            variable.name = "age_idx",
            value.name = "value"
          )
          dt_metric[, age_i := as.integer(stringi::stri_replace_all_fixed(age_idx, "V", ""))]
          dt_metric[, age := ages[age_i]]
          dt_metric[, c("age_idx", "age_i") := NULL]
        }

        dt_metric[, metric := stringi::stri_replace_all_fixed(metric_nm, "_p_a", "", vectorize_all = FALSE)]
        dt_metric
      }),
      use.names = TRUE
    )

    out_long[, critter := critter_name]

    out_wide <- data.table::dcast(
      out_long,
      critter + patch + age ~ metric,
      value.var = "value"
    )

    out_wide <- out_wide[coords, on = "patch"]

    metric_cols <- intersect(c("n", "b", "ssb", "c"), names(out_wide))
    data.table::setcolorder(out_wide, c("critter", "patch", "x", "y", "age", metric_cols))
    out_wide
  }

  tidy_fauna_one_step <- function(step_list) {
    data.table::rbindlist(
      lapply(names(step_list), function(critter_name) {
        tidy_pop_one_critter(step_list[[critter_name]], critter_name)
      }),
      use.names = TRUE,
      fill = TRUE
    )
  }

  fauna_dt <- data.table::rbindlist(
    lapply(names(sim), function(step_nm) {
      tmp <- tidy_fauna_one_step(sim[[step_nm]])
      tmp[, step := step_nm]
      tmp
    }),
    use.names = TRUE,
    fill = TRUE
  )

  fauna_dt[, year   := as.integer(stringi::stri_extract_first_regex(step, "^[0-9]+"))]
  fauna_dt[, season := as.integer(stringi::stri_extract_last_regex(step, "[0-9]+$"))]
  fauna_dt[, step   := year + (season * time_step - time_step)]

  if (keep_age) {
    fauna_dt <- fauna_dt[length_lookup, on = .(critter, age)]
  } else {
    fauna_dt[, mean_length := as.numeric(NA)]
  }

  # ===========================
  # Fleets (UPDATED + faster)
  # ===========================
  tidy_fleets_one_critter <- function(obj, critter_name) {
    # Catch & revenue arrays have dimnames:
    #   [patch, age, fleet] with patch/age as character dimnames (e.g. "1","2",... and "0","1",...)
    catch_dt <- dt(as.data.table(as.table(obj$c_p_a_fl)))
    data.table::setnames(catch_dt, c("patch", "age", "fleet", "catch"))
    catch_dt[, `:=`(patch = as.numeric(patch), age = as.numeric(age), fleet = as.character(fleet))]

    rev_dt <- dt(as.data.table(as.table(obj$r_p_a_fl)))
    data.table::setnames(rev_dt, c("patch", "age", "fleet", "revenue"))
    rev_dt[, `:=`(patch = as.numeric(patch), age = as.numeric(age), fleet = as.character(fleet))]

    if (!keep_age) {
      catch_dt <- catch_dt[, .(catch = sum(catch, na.rm = TRUE)), by = .(patch, fleet)]
      catch_dt[, age := "all"]
      rev_dt <- rev_dt[, .(revenue = sum(revenue, na.rm = TRUE)), by = .(patch, fleet)]
      rev_dt[, age := "all"]
    }

    # effort: e_p_fl is a data.frame with columns = fleets, rows = patches
    eff_dt <- dt(obj$e_p_fl)
    eff_dt[, patch := .I]
    eff_dt <- data.table::melt(
      eff_dt,
      id.vars = "patch",
      variable.name = "fleet",
      value.name = "effort"
    )
    eff_dt[, fleet := as.character(fleet)]

    # add coords
    catch_dt <- catch_dt[coords, on = "patch"]
    rev_dt   <- rev_dt[coords, on = "patch"]

    # combine catch + revenue + effort
    out <- catch_dt[rev_dt, on = .(patch, age, fleet, x, y)]
    out <- out[eff_dt, on = .(patch, fleet)]

    out[, `:=`(
      critter = critter_name,
      cpue = catch / effort
    )]

    # column order similar to original
    preferred <- c("critter", "patch", "x", "y", "age", "fleet",
                   "catch", "revenue", "effort", "cpue")
    keep <- intersect(preferred, names(out))
    data.table::setcolorder(out, c(keep, setdiff(names(out), keep)))

    out
  }

  tidy_fleets_one_step <- function(step_list) {
    data.table::rbindlist(
      lapply(names(step_list), function(critter_name) {
        tidy_fleets_one_critter(step_list[[critter_name]], critter_name)
      }),
      use.names = TRUE,
      fill = TRUE
    )
  }

  fleets_dt <- data.table::rbindlist(
    lapply(names(sim), function(step_nm) {
      tmp <- tidy_fleets_one_step(sim[[step_nm]])
      tmp[, step := step_nm]
      tmp
    }),
    use.names = TRUE,
    fill = TRUE
  )

  fleets_dt[, year   := as.integer(stringi::stri_extract_first_regex(step, "^[0-9]+"))]
  fleets_dt[, season := as.integer(stringi::stri_extract_last_regex(step, "[0-9]+$"))]
  fleets_dt[, step   := year + (season * time_step - time_step)]

  if (keep_age) {
    # keep numeric age join
    fleets_dt <- fleets_dt[length_lookup, on = .(critter, age)]
  } else {
    fleets_dt[, mean_length := as.numeric(NA)]
  }

  list(
    fauna  = tibble::as_tibble(fauna_dt),
    fleets = tibble::as_tibble(fleets_dt)
  )
}

#' Plot Marlin
#' generates time series, spatial, or length composition plots from marlin run
#' @param ... named or unnamed outputs of process_marlin
#' @param steps_to_plot  which steps to include in the plot
#' @param plot_var the primary plotting variable, one of "c" (catch) "n" (numbers),"b" (biomass),"ssb" (spawning stock biomass)
#' @param plot_type the plot type, one of "time", "length", or "space"
#' @param fauna a fauna object corresponding to the critters in process_marlin outputs
#' @param drop_recruits TRUE to drop recruits from length composition plots
#'
#' @return a plot
#' @export
#'
#' @examples
#' \dontrun{
#'
#' plot_marlin(
#' `MPA: Sharks Offshore` = proc_offshore_mpa,
#' `No MPA` = proc_nearshore,
#' `MPA: Sharks Nearshore` = proc_nearshore_mpa,
# 'steps_to_plot = NA,
# 'plot_var = "b"
#' )
#'
#' }
#'
plot_marlin <- function(...,
                        steps_to_plot = NA,
                        plot_var = "ssb",
                        plot_type = "time",
                        fauna = NULL,
                        drop_recruits = TRUE,
                        plots = "fauna",
                        max_scale = TRUE) {
  fit_names <-
    names(list(...)) # allows users to pass and plot arbitrary numbers of objects from `process_marlin`
  if (is.null(fit_names)) {
    fit_names <- letters[seq_along(list(...))]

  }

  fit_names[fit_names == ''] <-
    letters[!letters %in% fit_names][seq_along(fit_names[fit_names == ''])] # add in names if only some are named for some reason

  fits <- list(...)  %>%
    purrr::set_names(fit_names)


  if (plots == "fauna") {
    fits <- purrr::map(fits, "fauna")

  }

  fit_frame <- dplyr::tibble(fit = fit_names,
                             temp = fits) %>%
    tidyr::unnest(cols = temp)
  if (all(is.na(steps_to_plot))) {
    steps_to_plot <- unique(fit_frame$step)
  }

  if (plot_type == "time") {
    out  <-  fit_frame %>%
      dplyr::group_by(step, critter, fit) %>%
      dplyr::summarise(across({
        {
          plot_var
        }
      }, ~ sum(., na.rm = TRUE))) %>%
      # ssb = sum(.data[[plot_var]])) %>%
      # dplyr::summarise(ssb = sum(.data[[plot_var]])) %>%
      dplyr::group_by(critter, fit) %>% {
        if (max_scale == TRUE) {
          dplyr::mutate(., dplyr::across({
            {
              plot_var
            }
          }, ~ (. / max(., na.rm = TRUE))))

        } else {
          .
        }

      } %>%
      dplyr::ungroup() %>%
      dplyr::filter(step %in% steps_to_plot) %>%
      ggplot(aes(step, .data[[plot_var]], color = fit)) +
      ggplot2::geom_hline(aes(yintercept = 0)) +
      ggplot2::geom_line(size = 3) +
      ggplot2::facet_wrap(~ critter, scales = "free_y") + {
        if (max_scale) {
          ggplot2::scale_y_continuous(
            limits = c(0, NA),
            labels = scales::percent,
            expand = ggplot2::expansion(mult = c(0, .1))
          )
        }
        else {
          ggplot2::scale_y_continuous(limits = c(0, NA),
                                      expand = ggplot2::expansion(mult = c(0, .1)))

        }
      } + ggplot2::scale_x_continuous(name = "Year") +
      ggplot2::scale_color_manual(
        name = "Fit",
        values = marlin::marlin_pal(palette = "diverging_fish")(dplyr::n_distinct(fit_frame$fit))
      ) +
      marlin::theme_marlin() +
      ggplot2::theme(legend.position = "top")

  } else if (plot_type == "length") {
    if (is.null(fauna)) {
      stop("plotting length compositions requires a supplied fauna object")
    }

    if (dplyr::n_distinct(steps_to_plot) > 10) {
      warning("trying to plot too many steps at once, cutting down to 10")

      steps_to_plot <-
        floor(seq(min(steps_to_plot), max(steps_to_plot), length.out = 10))

    }


    make_lcomps <- function(x, critter, plot_var = "n", fauna) {
      lkey <- fauna[[critter]]$length_at_age_key

      tallies <- x %>%
        dplyr::group_by(step, age) %>%
        dplyr::summarise(across({
          {
            plot_var
          }
        }, ~ sum(., na.rm = TRUE))) %>%
        dplyr::group_by(step) %>%
        tidyr::nest()


      age_to_length <- function(z, lkey, plot_var = "n") {
        thing_at_l <- z[[plot_var]] %*% as.matrix(lkey)

        out <-
          data.frame(length = as.numeric(colnames(lkey)),
                     thing = as.numeric(thing_at_l))

      }

      tallies <- tallies %>%
        dplyr::ungroup() %>%
        dplyr::mutate(tally = purrr::map(
          data,
          age_to_length,
          lkey = lkey,
          plot_var = plot_var
        )) %>%
        dplyr::select(step, tally) %>%
        tidyr::unnest(cols = tally)

    }

    tmp <- fit_frame %>%
      dplyr::ungroup() %>%
      dplyr::filter(step %in% steps_to_plot) %>%
      dplyr::group_by(critter, fit) %>%
      tidyr::nest() %>%
      dplyr::mutate(lcomps = purrr::map2(
        data,
        critter,
        make_lcomps,
        plot_var = plot_var,
        fauna = fauna
      )) %>%
      dplyr::select(critter, fit, lcomps) %>%
      tidyr::unnest(cols = lcomps)

    if (drop_recruits) {
      message("dropping recruits from plot since drop_recruits = TRUE")
      tmp <- tmp %>%
        dplyr::group_by(critter) %>%
        dplyr::mutate(lrank = percent_rank(length)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(lrank > 0.1)
    }

    out <- tmp %>%
      ggplot(aes(length, thing, fill = fit)) +
      ggplot2::geom_density(stat = "identity", color = "transparent") +
      ggplot2::facet_grid(step ~ critter, scales = "free_x") +
      ggplot2::scale_y_continuous(name = plot_var) +
      ggplot2::scale_x_continuous(name = "Length") +
      ggplot2::scale_fill_manual(
        name = "Fit",
        values = marlin::marlin_pal(palette = "diverging_fish")(dplyr::n_distinct(fit_frame$fit))
      ) +
      marlin::theme_marlin()

  } else if (plot_type == "age") {

    tmp <- fit_frame %>%
      dplyr::ungroup() %>%
      dplyr::filter(step %in% steps_to_plot)

    if (drop_recruits) {
      message("dropping recruits from plot since drop_recruits = TRUE")
      tmp <- tmp %>%
        dplyr::group_by(critter) %>%
        dplyr::filter(age != min(age))
    }

    tmp$thing <- tmp[plot_var]

    tmp <- tmp |>
      dplyr::filter(step == max(steps_to_plot)) |>
      dplyr::group_by(age, thing, step, critter,fit) |>
      dplyr::summarise(thing = sum(thing, na.rm = TRUE)) |>
      dplyr::ungroup()

    out <- tmp %>%
      ggplot(aes(age, thing, fill = fit)) +
      ggplot2::geom_density(stat = "identity", color = "transparent") +
      ggplot2::facet_grid(step ~ critter, scales = "free_x") +
      ggplot2::scale_y_continuous(name = plot_var) +
      ggplot2::scale_x_continuous(name = "Age") +
      ggplot2::scale_fill_manual(
        name = "Fit",
        values = marlin::marlin_pal(palette = "diverging_fish")(dplyr::n_distinct(fit_frame$fit))
      ) +
      marlin::theme_marlin()

  } else if (plot_type == "space") {
    if (dplyr::n_distinct(steps_to_plot) > 1) {
      warning(
        "Can only plot one time step for spatial plots, defaulting to last of the supplied steps"
      )
    }


    out <- fit_frame %>%
      dplyr::filter(step == max(steps_to_plot))


    if (max_scale) {
      out <- out %>%
        dplyr::group_by(x, y, critter, fit) %>%
        dplyr::summarise(across({
          {
            plot_var
          }
        }, ~ sum(., na.rm = TRUE))) %>%
        dplyr::group_by(critter) %>%
        dplyr::mutate(across({
          {
            plot_var
          }
        }, ~ (. / max(., na.rm = TRUE)))) %>%
        dplyr::ungroup()
    }

    out <- out %>%
      ggplot2::ggplot(aes(x, y, fill = round(.data[[plot_var]],2))) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(
        name = plot_var,
        guide = ggplot2::guide_colorbar(
          frame.colour = "black",
          tick.colour = "black",
          nbin = 1000,
          barwidth = unit(15, "lines")
        )
      ) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::facet_grid(critter ~ fit) +
      marlin::theme_marlin() +
      ggplot2::theme(legend.position = "top")



  }
  return(out)


}

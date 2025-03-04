#' R6 Class Representing a fishing fleet
#'
#' @description
#' A fleet object has all the required characteristics of a fishing fleet
#'
#' @details
#' creates fleet object with spaces for selectivity, mpa response, etc.

Metier <- R6::R6Class("metier",
  lock_objects = FALSE,
  public = list(
    #' @description
    #' Parameters for specific metier of a given fleet
    #'
    #' @param critter the name of the critter in the fauna object this metier applies to
    #' @param price the price per unit weight of the critter in question
    #' @param sel_form the selectivity form, one of "logistic", "dome","uniform", or "manual"
    #' @param sel_unit the unit of selectivity, one of "p_of_mat" which means selectivity is in proportion of age at maturity, or "length" where selectivity is in units of length
    #' @param sel_start the value of sel_unit at which selectivity "starts"
    #' @param sel_delta the delta parameter in the selectivity function
    #' @param catchability the catchability per uni effort pararmeter, generally overwritten by tune_fleet
    #' @param spatial_catchability a matrix of spatial q
    #' @param sel_at_age a manual vector of gear (contact) selectivity at age, where values are between 0 and 1
    #' @param sel05_anchor lower anchor at which (contact) selectivity is 0.05, either length or p_of_mat
    #' @param sel_at_linf contact selectivity at linf
    #' @param p_explt the proportion of total exploitation for a given critter coming from this metier. This value is relative to all other p_explt values for the critter in question. Set to 0 to have metier not catch critter at all
    initialize = function(critter = NA, # this might be redundant
                          price = 10,
                          sel_form = "logistic",
                          sel_unit = "p_of_mat",
                          sel_start = 1,
                          sel_delta = .1,
                          sel05_anchor = NULL,
                          sel_at_linf = NULL,
                          catchability = 0.2,
                          spatial_catchability = NA,
                          p_explt = 1,
                          sel_at_age = NULL) {
      catchability <- pmax(1e-9, catchability)

      self$price <- price

      self$sel_form <- sel_form

      self$sel_start <- sel_start

      self$sel_delta <- sel_delta

      self$catchability <- catchability

      self$p_explt <- p_explt

      self$port_distance <- NA

      self$sel_at_length <- NULL


      length_bins <-
        as.numeric(colnames(critter$length_at_age_key))

      ages <- length(critter$ages)

      if (sel_unit == "p_of_mat") {
        l_50_sel <-
          critter$length_50_mature * sel_start

        l_95_sel <-
          critter$length_50_mature * (sel_start + sel_delta)

        length_at_sel05 <- critter$length_50_mature * sel05_anchor
      } else if (sel_unit == "length") {
        l_50_sel <- sel_start

        l_95_sel <- sel_start + sel_delta

        length_at_sel05 <- sel05_anchor
      }

      if (sel_form == "logistic") {
        sel_at_bin <- ((1 / (1 + exp(
          -log(19) * ((length_bins - l_50_sel) / (l_95_sel - l_50_sel))
        ))))

        p_sel_at_age <-
          as.matrix(critter$length_at_age_key) %*% sel_at_bin

        sel_at_age <- p_sel_at_age

        self$sel_at_age <- as.numeric(sel_at_age)

        self$sel_at_length <- sel_at_bin
      } else if (sel_form == "dome") { # close logistic form if

        sel_at_bin <- dnorm(length_bins, l_50_sel, sd = (l_95_sel - l_50_sel))

        p_sel_at_age <-
          (as.matrix(critter$length_at_age_key) %*% sel_at_bin)

        sel_at_age <- p_sel_at_age / max(p_sel_at_age)

        self$sel_at_age <- as.numeric(sel_at_age)

        self$sel_at_length <- sel_at_bin / max(sel_at_bin)
      } else if (sel_form == "double_normal") {
        # Based on carruthers et al. 2014


        tune_double_normal <- function(log_sigmas,
                                       l_95_sel,
                                       ls,
                                       len_sel05 = 1,
                                       linf_sel = 0,
                                       output = "sigmas") {
          # sigma_asc = 0.2
          #
          # sigma_dsc = 10
          #
          # ls <- seq(0,100)
          #
          # len_sel05 = 1
          #
          # linf_sel = 0
          #
          # l_95_sel <- 42

          sigma_asc <- exp(log_sigmas[1])

          sigma_dsc <- exp(log_sigmas[2])

          asc <- dnorm(ls, l_95_sel, sigma_asc)

          asc <- asc / max(asc)

          dsc <- dnorm(ls, l_95_sel, sigma_dsc)

          dsc <- dsc / max(dsc)

          sel <- ifelse(ls <= l_95_sel, asc, dsc)

          out <- data.frame(length = ls, selectivity = sel)

          sel05_hat <- sel[which.min((ls - len_sel05)^2)]

          linf_sel_hat <- sel[length(sel)]

          if (output == "sigmas") {
            out <- (sel05_hat - 0.05)^2 + (linf_sel_hat - linf_sel)^2
          }

          return(out)
        }

        if (length_at_sel05 >= l_95_sel) {
          stop("length at sel05 must by less than length at peak selectivity")
        }
        tuned_sigmas <- optim(
          log(c(100, 100)),
          tune_double_normal,
          l_95_sel = l_95_sel,
          len_sel05 = length_at_sel05,
          linf_sel = sel_at_linf,
          ls = length_bins
        )

        sel_at_bin <- tune_double_normal(tuned_sigmas$par, l_95_sel = l_95_sel, len_sel05 = length_at_sel05, linf_sel = sel_at_linf, ls = length_bins, output = "sels")

        p_sel_at_age <-
          (as.matrix(critter$length_at_age_key) %*% sel_at_bin$selectivity)

        sel_at_age <- p_sel_at_age / max(p_sel_at_age)

        self$sel_at_age <- as.numeric(sel_at_age)

        self$sel_at_length <- sel_at_bin$selectivity
      } else if (sel_form == "manual") {
        if (is.null(sel_at_age)) {
          stop("sel_form = 'manual' but no manual set_at_age provided")
        }

        self$sel_at_age <- sel_at_age
      } else if (sel_form == "uniform") {
        self$sel_at_age <- rep(1, ages)

        self$sel_at_length <- rep(1, length(length_bins))
      } # close sel_form things

      if (all(is.na(spatial_catchability))) {
        self$spatial_catchability <- rep(catchability, critter$patches)
      } else {
        if (unique(dim(spatial_catchability)) != sqrt(critter$patches)) {
          stop(
            glue::glue(
              "spatial_catchability must either be NA or a {sqrt(critter$patches}) by {sqrt(critter$patches)} matrix"
            )
          )
        } # close dim check


        if (any(spatial_catchability < 0)) {
          spatial_catchability <- spatial_catchability - min(spatial_catchability)
        }

        tmp <- spatial_catchability %>%
          as.data.frame() %>%
          dplyr::mutate(x = 1:nrow(.)) %>%
          tidyr::pivot_longer(
            -x,
            names_to = "y",
            values_to = "catchability",
            names_prefix = "V"
          ) %>%
          dplyr::mutate(catchability = catchability / mean(catchability))

        self$spatial_catchability <- tmp$catchability * catchability
      } # close deal with spatial q
    } # close initialize
  ) # close public
) # close class

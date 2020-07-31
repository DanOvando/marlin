#' Process Marlin
#' 
#' tidy's arrays of outputs from marlin::simmar
#'
#' @param sim # the output from simmar
#' @param steps_to_keep # which steps you'd like to keep
#' @param time_step # the time step interval, as fractions of a year
#'
#' @return a tidy dataframe of population results
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' processed_marlin <- process_marlin(sim, steps_to_keep = 100, time_step = 0.25)
#' 
#' }
#' 
process_marlin <- function(sim,
                           steps_to_keep = NA,
                           time_step = 1) {
  if (is.na(steps_to_keep)) {
    steps_to_keep <-  seq_along(sim)
    
  }
  
  
  sim <-
    sim[steps_to_keep] # option to only select some years for memory's sake
  
  stepper <- function(x) {
    # x <- sim[[1]]
    
    
    tidy_marlin <- function(y, z) {
      pop <- y[c("n_p_a", "b_p_a", "ssb_p_a")]
      
      # create coordinates for each location
      tmp <-
        purrr::imap_dfr(pop,
                        ~ cbind(
                          patch = 1:nrow(.x),
                          as.data.frame(.x),
                          tidyr::expand_grid(x = 1:sqrt(nrow(.x)), y = 1:sqrt(nrow(.x)))
                        ),
                        .id = "metric") %>%
        tidyr::pivot_longer(
          # tidy ages
          tidyselect::contains("V"),
          names_to = "age",
          values_to = "value",
          names_prefix = "V",
          names_ptypes = list(value = integer())
        ) %>%
        dplyr::mutate(metric = gsub("_p_a", "", metric),
                      critter = z) %>%
        tidyr::pivot_wider(names_from = metric, values_from = value) %>%  # spread out metrics
        dplyr::select(critter, dplyr::everything())
      
      return(tmp)
      
    }
    
    out <- purrr::map2_dfr(x, names(x), tidy_marlin)
    
  }
  
  tidy_sim <-  purrr::imap_dfr(sim, ~ stepper(.x), .id = "step") %>%
    mutate(step = as.integer(step) * time_step)
  
  
}
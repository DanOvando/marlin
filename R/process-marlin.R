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
                           time_step = 1,
                           keep_age = TRUE) {
  if (all(is.na(steps_to_keep))) {
    steps_to_keep <-  seq_along(sim)
    
  }
  
  
  sim <-
    sim[steps_to_keep] # option to only select some years for memory's sake
  
  stepper <- function(x) {
    # x <- sim[[1]]
    
    
    tidy_marlin <- function(y, z) {
      ages <- y$ages
      
      pop <- y[c("n_p_a", "b_p_a", "ssb_p_a", "c_p_a")]
      
      # create coordinates for each location
      
      tmp <-
        purrr::imap_dfr(pop,
                        ~ cbind(
                          patch = 1:nrow(.x),
                          as.data.frame(.x),
                          tidyr::expand_grid(x = 1:sqrt(nrow(.x)), y = 1:sqrt(nrow(.x)))
                        ) %>% {
                          if (keep_age == FALSE) {
                            dplyr::rowwise(.) %>%
                              dplyr::mutate(tmp = sum(dplyr::c_across(tidyselect::contains("V")))) %>%
                              dplyr::select(-contains("V")) %>%
                              dplyr::rename(V0 = tmp)
                          } else {
                            .
                          }
                        },
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
        dplyr::select(critter, dplyr::everything()) %>% {
          if (keep_age == TRUE) {
            dplyr::mutate(., age =  rep(ages, n_distinct(.$patch)))
          } else {
            dplyr::mutate(., age = "all")
            
          }
        }
      
      return(tmp)
      
    }
    
    out <- purrr::map2_dfr(x, names(x), tidy_marlin)
    
  }
  
  
  tidy_sim <-  purrr::imap_dfr(sim, ~ stepper(.x), .id = "step") %>%
    dplyr::mutate(step = as.integer(step) * time_step,
           year = floor(as.integer(step) * time_step))
  
  # process fleets
  # ok this is the concept, but you need to wrap this in a map function to tidy by species
  
  
  fleet_stepper <- function(tmp) {
    get_fleet <- function(x, z) {
      #
      # tidy_catch <-  reshape::melt(sim[[1]]$bigeye$c_p_a_fl) %>%
      #   purrr::set_names("patch","age","fleet","catch")
      #
      tidy_catch <-  reshape::melt(x$c_p_a_fl) %>%
        purrr::set_names("patch", "age", "fleet", "catch")
      
      if (keep_age == FALSE){
        
      
        tidy_catch <-  tidy_catch %>% 
          dplyr::group_by(patch, fleet) %>% 
          dplyr::summarise(catch = sum(catch)) %>% 
          dplyr::mutate(age = "all")
        
        
      }
      
      coords <-
        tidyr::expand_grid(x = 1:sqrt(max(tidy_catch$patch)), y = 1:sqrt(max(tidy_catch$patch))) %>%
        dplyr::mutate(patch = 1:nrow(.))
      
      tidy_catch <- tidy_catch %>%
        left_join(coords, by = "patch")
      
      tidy_effort <- x$e_p_fl %>%
        purrr::set_names(paste0(colnames(.), "_effort")) %>%
        dplyr::mutate(patch = 1:nrow(.))
      
      tidy_fleet <- tidy_catch %>%
        dplyr::left_join(tidy_effort, by = "patch") %>%
        dplyr::mutate(critter = z)
      
    }
    
    # tmp <- sim[[1]]
    
    step_fleet <- purrr::imap_dfr(tmp, get_fleet)
    
    
  }
  
  tidy_sim_fleet <-
    purrr::imap_dfr(sim, ~ fleet_stepper(.x), .id = "step") %>%
    dplyr::mutate(step = as.integer(step) * time_step,
                  year = floor(as.integer(step) * time_step))
  
  out <- list(fauna = tidy_sim,
              fleets = tidy_sim_fleet)
  
  return(out)
  
}
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
                           time_step = NA,
                           keep_age = TRUE) {
  if (all(is.na(steps_to_keep))) {
    
    steps_to_keep <-  names(sim)
    
  }
  if (is.na(time_step)){
    if (length(sim) > 1){
      time_step <- as.numeric(names(sim))
      time_step <- time_step[2] - time_step[1]
    } else {
      time_step <-  1
      warning("unclear what time step length is; assuming it is 1. Consider setting manually with time_step parameter")
    }
  }
  
  sim <-
    sim[as.character(steps_to_keep)] # option to only select some years for memory's sake
  
  resolution <- sqrt(nrow(sim[[1]][[1]]$n_p_a))
  
  grid <- tidyr::expand_grid(x = 1:resolution, y = 1:resolution)
  
  stepper <- function(x, grid) {
    # x <- sim[[1]]
    
    
    tidy_marlin <- function(y, z, grid) {
      ages <- y$ages
      
      pop <- y[c("n_p_a", "b_p_a", "ssb_p_a", "c_p_a")]
      
      # create coordinates for each location
  
      bindfoo <- function(x, grid){
        
        x <- as.data.frame(x)
        
        if (keep_age == FALSE){
          
          x <- data.frame(V0 = rowSums(x))
          
        }
        
        out <- cbind(patch = 1:nrow(x),
                     x,
                     grid)
        
      }
      
      tmp <-
        purrr::map(pop,bindfoo,
                       grid= grid) %>%
        purrr::list_rbind(names_to = "metric") |> 
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
            dplyr::mutate(., age =  rep(ages, dplyr::n_distinct(.$patch)))
          } else {
            dplyr::mutate(., age = "all")
            
          }
        }
      
      return(tmp)
      
    }
    

    out <- purrr::map2(x, names(x), tidy_marlin, grid= grid) |> 
      purrr::list_rbind()
    
  }
  
  tidy_sim <-  purrr::imap(sim, ~ stepper(.x, grid = grid), .id = "step") %>%
    purrr::list_rbind(names_to = "step") |> 
    dplyr::mutate(step = as.numeric(step),
                  year = floor(as.numeric(step)))


  # process fleets
  # ok this is the concept, but you need to wrap this in a map function to tidy by species
  
  
  fleet_stepper <- function(tmp, grid) {
    get_fleet <- function(x, z, grid) {

      tidy_catch <-  data.frame(expand.grid(dimnames(x$c_p_a_fl)), value = as.vector(x$c_p_a_fl)) |> 
        purrr::set_names("patch", "age", "fleet", "catch") |> 
        dplyr::mutate(across(patch:age, ~as.numeric(as.character(.x))))
      
      tidy_rev <-  data.frame(expand.grid(dimnames(x$r_p_a_fl)), value = as.vector(x$r_p_a_fl)) |> 
        purrr::set_names("patch", "age", "fleet", "revenue") |> 
        dplyr::mutate(across(patch:age, ~as.numeric(as.character(.x))))
      
      
      if (keep_age == FALSE){
        
        tidy_catch <-  tidy_catch %>% 
          dplyr::group_by(patch, fleet) %>% 
          dplyr::summarise(catch = sum(catch, na.rm = TRUE)) %>% 
          dplyr::mutate(age = "all")
        
        tidy_rev <-  tidy_rev %>% 
          dplyr::group_by(patch, fleet) %>% 
          dplyr::summarise(revenue = sum(revenue, na.rm = TRUE)) %>% 
          dplyr::mutate(age = "all") 
        

      }
      
      coords <- grid %>%
        dplyr::mutate(patch = 1:nrow(.))
      
      tidy_catch <- tidy_catch %>%
        dplyr::left_join(coords, by = "patch")
      
      tidy_rev <- tidy_rev %>%
        dplyr::left_join(coords, by = "patch")
      
      
      tidy_effort <- x$e_p_fl %>%
        purrr::set_names(paste0(colnames(.), "_effort")) %>%
        dplyr::mutate(patch = 1:nrow(.))
      
      # hello? is it me you're looking for?
      tidy_fleet <- tidy_catch %>%
        dplyr::left_join(tidy_rev, by = c("patch", "age","fleet","x","y")) %>%
        dplyr::left_join(tidy_effort, by = "patch") %>%
        dplyr::mutate(critter = z)
      
    }
    
    # tmp <- sim[[1]]
    
    step_fleet <- purrr::imap(tmp, get_fleet, grid = grid) |> 
      purrr::list_rbind()
    
    
  }
  
  tidy_sim_fleet <-
    purrr::imap(sim, ~ fleet_stepper(.x, grid = grid)) %>%
    purrr::list_rbind(names_to = "step") |> 
    dplyr::mutate(step = as.numeric(step),
      year = floor(as.numeric(step)))
  

  out <- list(fauna = tidy_sim,
              fleets = tidy_sim_fleet)
  
  return(out)
  
}
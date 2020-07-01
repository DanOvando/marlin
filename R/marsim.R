marsim <- function(fauna = list(),
                   steps = 100,
                   n_cores = 1) {

  # asomething about movement and initial conditions isn't quite right. It's a problem with the
  # reslting movement matrix from create critter
# library(tidyverse)
#   resolution <- 25
#   habitat <- expand_grid(x =1:resolution, y = 1:resolution) %>%
#     mutate(habitat =  dnorm((x^2 + y^2), 600,100))
# 
#   habitat_mat <-
#     matrix(
#       rep(habitat$habitat, resolution),
#       nrow = resolution^2,
#       ncol = resolution^2,
#       byrow = TRUE
#     )
# 
#   skj_hab <- habitat_mat / rowSums(habitat_mat)
# 
#   habitat <- expand_grid(x =1:resolution, y = 1:resolution) %>%
#     mutate(habitat =  dnorm((x^2 + y^2), 100,2))
# 
#   habitat_mat <-
#     matrix(
#       rep(habitat$habitat, resolution),
#       nrow = resolution^2,
#       ncol = resolution^2,
#       byrow = TRUE
#     )
# 
#   bet_hab <- habitat_mat / rowSums(habitat_mat)
# 
# 
#   fauna <-
#     list("skipjack" = create_critter(
#       scientific_name = "Katsuwonus pelamis",
#       habitat = skj_hab,
#       adult_movement = 2
#     ),
#     "bigeye" = create_critter(
#       common_name = "bigeye tuna",
#       habitat = bet_hab,
#       adult_movement = 10
#     ))
# 
#   steps = 100
# 
#   n_cores = 1


  fauni <- names(fauna)

  patches <- unique(purrr::map_dbl(fauna, "patches"))

  initial_conditions <- purrr::map(fauna, c("unfished"))

  # a <- matrix(0.01, nrow = patches, ncol = length(fauna[[1]]$length_at_age))
  #
  # a[,1] <- fauna[[1]]$r0 / (patches)
  #
  # initial_conditions$`white seabass`$n_p_a <- a

  if (length(patches) > 1) {
    stop(
      "fauna have different habitat resolutions: set resolution to same number for all species!"
    )
  }

  # you'll need to fish down to joint EQ here


  # forecast from joint EQ, basically start MPA here


  storage <- vector("list", steps)

  storage[[1]] <- initial_conditions


  # a <- Sys.time()

  for (s in 2:steps) {
    for (f in seq_along(fauni)) {
      last_n_p_a <- storage[[s - 1]][[f]]$n_p_a

      # you can build a series of if statements here to sub in the correct species module

      pop <- marlin::sim_fish_pop(
        length_at_age = fauna[[f]]$length_at_age,
        weight_at_age = fauna[[f]]$weight_at_age,
        maturity_at_age = fauna[[f]]$maturity_at_age,
        steepness = fauna[[f]]$steepness,
        m = fauna[[f]]$m,
        patches = patches,
        burn_steps = 0,
        r0 = fauna[[f]]$r0,
        ssb0 = fauna[[f]]$ssb0,
        movement = fauna[[f]]$move_mat,
        last_n_p_a = last_n_p_a,
        tune_unfished = 0
      )

      storage[[s]][[f]] <- pop

    } # close fauni, much faster this way than dopar

  } #close steps

  # Sys.time() - a

  storage <- purrr::map(storage, ~ rlang::set_names(.x, fauni))

  # rec <- map(storage, ~.x[[2]]$n_p_a) %>%
  #   map_df(~tibble(rec = .x[,1]),.id = "i") %>%
  #   mutate(i = as.numeric(i)) %>%
  #   filter(i > 1) %>%
  #   group_by(i) %>%
  #   summarise(recs = sum(rec))
  # 
  # ssb <- map(storage, ~.x[[2]]$ssb_p_a)%>%
  #   map_df(~tibble(ssb = rowSums(.x)),.id = "i") %>%
  #   mutate(i = as.numeric(i)) %>%
  #   filter(i > 1) %>%
  #   group_by(i) %>%
  #   summarise(ssb = sum(ssb))
  # 
  # plot(ssb$ssb, rec$recs)
  # 
  # plot(ssb$ssb)
  # 
  # last(ssb$ssb) / fauna[[1]]$ssb0
  # 
  # ya <- rowSums(storage[[s]][[1]]$ssb_p_a)
  # 
  # check <- expand_grid(x = 1:sqrt(patches), y = 1:sqrt(patches)) %>%
  #   mutate(ssb = ya)
  # 
  #  ggplot(check, aes(x,y, fill = ssb)) +
  #   geom_tile() +
  #    scale_fill_viridis_c()




} # close function


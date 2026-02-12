go_fish <- function(e_p_f, fauna, fleets, patches, ages) {

# setup total f and f by fleet --------------------------------------------

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

  tmp <- fleets[[fleet_names[l]]]$metiers[[critter]]$vul_p_a

  f_p_a <-
    f_p_a + fleets[[l]]$e_p_s[, s] * tmp

  f_p_a_fl[, , l] <-
    fleets[[l]]$e_p_s[, s] * tmp

  p_p_a_fl[, , l] <- fleets[[l]]$metiers[[critter]]$price
} # calculate cumulative f at age by patch

f_p_a_fl <-
  f_p_a_fl / array(
    f_p_a,
    dim = c(patches, ages, length(fleets)),
    dimnames = list(1:patches, fauna[[critter]]$ages, names(fleets))
  ) # f by patch, age, and fleet


if (!(current_season %in% fishing_seasons[[critter]])) {
  f_p_a <- f_p_a * 0
}



# apply fishing mortality to population -----------------------------------


trial_run <- fauna[[critter]]$swim(
  season = current_season,
  adult_movement = movement,
  f_p_a = f_p_a,
  last_n_p_a = last_n_p_a,
  rec_devs = fauna_rec_devs
)

}



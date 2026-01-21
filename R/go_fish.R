go_fish <- function(effort,fauni, fleet){

# some weirdness. Fishing mortality happens post-movement at the movement.
# That means that there is no way at the moment for the fleet to know
# the actual biomass that will be available to them, since they can only observe the
# spatial biomass pre-movement.
# There are two solutions. The simplest is to test out the implications of simply moving the movement to post-fishing.
# Let's start there. If that works fine, then things are much simpler to implement here
# Whoops would also need to be pre-growth
# sigh that creates accounting problems though, since the catch by patch would be the catch that came from the biomass
# pre-movement, but the biomass and numbers etc. by patch would be post movement.
# but let's see the conequences. The basic question is do you allow fishers to see into
# the future or now. The argument for yes is that implicitly there is some information
# or experimental fishing going on that is too minimcal to model but helps them
# know where to go if distributions changed from one time step to another. It's
# where the meaningful effort setles out for the year
# The current approach also works, but means that allocations in space will be based on last year, which mi
#   might be subtopimal
#   Interesting, that doesn't seem to cause any problems, so might be the best
#   approach going forward since it makes perfect information possible, and you
#   can always index based on one more time step back to get the non-psychic
#   fishermen behavior
#   Too complicated for t-rex coding, need full screen and time to place this in context
#   of the IFD calculations.
#   Basically, you. need to account for spatial catchability and age structure inside the
#   IFD calculations loop you currently have




# calculate f  -----------------------------------------------------------

  for (f in seq_along(fauni)) {
    last_b_p_a <- storage[[s - 1]][[f]]$b_p_a

    last_e_p <- fleets[[l]]$e_p_s[, s - 1]

    # calculate fishable biomass in each patch for each species for that fleet

    # account for spatial catchability
    tmp <- 1 - exp(-(
      matrix(
        fleets[[l]]$metiers[[fauni[f]]]$spatial_catchability,
        nrow = nrow(last_b_p_a),
        ncol = ncol(last_b_p_a),
        byrow = FALSE
      ) *
        matrix(
          fleets[[l]]$metiers[[fauni[f]]]$sel_at_age,
          nrow = nrow(last_b_p_a),
          ncol = ncol(last_b_p_a),
          byrow = TRUE
        )
    ))

    last_b_p <-
      rowSums(last_b_p_a * tmp * (current_season %in% fishing_seasons[[f]])) * fleet_fishable

    r_p_f[, f] <-
      last_b_p * fleets[[l]]$metiers[[fauni[f]]]$price

    f_q[f] <- fleets[[l]]$metiers[[fauni[f]]]$catchability
  } # close fauni loop



# calculate catch / revenue (make toggle?) ---------------------------------------------------------



# return catch by patch across all critters -------------------------------




}

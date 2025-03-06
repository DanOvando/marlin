#' Prepare movement matrix
#'
#' @param multiplier multiplier for adjacency matrix
#' @param time_step time step in question
#' @param resolution spatial resolution
#'
#' @return a prepared movement matrix
#' @export
#'
prep_movement <-
  function(multiplier,
           resolution) {
    # reminder,
    #
    # .2^2 == exp(2 * log(.2)), hence notation in Thorson et al. going back and forth between different structures


    # set up spatial grid
    adjacent <-
      tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2]) %>%
      dist() %>%
      as.matrix()

    # Mark adjacent cells
    adjacent[adjacent != 1] <- 0

    # mark off water

    water <- !is.na(multiplier)

    multiplier[is.na(multiplier)] <- 0

    adjacent <- adjacent * water # set as non-adjacent any patches that contain land

    # generate base movement matrix
    move_mat <- adjacent * multiplier

    # fill in diagonal
    diag(move_mat) <- -1 * colSums(move_mat)

    # ensure matrix class
    move_mat <- as.matrix(move_mat)

    return(move_mat)
  } # close calc_move_mat

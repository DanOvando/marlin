#' Enforce Critter Names on Storage Steps
#'
#' Ensures that a storage step list is named and ordered consistently with the
#' fauna names used elsewhere in the simulation.
#'
#' @param x List representing one time-step of simulation storage.
#' @param fauni Character vector of expected fauna names.
#'
#' @return The input list \code{x}, named and reordered to match \code{fauni}.
#'
#' @keywords internal
enforce_critter_names <- function(x, fauni) {
  nm <- names(x)

  # if unnamed, we can only *assign* names (cannot verify mapping)
  if (is.null(nm) || anyNA(nm) || any(nm == "")) {
    if (length(x) != length(fauni)) {
      stop("storage step has wrong length; cannot assign critter names safely.")
    }
    names(x) <- fauni
    return(x)
  }

  # if named, verify exactly the same set and reorder
  if (!setequal(nm, fauni)) {
    stop("storage critter names don't match fauna names.")
  }
  x[fauni]
}

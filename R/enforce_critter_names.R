# ensure storage step is named/reordered to match fauna names
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

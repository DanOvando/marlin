bird <- R6::R6Class("Bird",
                public = list(
                  name = NULL,
                  hair = NULL,
                  initialize = function(name = NA, hair = NA) {
                    self$name <- name
                    self$hair <- hair
                    self$greet()
                  },
#' Title
#'
#' @param val 
#'
#' @return
#' @export
#'
#' @examples
                  set_hair = function(val) {
                    self$hair <- val
                  },
                  greet = function() {
                    cat(paste0("Hello, my name is ", self$name, ".\n"))
                  }
                )
)
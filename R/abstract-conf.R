#' Base class for tool settings
#'
#' Tool settings are checked when they are created.
#'
#' @keywords internal
AbstractConfig <- R6::R6Class(
  classname = "AbstractConfig",
  cloneable = FALSE,
  public = list(
    #' @description
    #' Store checked settings.
    #' @param values A named list of settings.
    initialize = function(values) {
      checkmate::assert_list(values, names = "unique")

      unknown <- setdiff(names(values), names(self))
      if (length(unknown) > 0L) {
        rlang::abort(sprintf(
          "Unknown configuration field(s): %s.",
          paste(unknown, collapse = ", ")
        ))
      }

      for (field in names(values)) {
        self[[field]] <- values[[field]]
      }
      private$.fields <- names(values)

      for (field in private$.fields) {
        lockBinding(field, self)
      }
    },

    #' @description
    #' Return a list of settings.
    #' @return A named list containing the current settings.
    as_list = function() {
      stats::setNames(
        lapply(private$.fields, function(field) self[[field]]),
        private$.fields
      )
    }
  ),
  private = list(
    .fields = character()
  )
)
